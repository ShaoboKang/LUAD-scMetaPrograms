#!/usr/bin/env python3
"""Batch AutoDock Vina v1.2.x for six GraphBAN target genes."""

import argparse
import os
import csv
import re
import shutil
import subprocess
import sys

from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


# ============================================================
# 1. Paths and target genes
# ============================================================

GENES = (
    "FOSL2",
    "ATF1",
    "CEBPG",
    "FOXM1",
    "KLF4",
    "TFDP1",
)

ROOT = Path(os.environ.get("LUAD_DOCKING_DIR", "~/LUAD_project/GraphBAN_teach/docking")).expanduser()

INPUT = ROOT / "input" / "GraphBAN_candidates.csv"
MANIFEST = ROOT / "input" / "receptors.tsv"

RESULT_DIR = ROOT / "results"
SUMMARY = RESULT_DIR / "vina_batch_summary.csv"
PREFLIGHT = RESULT_DIR / "vina_preflight.csv"

OUTPUT_COLUMNS = (
    "Gene_Name",
    "id",
    "SMILES",
    "bindingdb_pred",
    "biosnap_pred",
    "kiba_pred",
    "receptor_pdbqt",
    "box_config",
    "vina_affinity_kcal_mol",
    "status",
    "pose_file",
    "log_file",
)


# ============================================================
# 2. Create receptor manifest
# ============================================================

def create_manifest_if_absent():

    MANIFEST.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if MANIFEST.exists():
        print(
            f"Receptor manifest already exists: {MANIFEST}",
            flush=True
        )
        return

    with MANIFEST.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as fh:

        writer = csv.writer(
            fh,
            delimiter="\t"
        )

        writer.writerow((
            "Gene_Name",
            "receptor_pdbqt",
            "box_config"
        ))

        for gene in GENES:

            if gene == "FOXM1":

                writer.writerow((
                    gene,
                    "receptors_prepared/FOXM1_3G73.pdbqt",
                    "receptors_prepared/FOXM1_3G73_test.box.txt"
                ))

            else:
                writer.writerow((
                    gene,
                    "",
                    ""
                ))

    print(
        f"Created: {MANIFEST}",
        flush=True
    )

    print(
        "Other targets have blank receptor and box entries.",
        flush=True
    )


# ============================================================
# 3. Resolve file paths
# ============================================================

def resolve(value):

    if not value:
        return None

    path = Path(value).expanduser()

    if path.is_absolute():
        return path

    return ROOT / path


def valid_file(path):

    return (
        path is not None
        and path.is_file()
        and path.stat().st_size > 0
    )


def ligand_file(drug):

    return (
        ROOT
        / "ligands_pdbqt"
        / f"{drug}.pdbqt"
    )


# ============================================================
# 4. Read receptor manifest
# ============================================================

def read_manifest():

    mappings = {}

    with MANIFEST.open(
        newline="",
        encoding="utf-8-sig"
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        required = {
            "Gene_Name",
            "receptor_pdbqt",
            "box_config"
        }

        if not required.issubset(
            reader.fieldnames or []
        ):
            raise ValueError(
                f"Invalid manifest headers: {reader.fieldnames}"
            )

        for row in reader:

            gene = (
                row.get("Gene_Name") or ""
            ).strip().upper()

            if gene not in GENES:
                print(
                    f"WARNING: ignored unknown gene: {gene}"
                )
                continue

            if gene in mappings:
                raise ValueError(
                    f"Duplicate manifest gene: {gene}"
                )

            receptor = (
                row.get("receptor_pdbqt") or ""
            ).strip()

            box = (
                row.get("box_config") or ""
            ).strip()

            mappings[gene] = (
                resolve(receptor),
                resolve(box)
            )

    return mappings


# ============================================================
# 5. Read GraphBAN candidates
# ============================================================

def read_candidates():

    if not INPUT.is_file():
        raise FileNotFoundError(
            f"Cannot find GraphBAN candidate CSV: {INPUT}"
        )

    rows = []
    seen = {}

    with INPUT.open(
        newline="",
        encoding="utf-8-sig"
    ) as fh:

        reader = csv.DictReader(fh)

        required = {
            "Gene_Name",
            "id"
        }

        if not required.issubset(
            reader.fieldnames or []
        ):
            raise ValueError(
                f"CSV requires Gene_Name,id; "
                f"found: {reader.fieldnames}"
            )

        for row in reader:

            gene = (
                row.get("Gene_Name") or ""
            ).strip().upper()

            if gene not in GENES:
                continue

            drug = (
                row.get("id") or ""
            ).strip()

            if not re.fullmatch(
                r"[A-Za-z0-9_-]+",
                drug
            ):
                raise ValueError(
                    f"Invalid drug ID {drug!r} for {gene}"
                )

            row["Gene_Name"] = gene
            row["id"] = drug

            key = (gene, drug)

            smiles = (
                row.get("SMILES") or ""
            ).strip()

            if key in seen:

                if smiles != seen[key]:
                    raise ValueError(
                        f"Conflicting SMILES for "
                        f"{gene} and {drug}"
                    )

                continue

            seen[key] = smiles
            rows.append(row)

    return rows


# ============================================================
# 6. Check Vina docking box
# ============================================================

def check_box(box):

    required = {
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z"
    }

    entries = {}

    with box.open(
        encoding="utf-8"
    ) as fh:

        for line in fh:

            line = line.split(
                "#", 1
            )[0].strip()

            if "=" not in line:
                continue

            key, value = line.split(
                "=", 1
            )

            key = key.strip()

            if key in required:
                entries[key] = float(
                    value.strip()
                )

    if not required.issubset(entries):
        raise ValueError(
            f"Box missing fields "
            f"{required - entries.keys()}: {box}"
        )

    sizes = (
        entries["size_x"],
        entries["size_y"],
        entries["size_z"]
    )

    if min(sizes) <= 0:
        raise ValueError(
            f"Box dimensions must be positive: {box}"
        )

    volume = (
        entries["size_x"]
        * entries["size_y"]
        * entries["size_z"]
    )

    return volume


# ============================================================
# 7. Preflight checks
# ============================================================

def inspect(rows, mappings):

    counts = Counter(
        row["Gene_Name"]
        for row in rows
    )

    tasks = []
    report = []

    for gene in GENES:

        receptor, box = mappings.get(
            gene,
            (None, None)
        )

        pairs = [
            row for row in rows
            if row["Gene_Name"] == gene
        ]

        if not pairs:

            status = "no_candidate_pairs"

        elif receptor is None or box is None:

            status = "missing_receptor_or_box_mapping"

        elif not valid_file(receptor) or not valid_file(box):

            status = "missing_receptor_or_box_file"

        else:

            try:
                volume = check_box(box)

            except (ValueError, OSError) as exc:

                status = f"invalid_box: {exc}"

            else:

                if volume > 27000:

                    print(
                        f"WARNING: {gene} search box = "
                        f"{volume:.1f} A^3; exploratory only."
                    )

                status = "ready"

                for row in pairs:

                    ligand = ligand_file(
                        row["id"]
                    )

                    if valid_file(ligand):

                        tasks.append((
                            row,
                            receptor,
                            box,
                            ligand
                        ))

                    else:

                        status = "some_ligands_missing"

        ready_n = sum(
            1 for task in tasks
            if task[0]["Gene_Name"] == gene
        )

        report.append({
            "Gene_Name": gene,
            "candidate_pairs": counts[gene],
            "ready_pairs": ready_n,
            "status": status,
            "receptor_pdbqt": (
                "" if receptor is None
                else str(receptor)
            ),
            "box_config": (
                "" if box is None
                else str(box)
            )
        })

        print(
            f"{gene:7s} "
            f"candidates={counts[gene]:3d} "
            f"ready={ready_n:3d} "
            f"status={status}",
            flush=True
        )

    with PREFLIGHT.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as fh:

        columns = (
            "Gene_Name",
            "candidate_pairs",
            "ready_pairs",
            "status",
            "receptor_pdbqt",
            "box_config"
        )

        writer = csv.DictWriter(
            fh,
            fieldnames=columns
        )

        writer.writeheader()
        writer.writerows(report)

    return tasks, report


# ============================================================
# 8. Extract docking affinity
# ============================================================

def score_from_pose(pose):

    if not valid_file(pose):
        return None

    with pose.open(
        errors="replace"
    ) as fh:

        for line in fh:

            if line.startswith(
                "REMARK VINA RESULT:"
            ):

                value = (
                    line.split(":", 1)[1]
                    .split()[0]
                )

                return float(value)

    return None


# ============================================================
# 9. Execute one docking task
# ============================================================

def run_one(task, args):

    row, receptor, box, ligand = task

    gene = row["Gene_Name"]
    drug = row["id"]

    pose_dir = (
        ROOT
        / "poses"
        / "batch"
        / gene
    )

    log_dir = (
        ROOT
        / "logs"
        / "batch"
        / gene
    )

    pose_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    log_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    stem = (
        f"{gene}_{drug}_{receptor.stem}"
    )

    pose = pose_dir / f"{stem}.pdbqt"
    log = log_dir / f"{stem}.log"

    command = [
        "vina",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--config", str(box),
        "--exhaustiveness", str(args.exhaustiveness),
        "--cpu", str(args.cpu_per_job),
        "--num_modes", "9",
        "--seed", str(args.seed),
        "--out", str(pose)
    ]

    # Avoid reading an old pose from a previous run
    pose.unlink(missing_ok=True)

    with log.open(
        "w",
        encoding="utf-8"
    ) as fh:

        result = subprocess.run(
            command,
            stdout=fh,
            stderr=subprocess.STDOUT,
            check=False
        )

    if result.returncode == 0:

        score = score_from_pose(pose)

    else:

        score = None

    if result.returncode != 0:

        status = (
            f"vina_failed_exit_{result.returncode}"
        )

    elif score is None:

        status = "missing_vina_score"

    else:

        status = "success"

    data = {
        name: row.get(name, "")
        for name in OUTPUT_COLUMNS
    }

    data.update({
        "receptor_pdbqt": str(receptor),
        "box_config": str(box),
        "vina_affinity_kcal_mol": (
            "" if score is None else score
        ),
        "status": status,
        "pose_file": str(pose),
        "log_file": str(log)
    })

    return data


# ============================================================
# 10. Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Six-target GraphBAN Vina batch"
        )
    )

    parser.add_argument(
        "--run",
        action="store_true",
        help="Execute docking (default: preflight only)"
    )

    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Dock available pairs even if other targets are missing"
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help="Simultaneous docking processes"
    )

    parser.add_argument(
        "--cpu-per-job",
        type=int,
        default=2,
        help="Vina CPU threads per process"
    )

    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=8
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=2026
    )

    args = parser.parse_args()

    if min(
        args.workers,
        args.cpu_per_job,
        args.exhaustiveness
    ) < 1:

        parser.error(
            "workers, cpu-per-job, and "
            "exhaustiveness must be >= 1"
        )

    print(
        "Working directory:",
        ROOT,
        flush=True
    )

    RESULT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    create_manifest_if_absent()

    rows = read_candidates()
    mapping = read_manifest()

    print(
        f"Loaded {len(rows)} distinct drug-target pairs",
        flush=True
    )

    tasks, report = inspect(
        rows,
        mapping
    )

    print(
        "Preflight saved:",
        PREFLIGHT,
        flush=True
    )

    # Preflight only
    if not args.run:

        print(
            "CHECK ONLY: No Vina calculations launched. "
            "Add --run when ready.",
            flush=True
        )

        return

    if not shutil.which("vina"):

        raise RuntimeError(
            "Command 'vina' not found. "
            "Activate the docking environment."
        )

    incomplete = [
        row for row in report
        if row["status"] != "ready"
    ]

    if incomplete and not args.allow_partial:

        raise SystemExit(
            "Not all 7 targets are ready. "
            "Edit input/receptors.tsv and prepare ligands, "
            "or use --allow-partial to process ready pairs only."
        )

    if not tasks:

        raise SystemExit(
            "No ready docking jobs. "
            "Inspect vina_preflight.csv"
        )

    print(
        f"Launching {len(tasks)} docking jobs with "
        f"{args.workers} workers * "
        f"{args.cpu_per_job} CPU/job",
        flush=True
    )

    with SUMMARY.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=OUTPUT_COLUMNS
        )

        writer.writeheader()

        with ThreadPoolExecutor(
            max_workers=args.workers
        ) as pool:

            futures = {
                pool.submit(
                    run_one,
                    task,
                    args
                ): task

                for task in tasks
            }

            for future in as_completed(futures):

                try:

                    record = future.result()

                except Exception as exc:

                    row = futures[future][0]

                    record = {
                        k: row.get(k, "")
                        for k in OUTPUT_COLUMNS
                    }

                    record["status"] = (
                        f"python_error: {exc}"
                    )

                writer.writerow(record)
                fh.flush()

                print(
                    f"{record['Gene_Name']:7s} "
                    f"{record['id']:12s} "
                    f"{record['status']:24s} "
                    f"{record['vina_affinity_kcal_mol']}",
                    flush=True
                )

    print(
        "Done. Summary:",
        SUMMARY,
        flush=True
    )

    if incomplete:

        print(
            "PARTIAL RUN: Some targets or ligand files "
            "remain unprepared; see preflight CSV."
        )


# ============================================================
# 11. Entry point — do not delete
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except (
        FileNotFoundError,
        ValueError,
        RuntimeError
    ) as exc:

        print(
            f"ERROR: {exc}",
            file=sys.stderr,
            flush=True
        )

        sys.exit(1)
