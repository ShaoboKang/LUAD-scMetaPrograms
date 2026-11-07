#!/usr/bin/env bash
set -euo pipefail

D="${LUAD_DOCKING_DIR:-$HOME/LUAD_project/GraphBAN_teach/docking}"

mkdir -p \
    "$D/receptors_prepared" \
    "$D/logs" \
    "$D/results"

GENES=(FOSL2 ATF1 CEBPG KLF4 TFDP1)

echo "===== STEP 1: Prepare receptor PDBQT ====="

for gene in "${GENES[@]}"; do

    INPUT="$D/receptors_raw/${gene}_AlphaFold.pdb"
    OUTPUT="$D/receptors_prepared/${gene}_AlphaFold"
    LOG="$D/logs/${gene}_AlphaFold_prepare.log"

    echo ""
    echo "Processing ${gene} ..."

    if [ ! -s "$INPUT" ]; then
        echo "[FAILED] Missing PDB: $INPUT"
        continue
    fi

    rm -f "${OUTPUT}.pdbqt"

    if mk_prepare_receptor.py \
        --read_pdb "$INPUT" \
        -o "$OUTPUT" \
        -p > "$LOG" 2>&1 \
        && [ -s "${OUTPUT}.pdbqt" ]; then

        echo "[SUCCESS] ${gene} receptor"
        ls -lh "${OUTPUT}.pdbqt"

    else
        echo "[FAILED] ${gene} receptor"
        echo "See: $LOG"
        tail -n 12 "$LOG"
        rm -f "${OUTPUT}.pdbqt"
    fi

done


echo ""
echo "===== STEP 2: Generate exploratory docking boxes ====="

python - "$D" <<'PY'
from pathlib import Path
import csv
import sys

root = Path(sys.argv[1])

genes = [
    "FOSL2",
    "ATF1",
    "CEBPG",
    "KLF4",
    "TFDP1",
]

manifest = root / "input" / "receptors.tsv"
summary_file = root / "results" / "AF_domain_box_summary.csv"

if not manifest.exists():
    raise FileNotFoundError(
        "receptors.tsv missing. Run run_batch_vina.py first."
    )

with manifest.open(newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    fields = reader.fieldnames
    mappings = list(reader)

required = {"Gene_Name", "receptor_pdbqt", "box_config"}

if not required.issubset(fields or []):
    raise ValueError("Invalid receptors.tsv columns")

summary = []

for gene in genes:

    pdb = root / "receptors_raw" / f"{gene}_AlphaFold.pdb"

    rec = (
        root / "receptors_prepared"
        / f"{gene}_AlphaFold.pdbqt"
    )

    box = (
        root / "receptors_prepared"
        / f"{gene}_AlphaFold_test.box.txt"
    )

    if not rec.exists() or rec.stat().st_size == 0:
        print(f"[SKIP] {gene}: receptor not prepared")
        continue

    # Read pLDDT from CA atoms
    residues = {}

    with pdb.open() as f:
        for line in f:

            if not line.startswith("ATOM  "):
                continue

            if line[12:16].strip() != "CA":
                continue

            try:
                resnum = int(line[22:26])
                score = float(line[60:66])
                residues[resnum] = score
            except ValueError:
                continue

    # Find consecutive regions with pLDDT >=70
    runs = []
    current = []

    for resnum in sorted(residues):

        if residues[resnum] >= 70:

            if current and resnum != current[-1] + 1:
                runs.append(current)
                current = []

            current.append(resnum)

        else:
            if current:
                runs.append(current)
                current = []

    if current:
        runs.append(current)

    if not runs:
        print(f"[SKIP] {gene}: no confident region")
        continue

    longest = max(runs, key=len)

    start = longest[0]
    end = longest[-1]

    if len(longest) < 30:
        print(f"[SKIP] {gene}: confident region too short")
        continue

    # Collect atom coordinates from this region
    coordinates = []

    with pdb.open() as f:
        for line in f:

            if not line.startswith("ATOM  "):
                continue

            try:
                resnum = int(line[22:26])

                if not start <= resnum <= end:
                    continue

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                coordinates.append((x, y, z))

            except ValueError:
                continue

    if not coordinates:
        print(f"[SKIP] {gene}: no coordinates")
        continue

    center = []
    size = []

    for axis in range(3):

        values = [xyz[axis] for xyz in coordinates]

        lower = min(values)
        upper = max(values)

        center.append((lower + upper) / 2)

        # 10 Angstrom padding on each side
        size.append(
            max(24.0, upper - lower + 20.0)
        )

    cx, cy, cz = center
    sx, sy, sz = size

    volume = sx * sy * sz

    with box.open("w") as f:
        f.write(f"center_x = {cx:.3f}\n")
        f.write(f"center_y = {cy:.3f}\n")
        f.write(f"center_z = {cz:.3f}\n")
        f.write(f"size_x = {sx:.3f}\n")
        f.write(f"size_y = {sy:.3f}\n")
        f.write(f"size_z = {sz:.3f}\n")

    # Update existing receptor mapping
    for row in mappings:

        if row["Gene_Name"] != gene:
            continue

        if (
            not row.get("receptor_pdbqt", "").strip()
            and not row.get("box_config", "").strip()
        ):
            row["receptor_pdbqt"] = (
                f"receptors_prepared/{gene}_AlphaFold.pdbqt"
            )

            row["box_config"] = (
                f"receptors_prepared/{gene}_AlphaFold_test.box.txt"
            )

    summary.append({
        "Gene_Name": gene,
        "Domain_start": start,
        "Domain_end": end,
        "Domain_residues": len(longest),
        "center_x": round(cx, 3),
        "center_y": round(cy, 3),
        "center_z": round(cz, 3),
        "size_x": round(sx, 3),
        "size_y": round(sy, 3),
        "size_z": round(sz, 3),
        "volume_A3": round(volume, 1),
        "receptor": str(rec),
        "box": str(box),
        "validation": "exploratory_unvalidated"
    })

    print(
        f"[OK] {gene}: region {start}-{end}; "
        f"box volume {volume:.1f} A^3"
    )

    if volume > 27000:
        print(
            f"[WARNING] {gene}: large search box; "
            "docking may require more sampling."
        )

# Save updated receptor manifest
with manifest.open("w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(mappings)

# Save box summary
if summary:

    with summary_file.open("w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=summary[0].keys()
        )

        writer.writeheader()
        writer.writerows(summary)

print("\nSaved:", manifest)
print("Saved:", summary_file)
print("Prepared:", len(summary), "of", len(genes))

PY

echo ""
echo "===== STEP 3: Check generated files ====="

ls -lh "$D/receptors_prepared/"*_AlphaFold.pdbqt \
       "$D/receptors_prepared/"*_AlphaFold_test.box.txt

echo ""
echo "Preparation completed."
