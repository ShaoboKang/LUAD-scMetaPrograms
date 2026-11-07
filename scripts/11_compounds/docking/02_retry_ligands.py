import os
import math
import pandas as pd

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

BASE = Path(os.environ.get("LUAD_DOCKING_DIR", "~/LUAD_project/GraphBAN_teach/docking")).expanduser()

INPUT = BASE / "input/GraphBAN_candidates.csv"
OUTDIR = BASE / "ligands_sdf"
BACKUP = OUTDIR / "unconverged_original"
SUMMARY = BASE / "results/ligand_retry_summary.csv"

OUTDIR.mkdir(parents=True, exist_ok=True)
BACKUP.mkdir(parents=True, exist_ok=True)
SUMMARY.parent.mkdir(parents=True, exist_ok=True)

# Retry the three compounds that failed initial optimization.
TARGETS = ["DB09372", "DB06206", "DB15008"]

df = pd.read_csv(INPUT, dtype={"id": str})

results = []

for drug_id in TARGETS:

    print(f"\n========== {drug_id} ==========", flush=True)

    matches = (
        df.loc[df["id"] == drug_id, "SMILES"]
        .dropna()
        .drop_duplicates()
    )

    if len(matches) != 1:
        print("ERROR: Missing or inconsistent SMILES")
        results.append({
            "id": drug_id,
            "status": "invalid_or_inconsistent_input"
        })
        continue

    smiles = matches.iloc[0]
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("ERROR: Invalid SMILES")
        results.append({
            "id": drug_id,
            "status": "invalid_smiles"
        })
        continue

    # Check molecular properties.
    mw = Descriptors.MolWt(mol)
    heavy = mol.GetNumHeavyAtoms()
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)

    print(f"Molecular weight: {mw:.2f}")
    print(f"Heavy atoms: {heavy}")
    print(f"Rotatable bonds: {rotb}")

    # Archive the unconverged DB06206 conformer before retrying.
    outfile = OUTDIR / f"{drug_id}.sdf"

    if drug_id == "DB06206" and outfile.exists():
        backup_file = BACKUP / f"{drug_id}_previous.sdf"
        outfile.replace(backup_file)
        print("Previous unvalidated SDF moved to backup")

    best_mol = None
    best_energy = float("inf")
    best_method = None

    # Try multiple random seeds and coordinate strategies.
    strategies = [
        (2026, False),
        (12345, False),
        (2026, True),
        (12345, True),
    ]

    for seed, random_coords in strategies:

        print(
            f"Trying seed={seed}, "
            f"useRandomCoords={random_coords}",
            flush=True
        )

        m = Chem.AddHs(Chem.Mol(mol))

        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.maxIterations = 500
        params.useRandomCoords = random_coords
        params.useMacrocycleTorsions = True

        try:
            status = AllChem.EmbedMolecule(m, params)

            if status != 0:
                print("  Embedding failed")
                continue

            # Prefer MMFF94 when parameters are available.
            if AllChem.MMFFHasAllMoleculeParams(m):

                method = "MMFF94"

                opt = AllChem.MMFFOptimizeMolecule(
                    m, maxIters=2000
                )

                if opt == 1:
                    opt = AllChem.MMFFOptimizeMolecule(
                        m, maxIters=2000
                    )

                if opt != 0:
                    print("  MMFF optimization not converged")
                    continue

                props = AllChem.MMFFGetMoleculeProperties(m)
                ff = AllChem.MMFFGetMoleculeForceField(
                    m, props
                )

            elif AllChem.UFFHasAllMoleculeParams(m):

                method = "UFF"

                opt = AllChem.UFFOptimizeMolecule(
                    m, maxIters=2000
                )

                if opt == 1:
                    opt = AllChem.UFFOptimizeMolecule(
                        m, maxIters=2000
                    )

                if opt != 0:
                    print("  UFF optimization not converged")
                    continue

                ff = AllChem.UFFGetMoleculeForceField(m)

            else:
                print("  No force field parameters")
                continue

            energy = ff.CalcEnergy()

            if not math.isfinite(energy):
                print("  Invalid force-field energy")
                continue

            print(
                f"  Success! {method} energy={energy:.3f}",
                flush=True
            )

            if energy < best_energy:
                best_energy = energy
                best_mol = Chem.Mol(m)
                best_method = method

        except Exception as e:
            print(f"  Error: {e}")

    if best_mol is not None:

        best_mol.SetProp("_Name", drug_id)
        best_mol.SetProp("DrugID", drug_id)
        best_mol.SetProp("SMILES", smiles)

        writer = Chem.SDWriter(str(outfile))
        writer.write(best_mol)
        writer.close()

        final_status = "success"

        print(f"Saved: {outfile}")

    else:

        final_status = "failed_after_retry"
        best_energy = None
        best_method = None

        print("No converged conformer found")

    results.append({
        "id": drug_id,
        "MW": mw,
        "HeavyAtoms": heavy,
        "RotatableBonds": rotb,
        "status": final_status,
        "forcefield": best_method,
        "energy": best_energy,
        "file": str(outfile) if best_mol is not None else ""
    })

pd.DataFrame(results).to_csv(SUMMARY, index=False)

print("\n========== FINAL SUMMARY ==========")
print(pd.DataFrame(results).to_string(index=False))
print(f"\nSummary saved: {SUMMARY}")
