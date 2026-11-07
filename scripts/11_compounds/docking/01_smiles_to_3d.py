import sys
import re
import pandas as pd

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

input_csv = Path(sys.argv[1])
outdir = input_csv.parent.parent / "ligands_sdf"
outdir.mkdir(parents=True, exist_ok=True)
(input_csv.parent.parent / "results").mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_csv, dtype={"id": str})

required = {"Gene_Name", "id", "SMILES"}
assert required.issubset(df.columns), (
    f"Missing columns: {required - set(df.columns)}"
)

# Reject inconsistent SMILES for a compound ID.
n_smiles = df.groupby("id")["SMILES"].nunique()

if (n_smiles > 1).any():
    raise ValueError(
        "Some drug IDs have multiple SMILES. "
        "Please check them before preparation."
    )

# Generate one conformer per unique compound.
drugs = df[["id", "SMILES"]].drop_duplicates()

records = []

for _, row in drugs.iterrows():

    drug_id = str(row["id"])
    smiles = str(row["SMILES"])
    status = "success"

    try:
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            raise ValueError("Invalid SMILES")

        if "." in smiles:
            raise ValueError(
                "Multiple fragments or salts: manual review needed"
            )

        mw = Descriptors.MolWt(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)

        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        params.randomSeed = 2026

        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError("3D embedding failed")

        if AllChem.MMFFHasAllMoleculeParams(mol):
            result = AllChem.MMFFOptimizeMolecule(
                mol, maxIters=500
            )
        elif AllChem.UFFHasAllMoleculeParams(mol):
            result = AllChem.UFFOptimizeMolecule(
                mol, maxIters=500
            )
        else:
            raise ValueError("No suitable force field parameters")

        if result != 0:
            status = "optimization_not_converged"

        safe_id = re.sub(r"[^A-Za-z0-9_-]", "_", drug_id)
        outfile = outdir / f"{safe_id}.sdf"

        mol.SetProp("_Name", drug_id)
        mol.SetProp("DrugID", drug_id)
        mol.SetProp("SMILES", smiles)

        writer = Chem.SDWriter(str(outfile))
        writer.write(mol)
        writer.close()

        records.append({
            "id": drug_id,
            "MW": mw,
            "RotatableBonds": rotb,
            "status": status,
            "file": str(outfile)
        })

    except Exception as e:
        records.append({
            "id": drug_id,
            "MW": None,
            "RotatableBonds": None,
            "status": str(e),
            "file": ""
        })

    print(drug_id, records[-1]["status"], flush=True)

summary = pd.DataFrame(records)

summary.to_csv(
    input_csv.parent.parent / "results" / "ligand_preparation.csv",
    index=False
)

print("\nLigand preparation finished")
print(summary["status"].value_counts())
