"""Prepare manuscript six-target x DrugBank inputs from local FASTA and TSV files.
Usage: python3 scripts/11_compounds/01_prepare_screen.py proteins.fasta drugs.tsv output.csv
Drug table has two tab-separated columns without header: DrugBank ID, SMILES.
No sequence or licensed compound library is downloaded or redistributed.
"""
import argparse
import csv
import re
from pathlib import Path

TFS = ('FOXM1', 'TFDP1', 'FOSL2', 'KLF4', 'CEBPG', 'ATF1')


def prepare(fasta, drugs, output):
    """Require the six screening targets confirmed by the author; exclude VPS72."""
    proteins = {}
    name = None
    for line in Path(fasta).read_text().splitlines():
        if line.startswith('>'):
            match = re.search(r'\bGN=(\w+)', line)
            name = match.group(1) if match else line[1:].split()[0]
            if name in proteins:
                raise ValueError('Duplicate protein: ' + name)
            proteins[name] = ''
        elif line.strip():
            if name is None:
                raise ValueError('FASTA sequence before header')
            proteins[name] += line.strip()
    if set(proteins) != set(TFS) or any(not proteins[g] for g in TFS):
        raise ValueError('FASTA must contain exactly the six screening TFs: ' + ', '.join(TFS))
    with open(drugs, newline='') as stream:
        compounds = list(csv.reader(stream, delimiter='\t'))
    if not compounds or any(len(row) != 2 or not all(row) for row in compounds):
        raise ValueError('Expected nonempty ID/SMILES TSV without header')
    if len({row[0] for row in compounds}) != len(compounds):
        raise ValueError('Duplicate compound IDs')
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['Gene_Name', 'id', 'SMILES', 'Protein', 'Y'])
        for gene in TFS:
            for compound, smiles in compounds:
                writer.writerow([gene, compound, smiles, proteins[gene], 1])
    print('Prepared {} TFs x {} compounds = {} pairs'.format(len(TFS), len(compounds), len(TFS)*len(compounds)))
    if len(compounds) != 11167:
        print('Note: manuscript library size is 11,167; supplied library differs.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fasta'); parser.add_argument('drugs'); parser.add_argument('output')
    args = parser.parse_args()
    prepare(args.fasta, args.drugs, args.output)
