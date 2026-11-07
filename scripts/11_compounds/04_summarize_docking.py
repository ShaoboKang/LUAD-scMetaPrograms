"""Select manuscript-TF interactions with docking energies strictly below -4 kcal/mol.
Usage: python3 .../04_summarize_docking.py Table_2_2_molecular_docking.csv output.csv
"""
import argparse
import csv
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec
spec = spec_from_file_location('prepare_screen', Path(__file__).with_name('01_prepare_screen.py'))
module = module_from_spec(spec); spec.loader.exec_module(module)


def summarize(source, output):
    with open(source, encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        rows = [r for r in reader if r['Protein'] in module.TFS and float(r['Affinity (kcal/mol)']) < -4.0]
    with open(output, 'w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    print('Retained {} interactions across {} signature TFs'.format(len(rows), len({r['Protein'] for r in rows})))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('source'); p.add_argument('output')
    a = p.parse_args(); summarize(a.source, a.output)
