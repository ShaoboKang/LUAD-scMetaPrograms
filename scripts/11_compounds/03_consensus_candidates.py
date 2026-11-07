"""Intersect three GraphBAN prediction tables using strict score > 0.9.
Usage: python3 .../03_consensus_candidates.py biosnap.csv kiba.csv bindingdb.csv out.csv
All historical targets remain in the audit; screening output uses the six author-confirmed TFs.
"""
import argparse
import csv
import math
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec

spec = spec_from_file_location('prepare_screen', Path(__file__).with_name('01_prepare_screen.py'))
module = module_from_spec(spec); spec.loader.exec_module(module)
TFS = set(module.TFS)


def load_predictions(path):
    result = {}
    with open(path, newline='') as stream:
        reader = csv.DictReader(stream)
        if not {'Gene_Name', 'id', 'SMILES', 'pred'}.issubset(reader.fieldnames or []):
            raise ValueError('Missing prediction columns: ' + str(path))
        for row in reader:
            key = (row['Gene_Name'], row['id'])
            if key in result:
                raise ValueError('Duplicate TF-compound pair: ' + str(key))
            score = float(row['pred'])
            if not math.isfinite(score):
                raise ValueError('Non-finite prediction: ' + str(key))
            result[key] = (row['SMILES'], score)
    return result


def consensus(paths, output):
    tables = [load_predictions(p) for p in paths]
    if any(set(t) != set(tables[0]) for t in tables[1:]):
        raise ValueError('Three models must predict the same input pairs')
    rows = []
    for key in sorted(tables[0]):
        values = [t[key] for t in tables]
        if len({v[0] for v in values}) != 1:
            raise ValueError('SMILES disagree across models: ' + str(key))
        if all(v[1] > .9 for v in values):
            rows.append(dict(Gene_Name=key[0], id=key[1], SMILES=values[0][0],
                biosnap_pred=values[0][1], kiba_pred=values[1][1], bindingdb_pred=values[2][1]))
    output = Path(output); output.parent.mkdir(parents=True, exist_ok=True)
    fields = ['Gene_Name', 'id', 'SMILES', 'biosnap_pred', 'kiba_pred', 'bindingdb_pred']
    for dest, subset in [(output, [r for r in rows if r['Gene_Name'] in TFS]),
                         (output.with_name(output.stem + '_all_targets.csv'), rows)]:
        with dest.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(subset)
    missing = TFS - {key[0] for key in tables[0]}
    if missing:
        print('Historical screen did not evaluate: ' + ', '.join(sorted(missing)))
    print('Saved {} manuscript-TF consensus pairs; {} all-target pairs'.format(
        sum(r['Gene_Name'] in TFS for r in rows), len(rows)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['biosnap', 'kiba', 'bindingdb', 'output']: p.add_argument(name)
    a = p.parse_args(); consensus([a.biosnap, a.kiba, a.bindingdb], a.output)
