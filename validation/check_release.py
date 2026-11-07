"""Check manuscript workflow paths, six-target screening boundaries and docking evidence."""
import ast
import csv
import importlib.util
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    files = list(ROOT.rglob('*.py'))
    for path in files:
        ast.parse(path.read_text(), filename=str(path))
    for path in ROOT.rglob('*.sh'):
        subprocess.run(['bash', '-n', str(path)], check=True)
    with (ROOT / 'workflow/manuscript_order.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert all((ROOT / row['script']).is_file() for row in rows)
    assert set(row['results_section'] for row in rows) == {'3.' + str(i) for i in range(1, 8)}
    prepare = load('prepare', 'scripts/11_compounds/01_prepare_screen.py')
    filterer = load('consensus', 'scripts/11_compounds/03_consensus_candidates.py')
    assert len(prepare.TFS) == 6 and 'VPS72' not in prepare.TFS
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fasta = tmp / 'proteins.fasta'; drugs = tmp / 'drugs.tsv'; out = tmp / 'input.csv'
        fasta.write_text(''.join('>{}\nMAAA\n'.format(g) for g in prepare.TFS))
        drugs.write_text('DB1\tCC\nDB2\tCCC\n')
        prepare.prepare(fasta, drugs, out)
        with out.open() as stream:
            assert len(list(csv.DictReader(stream))) == 12
        fasta.write_text(fasta.read_text() + '>VPS72\nMAAA\n')
        try:
            prepare.prepare(fasta, drugs, out)
        except ValueError:
            pass
        else:
            raise AssertionError('Unexpected VPS72 accepted')
        inputs = []
        for i in range(3):
            path = tmp / ('model{}.csv'.format(i)); inputs.append(path)
            with path.open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(['Gene_Name', 'id', 'SMILES', 'pred'])
                writer.writerows([
                    ['FOXM1', 'DB1', 'CC', .95],
                    ['FOXM1', 'DB2', 'CCC', .9 if i == 0 else .95],
                    ['VPS72', 'DB1', 'CC', .99]])
        selected = tmp / 'selected.csv'
        filterer.consensus(inputs, selected)
        with selected.open() as stream:
            result = list(csv.DictReader(stream))
        assert [(r['Gene_Name'], r['id']) for r in result] == [('FOXM1', 'DB1')]
        # Reject duplicate pairs instead of silently multiplying join results.
        with inputs[0].open('a') as stream:
            stream.write('FOXM1,DB1,CC,0.99\n')
        try:
            filterer.consensus(inputs, selected)
        except ValueError:
            pass
        else:
            raise AssertionError('Duplicate pair accepted')
    with (ROOT / 'evidence/manuscript_docking_pairs.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 8
    assert {r['Protein'] for r in rows} == set(prepare.TFS)
    assert len({r['Ligand'] for r in rows}) == 5
    assert all(float(r['Affinity (kcal/mol)']) < -4 for r in rows)
    print('PASS: {} Python files parsed; shell syntax valid; workflow paths valid.'.format(len(files)))
    print('PASS: six targets, strict >0.9 threshold, VPS72 exclusion and duplicate rejection.')
    print('PASS: eight docking pairs across six targets and five compounds, all below -4 kcal/mol.')


if __name__ == '__main__':
    main()
