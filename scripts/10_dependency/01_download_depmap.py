"""Download a pinned official DepMap release and verify every file's MD5.
Run before 09_FOXM1_DepMap_LUAD_CRISPR.R. Uses only Python's standard library
and curl. Full source matrix is retained for reproducibility, not portal scraping.
"""
from pathlib import Path
import os
import hashlib
import json
import subprocess
from datetime import datetime, timezone

ROOT = Path(os.environ.get("LUAD_ROOT", "~/LUAD_project")).expanduser()
OUT = ROOT / 'Result_8/results_data/FOXM1_DepMap/raw'
OUT.mkdir(parents=True, exist_ok=True)
API = 'https://api.figshare.com/v2/articles/27993248/versions/1'


def download(url, path):
    partial = path.with_name(path.name + '.partial')
    subprocess.run(['curl', '--location', '--fail', '--retry', '3',
                    '--connect-timeout', '30', '--silent', '--show-error',
                    url, '--output', str(partial)], check=True)
    partial.replace(path)


def md5(path):
    digest = hashlib.md5()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


metadata_path = OUT / 'release_metadata.json'
download(API, metadata_path)
metadata = json.loads(metadata_path.read_text())
assert metadata['title'] == 'DepMap 24Q4 Public'
assert metadata['version'] == 1
wanted = ['Model.csv', 'README.txt', 'CRISPRGeneEffect.csv']
files = {f['name']: f for f in metadata['files']}
manifest = []
for name in wanted:
    source = files[name]
    dest = OUT / name
    expected = source['computed_md5']
    if not dest.exists() or md5(dest) != expected:
        print(f'Downloading {name}: {source["size"] / 1024**2:.1f} MiB', flush=True)
        download(source['download_url'], dest)
    actual = md5(dest)
    if actual != expected or dest.stat().st_size != source['size']:
        raise RuntimeError(f'Checksum/size validation failed for {name}')
    manifest.append({'file': name, 'url': source['download_url'],
                     'bytes': dest.stat().st_size, 'expected_md5': expected,
                     'actual_md5': actual})
    print(f'Validated {name}', flush=True)
record = {'release': metadata['title'], 'doi': metadata['doi'],
          'version': metadata['version'], 'metadata_api': API,
          'retrieved_utc': datetime.now(timezone.utc).isoformat(),
          'license': metadata['license'], 'files': manifest}
(OUT / 'download_manifest.json').write_text(json.dumps(record, indent=2))
print('Download and checksum validation complete.', flush=True)
