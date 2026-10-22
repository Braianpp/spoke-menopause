# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Submit the Figure 4b gene union to Enrichr; retain raw responses and libraries.

Usage: bundled-python run_enrichment.py prepare|submit|fetch
Then run summarize_enrichment.py to reconstruct tables without web requests.
The submit step is explicit so offline reconstruction never creates a new list.
"""
import os
import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'figure4b_2026-09-30' / 'Figure4B_nodes.csv'
BASE = 'https://maayanlab.cloud/Enrichr'
LIBRARIES = ['GO_Biological_Process_2023', 'Reactome_2022', 'KEGG_2021_Human', 'WikiPathway_2023_Human']

def curl(path, *args):
    subprocess.run(['curl', '--fail', '--silent', '--show-error', '--max-time', '90', *args, '--output', str(path)], check=True)

def prepare():
    rows = list(csv.DictReader(SOURCE.open()))
    genes = sorted((r['name'] for r in rows if r['node_type'] == 'Gene'))
    assert len(genes) == len(set(genes)) == 15
    (ROOT / 'raw').mkdir(exist_ok=True)
    (ROOT / 'libraries').mkdir(exist_ok=True)
    (ROOT / 'Q3_input_genes.txt').write_text('\n'.join(genes) + '\n')
    manifest = {'prepared_at_utc': datetime.now(timezone.utc).isoformat(), 'source_file': str(SOURCE), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(), 'gene_count': len(genes), 'genes': genes, 'libraries': LIBRARIES, 'service': BASE, 'background': 'Enrichr default (no custom background supplied)', 'significance': 'Enrichr adjusted P-value < 0.05, separately per library', 'api_reference': 'https://maayanlab.cloud/Enrichr/help#api', 'input_scope': 'All 15 unique genes on direct insomnia/paresthesia–gene–dementia-related disease routes in Fig. 4b'}
    (ROOT / 'provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Prepared:', ', '.join(genes))

def submit():
    dest = ROOT / 'raw' / 'addList.json'
    if dest.exists():
        raise RuntimeError('Submission already exists; use fetch or summarize.')
    curl(dest, '-X', 'POST', '-F', 'list=<' + str(ROOT / 'Q3_input_genes.txt'), '-F', 'description=Q3 symptom-associated gene union (15 genes)', BASE + '/addList')
    reply = json.loads(dest.read_text())
    assert 'userListId' in reply
    provenance = json.loads((ROOT / 'provenance.json').read_text())
    provenance['submitted_at_utc'] = datetime.now(timezone.utc).isoformat()
    provenance['submission'] = reply
    (ROOT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print('Submitted gene list:', reply)

def fetch():
    uid = json.loads((ROOT / 'raw' / 'addList.json').read_text())['userListId']
    curl(ROOT / 'raw' / 'view.json', BASE + f'/view?userListId={uid}')
    curl(ROOT / 'raw' / 'datasetStatistics.json', BASE + '/datasetStatistics')
    for lib in LIBRARIES:
        curl(ROOT / 'raw' / f'{lib}.json', BASE + f'/enrich?userListId={uid}&backgroundType={lib}')
        curl(ROOT / 'libraries' / f'{lib}.gmt', BASE + f'/geneSetLibrary?mode=text&libraryName={lib}')
        data = json.loads((ROOT / 'raw' / f'{lib}.json').read_text())
        print(lib, 'returned', len(data[lib]), 'terms', flush=True)
    provenance = json.loads((ROOT / 'provenance.json').read_text())
    provenance['retrieved_at_utc'] = datetime.now(timezone.utc).isoformat()
    provenance['raw_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [*sorted((ROOT / 'raw').glob('*')), *sorted((ROOT / 'libraries').glob('*'))]}
    (ROOT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
if __name__ == '__main__':
    {'prepare': prepare, 'submit': submit, 'fetch': fetch}[sys.argv[1]]()
