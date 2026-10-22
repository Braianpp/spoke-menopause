# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Validate the fresh full-PD GO result against input membership and the GMT.

Run with /opt/anaconda3/bin/python3; requires scipy for an independent Fisher test.
Raw API responses are retained alongside this script.
"""
import os
from pathlib import Path
import csv
import hashlib
import json
import math
import re
from scipy.stats import hypergeom
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
LIBRARY = 'GO_Biological_Process_2023'

def bh(ps, n_tests=None):
    order = sorted(range(len(ps)), key=ps.__getitem__)
    qs = [1.0] * len(ps)
    running = 1.0
    for rank in range(len(ps), 0, -1):
        idx = order[rank - 1]
        running = min(running, ps[idx] * (n_tests or len(ps)) / rank)
        qs[idx] = running
    return qs

def csvout(name, rows):
    with (OUT / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
source = ROOT / 'audit_2026-09-30/excel_stage/Q2_IntermediateGenes_242.csv'
all_genes = list(csv.DictReader(source.open()))
pd_rows = [r for r in all_genes if 'Parkinson' in r['targets']]
genes = set((OUT / 'PD_165_genes.txt').read_text().split())
assert len(genes) == 165 == len(pd_rows)
assert genes == {r['gene_name'] for r in pd_rows}
shared = {r['gene_name'] for r in all_genes if int(r['n_targets']) >= 2}
assert len(shared) == 27
assert genes & shared == {'APOE', 'APP', 'MAPT', 'SNCA', 'SOD1', 'TARDBP'}
assert not genes & {'POLG', 'POLG2', 'TWNK'}
assert {g.upper() for g in json.loads((OUT / 'view.json').read_text())['genes']} == genes
gmt_path = OUT / f'{LIBRARY}.gmt'
cached = ROOT / 'figures/figure2d_2026-10-01/libraries' / f'{LIBRARY}.gmt'
assert gmt_path.read_bytes() == cached.read_bytes(), 'GO library changed since Figure 2 rerun'
gmt = {}
for line in gmt_path.read_text().splitlines():
    fields = line.split('\t')
    gmt[fields[0]] = {g.split(',')[0].upper() for g in fields[2:] if g}
raw = json.loads((OUT / 'GO_raw.json').read_text())[LIBRARY]
assert len({r[1] for r in raw}) == len(raw)
assert {r[1] for r in raw} == {term for term, symbols in gmt.items() if genes & symbols}
qs = bh([r[2] for r in raw])
full_qs = bh([r[2] for r in raw], len(gmt))
result = []
max_p_error = 0.0
max_q_error = 0.0
pattern = re.compile('mitochond|mitophagy|oxidative phosphorylation|respiratory.*(chain|electron)', re.I)
for r, q, full_q in zip(raw, qs, full_qs):
    rank, term, p, odds, combined, overlap, api_q, *_ = r
    expected = genes & gmt[term]
    assert expected == {g.upper() for g in overlap}, term
    independently_computed_p = float(hypergeom.sf(len(expected) - 1, 20000, len(gmt[term]), len(genes)))
    max_p_error = max(max_p_error, abs(independently_computed_p - p) / p)
    max_q_error = max(max_q_error, abs(api_q - q))
    result.append({'library': LIBRARY, 'term': term, 'input_genes': len(genes), 'overlap_count': len(expected), 'term_gene_count': len(gmt[term]), 'p_value': p, 'adjusted_p_value': api_q, 'full_library_BH_sensitivity': full_q, 'genes': ';'.join(sorted(expected)), 'overlap_genes_shown_in_Fig3C': ';'.join(sorted(expected & shared)), 'overlap_genes_excluded_from_Fig3C': ';'.join(sorted(expected - shared)), 'significant_with_2plus_genes': api_q < 0.05 and len(expected) >= 2, 'mitochondrial_term_name_screen': bool(pattern.search(term))})
assert max_p_error < 0.002, max_p_error
assert max_q_error < 1e-10, max_q_error
result.sort(key=lambda r: (r['adjusted_p_value'], r['p_value'], r['term']))
mito = [r for r in result if r['mitochondrial_term_name_screen']]
csvout('PD_165_GO_all_results.csv', result)
csvout('PD_165_mitochondrial_GO_terms.csv', mito)
claim_ids = ['GO:1903146', 'GO:0090140', 'GO:0051881']
claims = [r for r in result if any((go in r['term'] for go in claim_ids))]
assert len(claims) == 3
assert all((r['full_library_BH_sensitivity'] < 0.05 for r in claims))
autophagy = next((r for r in claims if 'GO:1903146' in r['term']))
assert not autophagy['overlap_genes_shown_in_Fig3C']
csvout('PD_165_selected_mitochondrial_terms.csv', claims)
report = {'PD_input_genes': len(genes), 'PD_genes_shared_with_other_targets': len(genes & shared), 'PD_genes_excluded_from_shared_gene_panels': len(genes - shared), 'input_symbols_accepted_exactly': True, 'library_identical_to_October_1_copy': True, 'library_sha256': hashlib.sha256(gmt_path.read_bytes()).hexdigest(), 'GO_library_terms': len(gmt), 'PD_genes_annotated_in_library': len(genes & set.union(*gmt.values())), 'returned_terms': len(raw), 'all_overlapping_terms_returned': True, 'all_gene_overlaps_verified': True, 'significant_terms_with_2plus_genes': sum((r['significant_with_2plus_genes'] for r in result)), 'mitochondrial_term_screen_regex': pattern.pattern, 'mitochondrial_screen_terms': len(mito), 'mitochondrial_screen_significant_2plus_genes': sum((r['significant_with_2plus_genes'] for r in mito)), 'max_relative_Fisher_p_error': max_p_error, 'max_absolute_BH_error': max_q_error, 'selected_claim_terms': claims, 'interpretation': 'New exploratory enrichment of the full saved PD set. The current Figure 3B enrichment remains the 27 shared genes. This analysis tests GO annotation overrepresentation against the Enrichr default background; it does not establish menopause-specific enrichment or a causal menopause-to-PD mechanism. The route extraction is capped at top five schemas and 200 paths/schema.'}
(OUT / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
provenance = json.loads((OUT / 'provenance.json').read_text())
provenance.update({'Enrichr_userListId': json.loads((OUT / 'addList.json').read_text())['userListId'], 'validated': True, 'GO_returned_terms': len(raw), 'BH_correction': 'All 2021 returned overlapping GO terms; additionally verified full-5407-term library BH sensitivity', 'validation_file': 'validation.json'})
(OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
print(json.dumps(report, indent=2))
