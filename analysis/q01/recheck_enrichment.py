# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Check the repeated Q1 Enrichr result and freeze the revised Figure 2D selection."""
import os
import csv, hashlib, json, math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
LIB = 'GO_Biological_Process_2023'
source_genes = set((OUT / 'Q1_input_genes_55.txt').read_text().split())
genes = {g.upper() for g in source_genes}
accepted = set(json.loads((OUT / 'raw/view.json').read_text())['genes'])
assert genes == accepted and len(genes) == 55
axes = json.loads((OUT / 'axis_keywords.json').read_text())
gmt = {}
for line in (OUT / f'libraries/{LIB}.gmt').read_text().splitlines():
    parts = line.split('\t')
    gmt[parts[0]] = {g.split(',')[0] for g in parts[2:] if g}
raw = json.loads((OUT / f'raw/{LIB}.json').read_text())[LIB]
old = {r['Term']: r for r in csv.DictReader((ROOT / 'audit_2026-09-30/excel_stage/Q1_GO_BP_55genes.csv').open())}

def bh(ps):
    order = sorted(range(len(ps)), key=ps.__getitem__)
    result = [1.0] * len(ps)
    current = 1.0
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        current = min(current, ps[i] * len(ps) / rank)
        result[i] = current
    return result

def fisher_tail(k, m, n, N=20000):
    return sum((math.comb(m, i) * math.comb(N - m, n - i) for i in range(k, min(m, n) + 1))) / math.comb(N, n)

def write(name, rows):
    with (OUT / name).open('w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
qs = bh([r[2] for r in raw])
rows = []
max_p_error = 0.0
max_q_error = 0.0
max_old_q_error = 0.0
for r, q in zip(raw, qs):
    rank, term, p, odds, combined, overlap, q_api, old_p, old_q = r
    overlap = set(overlap)
    assert overlap == genes & gmt[term], term
    k, m = (len(overlap), len(gmt[term]))
    max_p_error = max(max_p_error, abs(fisher_tail(k, m, 55) - p) / p)
    max_q_error = max(max_q_error, abs(q - q_api))
    tags = [a for a, keys in axes.items() if any((key in term.lower() for key in keys))] or ['Other']
    if term in old:
        saved = old[term]
        assert overlap == set(saved['Genes'].split(';'))
        assert k == int(saved['n_overlap']) == int(saved['Overlap'].split('/')[0])
        max_old_q_error = max(max_old_q_error, abs(q_api - float(saved['Adjusted P-value'])))
    rows.append({'Gene_set': LIB, 'Term': term, 'Overlap': f'{k}/{m}', 'P-value': p, 'Adjusted P-value': q_api, 'Odds Ratio': odds, 'Combined Score': combined, 'Genes': ';'.join(sorted(overlap)), 'n_overlap': k, 'significant (adj P<0.05)': q_api < 0.05, 'axis (keyword rule, dwpc.ipynb cell 32)': '; '.join(tags), 'at_least_2_overlap_genes': k >= 2})
assert max_p_error < 0.001 and max_q_error < 1e-10 and (max_old_q_error < 1e-12)
rows.sort(key=lambda r: (r['Adjusted P-value'], r['P-value'], r['Term']))
retained = [r for r in rows if r['n_overlap'] >= 2]
sig = [r for r in retained if r['Adjusted P-value'] < 0.05]
assert len(retained) == 654 and len(sig) == 583
assert {r['Term'] for r in retained} == set(old)
counts = Counter()
for r in sig:
    counts.update(r['axis (keyword rule, dwpc.ipynb cell 32)'].split('; '))
assert counts == Counter({'Neurodegenerative': 57, 'Inflammatory': 49, 'Vascular': 19, 'Endocrine': 4, 'Other': 456})
endo = [r for r in sig if 'Endocrine' in r['axis (keyword rule, dwpc.ipynb cell 32)']]
vascular = [r for r in sig if 'Vascular' in r['axis (keyword rule, dwpc.ipynb cell 32)']]
selected = {r['Term']: dict(r, selection_reason='Ten smallest adjusted P values overall') for r in sig[:10]}
for r in vascular[:2]:
    selected[r['Term']] = dict(r, selection_reason='Two smallest adjusted P values among vascular-tagged terms')
for r in endo:
    selected[r['Term']] = dict(r, selection_reason='All four significant endocrine-tagged terms')
selected = sorted(selected.values(), key=lambda r: (r['Adjusted P-value'], r['P-value'], r['Term']))
assert len(selected) == 16
palette = {'Neurodegenerative': '#598EAF', 'Inflammatory': '#B77798', 'Vascular': '#7B61A8', 'Endocrine': '#C59B55', 'Other': '#BAC3CA'}
labels = {'GO:0001932': 'Regulation of protein phosphorylation', 'GO:0001774': 'Microglial cell activation', 'GO:0048143': 'Astrocyte activation', 'GO:0051247': 'Positive regulation of\nprotein metabolic process', 'GO:0043523': 'Regulation of neuronal apoptosis', 'GO:0001934': 'Positive regulation of\nprotein phosphorylation', 'GO:1990000': 'Amyloid fibril formation', 'GO:0014002': 'Astrocyte development', 'GO:1901216': 'Positive regulation of neuron death', 'GO:0043066': 'Negative regulation of apoptosis', 'GO:0090276': 'Regulation of peptide hormone secretion', 'GO:0050796': 'Regulation of insulin secretion', 'GO:1901550': 'Regulation of endothelial cell development', 'GO:0031960': 'Response to corticosteroid', 'GO:1904994': 'Regulation of leukocyte adhesion\nto vascular endothelial cells', 'GO:0046626': 'Regulation of insulin receptor signaling'}
for r in selected:
    go = r['Term'].split('(')[-1].rstrip(')')
    r['GO_ID'] = go
    r['Plot_label'] = labels[go]
    r['Bar_colors'] = '; '.join((palette[a] for a in r['axis (keyword rule, dwpc.ipynb cell 32)'].split('; ')))
    r['minus_log10_adjusted_P'] = -math.log10(r['Adjusted P-value'])
selected_terms = {r['Term'] for r in selected}
for r in rows:
    r['shown_in_original_Fig2d_top20'] = r['Term'] in old and old[r['Term']]['shown_in_Fig2d (top 20)'] == 'True'
    r['shown_in_revised_Fig2d_representative'] = r['Term'] in selected_terms
write('Q1_Enrichr_all_1541_terms.csv', rows)
write('Q1_GO_BP_55genes_verified.csv', retained)
write('Q1_GO_BP_55genes_significant.csv', sig)
write('Q1_endocrine_terms.csv', endo)
write('Figure2D_plotted_terms.csv', selected)
write('Q1_axis_counts.csv', [{'Keyword_category': a, 'Significant_terms_n_ge_2': counts[a], 'Keywords': '; '.join(axes.get(a, [])), 'Nonexclusive_categories': True} for a in [*palette]])
validation = {'input_genes': 55, 'accepted_genes_identical_after_case_normalization': True, 'case_normalization': {g: g.upper() for g in source_genes if g != g.upper()}, 'raw_returned_terms': len(rows), 'significant_before_overlap_filter': sum((r['Adjusted P-value'] < 0.05 for r in rows)), 'retained_terms_n_ge_2': 654, 'significant_terms_n_ge_2': 583, 'axis_counts': dict(counts), 'multi_category_significant_terms': sum(('; ' in r['axis (keyword rule, dwpc.ipynb cell 32)'] for r in sig)), 'all_saved_GO_terms_and_gene_overlaps_reproduced': True, 'max_adjusted_P_difference_from_saved': max_old_q_error, 'all_overlaps_match_GMT': True, 'max_relative_Fisher_P_error_N_20000': max_p_error, 'max_absolute_BH_error_returned_terms': max_q_error, 'multiple_testing': 'Enrichr adjusted P values; BH verified across the 1,541 returned overlapping terms before the n>=2 display/summary filter.', 'plot_selection': '16 representative terms: ten smallest adjusted P values overall, two strongest vascular-tagged terms, all four significant endocrine-tagged terms; sorted by adjusted P then raw P then name.', 'category_note': 'Original term-name substring rules, with overlapping categories. Other means no keyword match, not absence of biological relevance.', 'library_terms': len(gmt), 'genes_with_library_annotation': len(genes & set().union(*gmt.values())), 'verified_at_utc': datetime.now(timezone.utc).isoformat()}
(OUT / 'verification.json').write_text(json.dumps(validation, indent=2) + '\n')
meta = json.loads((OUT / 'provenance.json').read_text())
meta.update({'submission': json.loads((OUT / 'raw/addList.json').read_text()), 'retrieved_at_utc': datetime.now(timezone.utc).isoformat(), 'raw_sha256': {str(p.relative_to(OUT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [*sorted((OUT / 'raw').glob('*')), *sorted((OUT / 'libraries').glob('*'))]}, 'reproduction': 'All 654 archived rows, overlaps and adjusted P values reproduced; figures use the repeated run.'})
(OUT / 'provenance.json').write_text(json.dumps(meta, indent=2) + '\n')
print(json.dumps(validation, indent=2))
