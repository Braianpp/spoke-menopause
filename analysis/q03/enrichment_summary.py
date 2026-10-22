# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Reconstruct and validate the Q3 enrichment tables from saved Enrichr responses."""
import os
import csv
import hashlib
import json
import math
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from run_enrichment import ROOT, LIBRARIES

def bh(ps):
    order = sorted(range(len(ps)), key=ps.__getitem__)
    q = [1.0] * len(ps)
    running = 1.0
    for rank in range(len(ps), 0, -1):
        i = order[rank - 1]
        running = min(running, ps[i] * len(ps) / rank)
        q[i] = running
    return q

def hypergeom(k, m, n, N=20000):
    return sum((math.comb(m, x) * math.comb(N - m, n - x) for x in range(k, min(m, n) + 1))) / math.comb(N, n)

def save_csv(name, rows):
    with (ROOT / name).open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def main():
    genes = set((ROOT / 'Q3_input_genes.txt').read_text().split())
    accepted = set(json.loads((ROOT / 'raw/view.json').read_text())['genes'])
    assert genes == accepted and len(genes) == 15
    rows, summaries, checks = ([], [], [])
    for lib in LIBRARIES:
        gmt = {}
        for line in (ROOT / 'libraries' / f'{lib}.gmt').read_text().splitlines():
            fields = line.split('\t')
            gmt[fields[0]] = {g.split(',')[0] for g in fields[2:] if g}
        data = json.loads((ROOT / 'raw' / f'{lib}.json').read_text())[lib]
        ps = [x[2] for x in data]
        qs = bh(ps)
        full_qs = bh(ps + [1.0] * (len(gmt) - len(data)))[:len(data)]
        max_p_error = 0.0
        max_q_error = 0.0
        for item, q, q_full in zip(data, qs, full_qs):
            rank, term, p, odds, combined, overlap, q_api, old_p, old_q = item
            overlap = set(overlap)
            assert overlap == genes & gmt[term], (lib, term, overlap, genes & gmt[term])
            m, k = (len(gmt[term]), len(overlap))
            p_check = hypergeom(k, m, len(genes))
            max_p_error = max(max_p_error, abs(p_check - p) / p)
            max_q_error = max(max_q_error, abs(q - q_api))
            rows.append({'Library': lib, 'Term': term, 'Enrichr_rank': rank, 'Overlap': f'{k}/{m}', 'Overlap_count': k, 'Term_gene_count': m, 'Input_gene_count': len(genes), 'P-value': p, 'Adjusted P-value': q_api, 'Odds Ratio': odds, 'Combined Score': combined, 'Genes': ';'.join(sorted(overlap)), 'Significant_Enrichr_adjusted_P_lt_0.05': q_api < 0.05, 'Full_library_BH_sensitivity': q_full})
        assert max_p_error < 0.001, (lib, max_p_error)
        assert max_q_error < 1e-10, (lib, max_q_error)
        coverage = genes & set().union(*gmt.values())
        summaries.append({'Library': lib, 'Total_library_terms': len(gmt), 'Returned_terms_with_overlap': len(data), 'Significant_terms_Enrichr': sum((x[6] < 0.05 for x in data)), 'Significant_terms_full_library_BH_sensitivity': sum((q < 0.05 for q in full_qs)), 'Input_genes_annotated_in_library': len(coverage), 'Input_genes_without_annotation': ';'.join(sorted(genes - coverage))})
        checks.append({'library': lib, 'all_overlap_genes_match_saved_GMT': True, 'max_relative_Fisher_P_error_N_20000': max_p_error, 'max_absolute_BH_error_over_returned_terms': max_q_error})
    rows.sort(key=lambda r: (r['Adjusted P-value'], r['Library'], r['Term']))
    save_csv('Q3_Enrichr_all_results.csv', rows)
    save_csv('Q3_Enrichr_significant_results.csv', [r for r in rows if r['Adjusted P-value'] < 0.05])
    save_csv('Q3_library_summary.csv', summaries)
    claim_patterns = 'neurodegeneration|microglial|astrocyte activation|innate immune system|neuron apoptotic process|post-translational protein modification|protein phosphorylation|positive regulation of inflammatory response|positive regulation of i-kappab kinase/nf-kappab|ubiquitin-dependent protein catabolic|mitophagy|amyloid fiber formation'
    save_csv('Q3_sentence_support.csv', [r for r in rows if re.search(claim_patterns, r['Term'], re.I)])
    select = [('KEGG_2021_Human', 'Pathways of neurodegeneration', 'Pathways of neurodegeneration'), ('GO_Biological_Process_2023', '(GO:0001934)', 'Positive regulation of\nprotein phosphorylation'), ('GO_Biological_Process_2023', '(GO:0031648)', 'Protein destabilization'), ('GO_Biological_Process_2023', '(GO:1990535)', 'Neuron projection maintenance'), ('GO_Biological_Process_2023', '(GO:1901214)', 'Regulation of neuron death'), ('GO_Biological_Process_2023', '(GO:1903146)', 'Regulation of mitochondrial\nautophagy'), ('GO_Biological_Process_2023', '(GO:2000060)', 'Positive regulation of ubiquitin-\ndependent protein catabolism'), ('GO_Biological_Process_2023', '(GO:0043523)', 'Regulation of neuronal apoptosis'), ('GO_Biological_Process_2023', '(GO:0001774)', 'Microglial cell activation'), ('GO_Biological_Process_2023', '(GO:0043123)', 'Positive regulation of\nIκB kinase/NF-κB signaling'), ('Reactome_2022', 'Amyloid Fiber Formation R-HSA-977225', 'Amyloid fiber formation'), ('GO_Biological_Process_2023', '(GO:0050729)', 'Positive regulation of\ninflammatory response')]
    selected = []
    for lib, match, label in select:
        hits = [r for r in rows if r['Library'] == lib and match in r['Term']]
        assert len(hits) == 1
        row = dict(hits[0], Plot_label=label)
        assert row['Adjusted P-value'] < 0.05 and row['Overlap_count'] >= 2
        selected.append(row)
    selected.sort(key=lambda r: r['Adjusted P-value'])
    save_csv('Figure4C_plotted_terms.csv', selected)
    (ROOT / 'validation.json').write_text(json.dumps({'input_equals_accepted_Enrichr_list': True, 'unique_genes': 15, 'libraries': checks, 'rows': len(rows), 'significant_rows': sum((r['Adjusted P-value'] < 0.05 for r in rows)), 'plot_selection': '12 representative significant terms, at least 2 overlapping input genes, sorted by Enrichr adjusted P-value; not the top 12 overall', 'BH_note': 'API adjusted P-values reproduce BH across returned overlapping terms per library. Full-library BH (zero-overlap terms P=1) is supplied as a separately labeled sensitivity calculation; it is not used for the main plot.', 'full_library_sensitivity_note': 'Microglial activation, neuronal apoptosis, protein phosphorylation and NF-kappaB signaling remain below 0.05. Positive regulation of inflammatory response and astrocyte activation do not.'}, indent=2) + '\n')
    for source, target in [('/private/tmp/q3_enrichr_libraries.json', 'raw/datasetStatistics.json'), ('/private/tmp/q3_enrichr_api.html', 'raw/api_documentation.html'), ('/private/tmp/q3_enrichr_background.html', 'raw/background_documentation.html')]:
        shutil.copyfile(source, ROOT / target)
    provenance = json.loads((ROOT / 'provenance.json').read_text())
    provenance.update({'submission': json.loads((ROOT / 'raw/addList.json').read_text()), 'submitted_at_utc': datetime.fromtimestamp((ROOT / 'raw/addList.json').stat().st_mtime, timezone.utc).isoformat(), 'retrieved_at_utc': datetime.fromtimestamp(max((p.stat().st_mtime for p in (ROOT / 'libraries').glob('*.gmt'))), timezone.utc).isoformat(), 'api_client': 'curl; direct Enrichr API (not gseapy)', 'background_size_verified_against_Fisher_P_values': 20000, 'correction': 'Enrichr BH across returned terms separately per library; all-library-term BH supplied as sensitivity', 'raw_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [*sorted((ROOT / 'raw').glob('*')), *sorted((ROOT / 'libraries').glob('*'))]}})
    (ROOT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(summaries, indent=2))
    print('Validated every overlap and BH value; independent Fisher P-values agree within 0.1% relative tolerance.')
if __name__ == '__main__':
    main()
