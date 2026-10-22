"""Offline exports and independent checks of the submitted numerical results."""
import csv
import json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.stats import fisher_exact
from .core import bh_adjust


def load_tables(root):
    return json.loads((Path(root)/'data/supplementary_tables.json').read_text())['tables']


def records(tables, number, block=0):
    b = next(t for t in tables if t['number'] == number)['blocks'][block]
    return [dict(zip(b['headers'], row)) for row in b['rows']]


def export_tables(root, out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    index = []
    for table in load_tables(root):
        for i, block in enumerate(table['blocks'], 1):
            name = f"S{table['number']:02d}_{i:02d}.csv"
            with (out/name).open('w', newline='') as f:
                w = csv.writer(f); w.writerow(block['headers']); w.writerows(block['rows'])
            index.append({'file': name, 'table': table['number'], 'question': table['question'],
                          'table_title': table['title'], 'block_title': block['title'],
                          'scope': table['scope'], 'notes': table['notes'], 'rows': len(block['rows'])})
    (out/'index.json').write_text(json.dumps(index, indent=2, ensure_ascii=False))
    print(f'Exported {len(index)} CSV blocks from 40 supplementary tables to {out}')


def validate(root, out):
    root, out = Path(root), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    tables = load_tables(root); checks = []
    def check(name, condition, detail=None):
        checks.append({'check': name, 'passed': bool(condition), 'detail': detail})
    def close(a, b): return bool(np.allclose(a, b, rtol=1e-8, atol=0))
    def rows(n, b=0): return records(tables, n, b)
    check('40 consecutively numbered supplementary tables', [t['number'] for t in tables] == list(range(1, 41)))
    check('Every table row matches its header width', all(len(r)==len(b['headers']) for t in tables for b in t['blocks'] for r in b['rows']))
    check('Q1: 55 genes; 1500 null pairs', len(rows(5))==55 and len(rows(6))==1500)
    check('Q1: all 55 gene BH values below 0.05', all(r['BH-adjusted P value'] < .05 for r in rows(5)))
    check('Q1: 34 distinct network gene nodes', len({r[k] for r in rows(8) for k in ('Source identifier','Target identifier') if str(r[k]).startswith('Gene:')})==34)
    q2=rows(10)
    shared={r['Gene symbol'] for r in q2 if r['Shared by at least two targets']}
    shown={r['Gene symbol'] for r in q2 if r['Shown in Figure 3C']}
    pd_only={'FBXO7','HTRA2','HTT','PARK7','PINK1','PRKN'}
    check('Q2: 27 shared genes and six additional PD mitophagy genes displayed',len(shared)==27 and shown==shared|pd_only)
    check('Q2: retained Figure 3B uses shared-set enrichment only',sum(bool(r['Shown in Figure 3B']) for r in rows(11))==20 and not any(r['Shown in Figure 3B'] for r in rows(11,1)))
    for summary in rows(9):
        terms=[r for r in rows(9,1) if r['Target identifier']==summary['Target identifier']]
        check('Q2 DWPC total: '+summary['Target disease'],close(sum(r['DWPC'] for r in terms),summary['Cumulative DWPC']))
    check('Q3: 15 genes, 36 routes, 8 disease endpoints',len({r['Gene symbol'] for r in rows(14)})==15 and len(rows(14))==36 and len({r['Disease identifier'] for r in rows(14)})==8)
    check('Q3: seven selected symptoms; 1891 total',len(rows(13))==7 and len(rows(13,1))==1891)
    check('Q4: 14 overlap genes, eight also in PM top 200',len(rows(16))==14 and sum(r['In premature menopause top 200']=='Yes' for r in rows(16))==8)
    check('Q4: generic dementia marker uses DOID:1307',len([r for r in rows(19) if str(r['Disease identifier']).endswith('DOID:1307')])==1)
    sets={}
    for r in rows(21): sets.setdefault(r['Seed group'],set()).add(r['Gene identifier'])
    check('Q5: 36 triple-overlap genes',len(set.intersection(*sets.values()))==36)
    for r in rows(22,2):
        odds,p=fisher_exact([[r['Intersection (a)'],r['First set only (b)']],[r['Second set only (c)'],r['Neither set (d)']]],alternative='greater')
        check('Q5 Fisher exact: '+r['Comparison'],close(odds,r['Odds ratio']) and close(p,r['P value']))
    raw=[json.loads(line) for line in (root/'data/q06/null_gene_counts.jsonl').read_text().splitlines()]
    check('Q6: all 1500 null pair identities and yields match Table 26', len(raw)==1500 and all(
        (x['pair'],x['source'],x['target'],x['unique_genes'])==(r['Random pair index'],r['Source disease identifier'],r['Target disease identifier'],r['Unique intermediate genes (n)'])
        for x,r in zip(raw,rows(26))))
    counts=Counter(n for path in json.loads((root/'data/q06/real_paths.json').read_text()) for n in path['nodes'][1:-1] if n.startswith('Gene:'))
    check('Q6: real path extraction reproduces 36 observed gene counts', len(counts)==36 and all(counts[r['Gene identifier']]==r['Observed recurrence'] for r in rows(25)))
    pvalues=[]
    for r in rows(25):
        c=np.array([x['counts'].get(r['Gene identifier'],0) for x in raw]); p=(1+np.count_nonzero(c>=r['Observed recurrence']))/(len(c)+1)
        pvalues.append(p)
        check('Q6 null replay: '+r['Gene symbol'],close(p,r['Empirical P value']) and close(c.mean(),r['Null mean']) and close(c.std(ddof=0),r['Null SD']))
    check('Q6: recomputed BH adjustment matches Table 25',close(bh_adjust(pvalues),[r['BH-adjusted P value'] for r in rows(25)]))
    check('Q6: 35 of 36 delirium–dementia genes occur in Q1',sum(bool(r['In Q1 55-gene set']) for r in rows(25))==35)
    for r in rows(28):
        # DES was not part of the refreshed five-compound query; its historical
        # score is retained in Table 28 without a matching route export.
        if r['Compound']=='Diethylstilbestrol': continue
        diseases={x['Disease identifier'] for x in rows(30) if x['Compound']==r['Compound'] and x['Disease identifier']}
        check('Q7 CGPS distinct disease count: '+r['Compound'],len(diseases)==r['Distinct cancer diseases (CGPS)'])
    check('Q7: 1000 random scores, 868 nonzero',len(rows(28,1))==1000 and sum(r['CGPS']>0 for r in rows(28,1))==868)
    drugsets={drug:{r['Gene symbol'] for r in rows(33,1) if r[drug]} for drug in ['Tamoxifen','Letrozole','Goserelin']}
    check('Q8: 927/174/237 drug-route gene memberships', [len(s) for s in drugsets.values()]==[927,174,237])
    check('Q8: 81 shared, 599/2/0 exclusive genes',len(set.intersection(*drugsets.values()))==81 and [len(s-set.union(*(v for k,v in drugsets.items() if k!=drug))) for drug,s in drugsets.items()]==[599,2,0])
    check('Q8: compound histogram accounts for 550596 nodes',sum(r['Compounds in bin (n)'] for r in rows(32,1))==550596)
    for r in rows(37):
        subset=[s for s in rows(37,1) if all(s[k]==r[k] for k in ['Source name','Target name','Route category'])]
        check('Q9: DWPC total '+r['Source name']+' / '+r['Target name']+' / '+r['Route category'],close(sum(s['DWPC'] for s in subset),r['Cumulative DWPC']))
    check('Q10: 14 queried pharmacogenes, eight connected, 28 edges',len(rows(39))==14 and len(rows(39,1))==8 and len(rows(40))==28)
    # Enrichr adjusts across returned terms. The larger full-library correction
    # is a separate recorded sensitivity check, not the displayed FDR.
    for q,number,block in [('q01',7,0),('q02',11,1),('q03',15,1)]:
        for path in sorted((root/'data/enrichment'/q).glob('*.json')):
            if path.name=='keywords.json': continue
            library=path.stem; rawterms=json.loads(path.read_text())[library]
            adjusted=bh_adjust([r[2] for r in rawterms]); lookup={r[1]:(r,adj) for r,adj in zip(rawterms,adjusted)}
            check(f'{q} {library}: Enrichr BH replay',close(adjusted,[r[6] for r in rawterms]))
            saved=[r for r in rows(number,block) if r['Library']==library]
            check(f'{q} {library}: published overlaps and adjusted P values match raw response',all(
                r['Term'] in lookup and r['Overlapping genes (n)']==len(lookup[r['Term']][0][5]) and close(r['Adjusted P value'],lookup[r['Term']][1]) for r in saved),{'published_terms':len(saved),'returned_terms':len(rawterms)})
    report={'checks':checks,'passed':sum(r['passed'] for r in checks),'failed':sum(not r['passed'] for r in checks),
            'scope':'Offline consistency and statistic replay; not a fresh SPOKE extraction.',
            'not_replayed':['Historical DES cancer-connectivity score: no route export in Table 30.',
                            'Q1 random-pair identities and Q7 random-compound identities were not saved.']}
    (out/'validation.json').write_text(json.dumps(report,indent=2))
    print(f"{report['passed']} checks passed; {report['failed']} failed")
    for c in checks:
        if not c['passed']: print('FAIL:',c['check'])
    if report['failed']: raise SystemExit(1)
    return report
