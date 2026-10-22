# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
Cancer Gene Proximity Score (CGPS) — final version.

Score = number of distinct cancer diseases reachable from a compound's
        protein targets, through any of four Gene–Disease edge types
        in SPOKE.

For each compound:
  1. Identify target genes: Compound –BINDS→ Protein ←ENCODES– Gene
  2. Collect all DISTINCT cancer diseases connected to those genes via:
       ASSOCIATES_DaG          (curated gene–disease associations)
       EXPRESSEDIN_GeiD        (gene expressed in disease context)
       MARKER_POS_GmpD         (positive disease biomarker)
       MARKER_NEG_GmnD         (negative disease biomarker)
  3. Score = count of distinct cancer diseases reached.

Cancer diseases are identified by case-insensitive substring matching
of disease names against: cancer, carcinoma, neoplasm, tumor, tumour,
lymphoma, leukemia, myeloma, sarcoma.

The score is contextualized against a background of N randomly sampled
compounds (with at least one protein-binding edge) scored identically.

Why this metric:
  - Single number, no arbitrary component weights.
  - Directly interpretable: "how many cancer diseases does this
    compound's target neighborhood touch?"
  - Compounds with more targets or targets more deeply embedded in
    cancer biology score higher — this is exactly the signal we want.
  - Background distribution provides empirical p-value.

Usage:
  python cancer_proximity_final.py     --compounds "Fezolinetant,Elinzanetant,Tamoxifen,..."     --random 1000     --out-json results/cgps_final.json     --out-fig  results/cgps_final.png     --out-legend results/cgps_final_legend.png
"""
import os
import os, sys, json, argparse, math
import numpy as np
CANCER_KWS = ['cancer', 'carcinoma', 'neoplasm', 'tumor', 'tumour', 'lymphoma', 'leukemia', 'myeloma', 'sarcoma']
GENE_DISEASE_EDGES = ['ASSOCIATES_DaG', 'EXPRESSEDIN_GeiD', 'MARKER_POS_GmpD', 'MARKER_NEG_GmnD']

def load_config(path):
    cfg = {}
    if not os.path.exists(path):
        return cfg
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if ':' in line:
                k, v = line.split(':', 1)
                cfg[k.strip()] = v.strip()
    return cfg

def get_driver(uri, user, password):
    from neo4j import GraphDatabase
    return GraphDatabase.driver(uri, auth=(user, password))

def resolve_compound(tx, name_or_id):
    q = 'MATCH (c:Compound) WHERE toLower(c.identifier) = $q OR toLower(c.name) = $q    OR any(s IN coalesce(c.synonyms,[]) WHERE toLower(s) = $q) RETURN c.identifier AS id, c.name AS name LIMIT 1'
    rec = tx.run(q, q=name_or_id.lower()).single()
    if rec:
        return dict(rec)
    q2 = 'MATCH (c:Compound) WHERE toLower(c.name) CONTAINS $q RETURN c.identifier AS id, c.name AS name LIMIT 1'
    rec2 = tx.run(q2, q=name_or_id.lower()).single()
    return dict(rec2) if rec2 else None

def count_cancer_diseases(tx, kws):
    q = 'WITH $kws AS kws MATCH (d:Disease) WHERE any(kw IN kws WHERE toLower(d.name) CONTAINS kw) RETURN count(d) AS n'
    return tx.run(q, kws=kws).single()['n']

def compute_score(tx, comp_id, kws):
    """
    Count distinct cancer diseases reachable from compound's target genes.
    Returns n_targets, n_cancer_genes, n_cancer_diseases (= the score),
    and the list of cancer disease names.
    """
    q = 'WITH $kws AS kws MATCH (c:Compound {identifier: $id})-[:BINDS_CbP]->(p:Protein) WITH p, kws MATCH (g:Gene)-[:ENCODES_GeP]->(p) WITH DISTINCT g, kws OPTIONAL MATCH (g)-[:ASSOCIATES_DaG|EXPRESSEDIN_GeiD|MARKER_POS_GmpD|MARKER_NEG_GmnD]-(d:Disease) WHERE any(kw IN kws WHERE toLower(d.name) CONTAINS kw) WITH g, collect(DISTINCT d.name) AS cancer_diseases RETURN g.name AS gname,        size(cancer_diseases) AS n_cd,        cancer_diseases'
    rows = [dict(r) for r in tx.run(q, id=comp_id, kws=kws)]
    if not rows:
        return {'n_targets': 0, 'n_cancer_genes': 0, 'n_cancer_diseases': 0, 'genes': []}
    all_diseases = set()
    gene_details = []
    n_cancer_genes = 0
    for r in rows:
        n_cd = r['n_cd']
        if n_cd > 0:
            n_cancer_genes += 1
            for d in r['cancer_diseases']:
                all_diseases.add(d)
        gene_details.append({'gname': r['gname'], 'n_cancer_diseases': n_cd})
    return {'n_targets': len(rows), 'n_cancer_genes': n_cancer_genes, 'n_cancer_diseases': len(all_diseases), 'genes': gene_details}

def compute_score_fast(tx, comp_id, kws):
    """Fast version for random background — single number."""
    q = 'WITH $kws AS kws MATCH (c:Compound {identifier: $id})-[:BINDS_CbP]->(p:Protein) WITH p, kws MATCH (g:Gene)-[:ENCODES_GeP]->(p) WITH DISTINCT g, kws MATCH (g)-[:ASSOCIATES_DaG|EXPRESSEDIN_GeiD|MARKER_POS_GmpD|MARKER_NEG_GmnD]-(d:Disease) WHERE any(kw IN kws WHERE toLower(d.name) CONTAINS kw) RETURN count(DISTINCT d) AS n_cancer_diseases'
    rec = tx.run(q, id=comp_id, kws=kws).single()
    return rec['n_cancer_diseases'] if rec else 0

def sample_random_compounds(tx, n, exclude_ids):
    q = 'MATCH (c:Compound)-[:BINDS_CbP]->(:Protein) WITH DISTINCT c WHERE NOT c.identifier IN $ex RETURN c.identifier AS id, c.name AS name ORDER BY rand() LIMIT $n'
    return [dict(r) for r in tx.run(q, n=n, ex=exclude_ids)]

def make_plots(panel_results, random_scores, out_fig, out_legend):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    COLORS = {'Fezolinetant': '#d62728', 'Elinzanetant': '#1f77b4', 'Goserelin': '#e377c2', 'Letrozole': '#9467bd', 'Tamoxifen': '#17becf', 'Diethylstilbestrol': '#2ca02c'}
    ORDER = ['Fezolinetant', 'Elinzanetant', 'Goserelin', 'Letrozole', 'Tamoxifen', 'Diethylstilbestrol']
    items = []
    for r in panel_results:
        name = r['compound']['name']
        score = r['score']
        key = None
        for k in COLORS:
            if k.lower() in name.lower():
                key = k
                break
        if key:
            items.append((key, score, COLORS[key]))
    seen = {}
    for k, s, c in items:
        if k not in seen:
            seen[k] = (s, c)
    ordered = [(k, seen[k][0], seen[k][1]) for k in ORDER if k in seen]
    rnd_nz = random_scores[random_scores > 0]
    n_nz = len(rnd_nz)
    TICK_SIZE = 14
    TICK_WEIGHT = 'bold'
    fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    ax.hist(rnd_nz, bins=40, color='lightgray', edgecolor='gray')
    nonzero_items = [(n, s, c) for n, s, c in ordered if s > 0]
    zero_items = [(n, s, c) for n, s, c in ordered if s == 0]
    for name, score, color in nonzero_items:
        ax.axvline(score, color=color, linewidth=3.0, alpha=0.95)
    if zero_items:
        ylim = ax.get_ylim()
        y_top = ylim[1]
        for i, (name, score, color) in enumerate(zero_items):
            y_pos = y_top * (0.92 - i * 0.08)
            ax.annotate(name, xy=(0.5, y_pos), fontsize=9, fontweight='bold', color=color, ha='left', va='center')
            ax.plot(0.5, y_pos, marker='<', markersize=7, color=color, clip_on=False, zorder=5)
    mu = float(rnd_nz.mean()) if n_nz > 0 else 0.0
    ax.axvline(mu, color='#f2c744', lw=3, ls='--')
    ax.set_xlabel('Distinct Cancer Diseases Reached', fontsize=13)
    ax.set_ylabel('Count', fontsize=13)
    ax.set_title('Cancer Gene Proximity Score — Compound Comparison\n(Active Compounds Only, n={})'.format(n_nz), fontsize=12)
    ax.tick_params(axis='both', labelsize=TICK_SIZE)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight(TICK_WEIGHT)
    fig.savefig(out_fig, dpi=220)
    plt.close(fig)
    handles = []
    labels_list = []
    handles.append(Line2D([0], [0], color='lightgray', lw=6))
    labels_list.append('Non-zero Random (n={})'.format(n_nz))
    for name, score, color in ordered:
        handles.append(Line2D([0], [0], color=color, lw=4))
        labels_list.append(name)
    handles.append(Line2D([0], [0], color='#f2c744', lw=4, ls='--'))
    labels_list.append('Mean')
    h = max(2.5, 0.5 * len(labels_list) + 1.0)
    fig2, ax2 = plt.subplots(figsize=(3.4, h), constrained_layout=True)
    ax2.legend(handles=handles, labels=labels_list, loc='center left', ncol=1, frameon=False, prop={'weight': 'bold', 'size': 14})
    ax2.axis('off')
    fig2.patch.set_alpha(0.0)
    fig2.savefig(out_legend, dpi=240, transparent=True)
    plt.close(fig2)

def main():
    ap = argparse.ArgumentParser(description='Cancer Gene Proximity Score (CGPS) via SPOKE')
    ap.add_argument('--compounds', default='Fezolinetant,inchikey:DWRIJNIPBUFCQS-DQEYMECFSA-N,Tamoxifen,Diethylstilbestrol,Letrozole,Goserelin')
    ap.add_argument('--random', type=int, default=1000)
    ap.add_argument('--out-json', default='')
    ap.add_argument('--out-fig', default='')
    ap.add_argument('--out-legend', default='')
    ap.add_argument('--config', default=os.path.join(os.path.dirname(__file__), '..', 'config', 'spoke.yaml'))
    args = ap.parse_args()
    names = [s.strip() for s in args.compounds.split(',') if s.strip()]
    kws = [k.lower() for k in CANCER_KWS]
    cfg = load_config(args.config)
    uri = os.getenv('SPOKE_URI') or cfg.get('bolt_uri')
    user = os.getenv('SPOKE_USER') or cfg.get('user')
    password = os.getenv('SPOKE_PASSWORD') or cfg.get('password')
    database = os.getenv('SPOKE_DATABASE') or cfg.get('database')
    driver = get_driver(uri, user, password)
    try:
        with driver.session(database=database) as sess:
            n_cancer_total = sess.execute_read(count_cancer_diseases, kws)
            panel_results = []
            panel_ids = []
            for nm in names:
                comp = sess.execute_read(resolve_compound, nm)
                if not comp:
                    print('WARN: not found:', nm, file=sys.stderr)
                    continue
                detail = sess.execute_read(compute_score, comp['id'], kws)
                score = detail['n_cancer_diseases']
                panel_results.append({'compound': comp, 'score': score, 'detail': {k: v for k, v in detail.items() if k != 'genes'}, 'genes': detail['genes']})
                panel_ids.append(comp['id'])
                print('  {:25s}  targets={:3d}  cancer_genes={:3d}  cancer_diseases={:4d}'.format(comp['name'], detail['n_targets'], detail['n_cancer_genes'], score))
            print('\nScoring {} random compounds...'.format(args.random))
            randoms = sess.execute_read(sample_random_compounds, args.random, panel_ids)
            random_scores = []
            for i, rc in enumerate(randoms):
                s = sess.execute_read(compute_score_fast, rc['id'], kws)
                random_scores.append(s)
                if (i + 1) % 200 == 0:
                    print('  ...{}/{}'.format(i + 1, len(randoms)))
            random_scores = np.array(random_scores, dtype=float)
            rnd_nz = random_scores[random_scores > 0]
            print('  Done. Non-zero: {}/{}'.format(len(rnd_nz), len(random_scores)))
            for r in panel_results:
                s = r['score']
                if len(rnd_nz) > 0:
                    pct = 100.0 * (rnd_nz < s).sum() / len(rnd_nz)
                else:
                    pct = float('nan')
                r['percentile_vs_nonzero'] = round(pct, 1)
            out_data = {'method': 'Distinct cancer diseases reachable from compound targets via Gene-Disease edges', 'cancer_keywords': kws, 'gene_disease_edges': GENE_DISEASE_EDGES, 'n_cancer_diseases_in_spoke': n_cancer_total, 'n_random_total': len(random_scores), 'n_random_nonzero': int((random_scores > 0).sum()), 'random_nonzero_mean': round(float(rnd_nz.mean()), 1) if len(rnd_nz) > 0 else None, 'random_nonzero_median': round(float(np.median(rnd_nz)), 1) if len(rnd_nz) > 0 else None, 'panel': [{'compound': r['compound'], 'score': r['score'], 'detail': r['detail'], 'percentile_vs_nonzero': r['percentile_vs_nonzero']} for r in panel_results], 'random_scores': [int(x) for x in random_scores]}
            if args.out_json:
                os.makedirs(os.path.dirname(args.out_json) or '.', exist_ok=True)
                with open(args.out_json, 'w') as f:
                    json.dump(out_data, f, indent=2)
                print('Wrote', args.out_json)
            if args.out_fig:
                os.makedirs(os.path.dirname(args.out_fig) or '.', exist_ok=True)
                legend_path = args.out_legend or args.out_fig.replace('.png', '_legend.png')
                make_plots(panel_results, random_scores, args.out_fig, legend_path)
                print('Wrote', args.out_fig, legend_path)
    finally:
        driver.close()
    print('\n' + '=' * 72)
    print('CANCER GENE PROXIMITY SCORE (CGPS) — SUMMARY')
    print('=' * 72)
    print('Score = distinct cancer diseases reachable from target genes')
    print('Gene-Disease edges: {}'.format(', '.join(GENE_DISEASE_EDGES)))
    print('Cancer diseases in SPOKE: {}'.format(n_cancer_total))
    print('Random background: {} non-zero / {} total'.format(out_data['n_random_nonzero'], out_data['n_random_total']))
    if out_data['random_nonzero_mean'] is not None:
        print('Random non-zero mean: {:.1f}  median: {:.1f}'.format(out_data['random_nonzero_mean'], out_data['random_nonzero_median']))
    print()
    print('{:25s} {:>8s} {:>8s} {:>10s} {:>8s}'.format('Compound', 'Targets', 'CncGene', 'CncDis', 'Pctile'))
    print('-' * 72)
    for r in sorted(out_data['panel'], key=lambda x: x['score'], reverse=True):
        d = r['detail']
        print('{:25s} {:8d} {:8d} {:10d} {:7.1f}%'.format(r['compound']['name'], d['n_targets'], d['n_cancer_genes'], r['score'], r.get('percentile_vs_nonzero', 0)))
    print('=' * 72)
if __name__ == '__main__':
    main()
