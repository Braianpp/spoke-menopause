# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import os, sys, time, pickle
from pathlib import Path
from importlib.machinery import SourceFileLoader
import pandas as pd
try:
    import networkx as nx
except Exception:
    nx = None
BASE = Path('${SPOKE_DATA_ROOT}')
LOGDIR = BASE / 'biorouter' / 'logs'
LOGDIR.mkdir(parents=True, exist_ok=True)
CANDS = [BASE / 'spoke_graph.gpickle', Path('${SPOKE_DATA_ROOT}/V1_homosapiens_1N_9E.gpickle'), Path('${SPOKE_DATA_ROOT}/extract_spoke/spokeML/extract_spoke/V1_homosapiens_1N_9E.gpickle')]
G = None
for gfile in CANDS:
    try:
        if nx is not None and hasattr(nx, 'read_gpickle'):
            G = nx.read_gpickle(gfile)
        else:
            raise AttributeError('no read_gpickle')
    except Exception:
        G = None
        if gfile.exists():
            try:
                with open(gfile, 'rb') as fh:
                    G = pickle.load(fh)
            except Exception:
                G = None
    if G is not None:
        print(f'[q1_run_null_v2] Loaded graph from {gfile}')
        break
if G is None:
    print('[q1_run_null_v2] ERROR: Could not load G from known paths', file=sys.stderr)
    sys.exit(2)
gene_csv = BASE / 'pm_dementia_dwpc_recurring_genes_v2.csv'
mp_csv = BASE / 'df_premature_menopause_to_dementia.csv'
if not gene_csv.exists() or not mp_csv.exists():
    print('[q1_run_null_v2] ERROR: Missing CSVs', file=sys.stderr)
    sys.exit(3)
gene_summary = pd.read_csv(gene_csv)
df_metapaths = pd.read_csv(mp_csv)
mod = SourceFileLoader('q1v2', str(BASE / 'biorouter' / 'scripts' / 'q1_null_recompute_v2.py')).load_module()
os.environ.setdefault('N_RANDOM_PAIRS', '1500')
os.environ.setdefault('MAX_PATHS', '50')
os.environ.setdefault('TOP_METAPATHS_LIMIT', '11')
os.environ.setdefault('ALPHA', '0.05')
print('[q1_run_null_v2] Settings:', {k: os.environ[k] for k in ['N_RANDOM_PAIRS', 'MAX_PATHS', 'TOP_METAPATHS_LIMIT', 'ALPHA']})
start = time.time()
null_df = mod.run_null(G=G, gene_summary=gene_summary, df_metapaths=df_metapaths)
elapsed = time.time() - start
print(f'[q1_run_null_v2] Completed null recompute in {elapsed / 60:.1f} min; rows={len(null_df)}')
