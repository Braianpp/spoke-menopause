import ast
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple
import networkx as nx
import numpy as np
import pytest
from scipy import sparse
from spoke_menopause.core import DWPC, bh_adjust, enumerate_schemas, extract_paths, personalized_pagerank

ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture
def graph():
    g=nx.MultiGraph()
    for n,t in [('s','Disease'),('a','Gene'),('b','Gene'),('t','Disease')]:g.add_node(n,node_type=t)
    for u,v,p in [('s','a','DaG'),('s','b','DaG'),('a','t','DaG'),('b','t','DaG'),('a','b','GrG'),('a','a','GrG'),('s','a','DaG')]:g.add_edge(u,v,predicate=p)
    return g

def test_bh_known_values_and_invalid():
    assert np.allclose(bh_adjust([.04,.001,.03,.02]),[.04,.004,.04,.04])
    assert len(bh_adjust([]))==0
    with pytest.raises(ValueError):bh_adjust([float('nan')])

def original_classes():
    env=dict(nx=nx,np=np,sparse=sparse,defaultdict=defaultdict,List=List,Dict=Dict,Set=Set,Tuple=Tuple)
    for name in ['dwpc_016.py','dwpc_014.py']:
        tree=ast.parse((ROOT/'analysis/shared/notebook_cells'/name).read_text())
        exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.ClassDef)],type_ignores=[]),name,'exec'),env)
    return env

def test_matches_original_matrix_products(graph):
    original=original_classes()['NetworkDWPC'](graph);new=DWPC(graph)
    for length in [1,2,3,4]:
        import itertools
        for schema in itertools.product(['DaG','GrG'],repeat=length):
            for source in graph:
                for target in graph:
                    assert new.score(source,target,schema)==pytest.approx(original.compute_dwpc(source,target,list(schema)))
                    assert new.score(source,target,schema,False)==original.compute_simple_path_count(source,target,list(schema))

def test_schema_enumeration_matches_original(graph):
    old=original_classes()['MetapathFinder'](graph)
    for source in ['Disease','Gene']:
        for target in ['Disease','Gene']:
            assert set(enumerate_schemas(graph,source,target,3))=={tuple(x) for x in old.find_metapaths(source,target,3)}

def test_weighted_walk_manual(graph):
    # Two 2-edge walks, each weighted by 2^-0.4 at all four endpoints.
    assert DWPC(graph).score('s','t',['DaG','DaG'])==pytest.approx(2*2**(-1.6))
    # A repeated node is permitted: s-a-s and s-b-s.
    assert DWPC(graph).score('s','s',['DaG','DaG'],False)==2

def test_parallel_edges_and_caps(graph):
    assert extract_paths(graph,'s','t',['DaG','DaG'],50)==[['s','a','t'],['s','b','t']]
    assert extract_paths(graph,'s','t',['DaG','DaG'],50,True)==[['s','a','t'],['s','a','t'],['s','b','t']]
    assert extract_paths(graph,'s','t',['DaG','DaG'],1)==[['s','a','t']]

def test_pagerank_stationary_distribution(graph):
    p=personalized_pagerank(graph,['s']);assert sum(p.values())==pytest.approx(1)
    # NetworkX adjacency counts self-loops once in its transition matrix.
    nodes=list(graph);a=nx.to_numpy_array(graph,nodelist=nodes);a=a/a.sum(axis=1,keepdims=True)
    vector=np.array([p[n] for n in nodes]);seed=np.array([float(n=='s') for n in nodes])
    assert np.allclose(vector,.85*vector@a+.15*seed,atol=1e-12)

def test_reject_directed_dwpc():
    with pytest.raises(ValueError):DWPC(nx.DiGraph())

def test_q6_reverse_pruning_preserves_order_and_cap(graph):
    tree=ast.parse((ROOT/'analysis/q06/rerun_matched_null.py').read_text())
    env={'G':graph,'CAP':2,'cache':{}}
    names={'adj','intersect_ordered','paths'}
    exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names],type_ignores=[]),'q6_functions','exec'),env)
    import itertools
    for length in [1,2,3]:
        for schema in itertools.product(['DaG','GrG'],repeat=length):
            for source in graph:
                for target in graph:
                    assert env['paths'](source,target,schema)==extract_paths(graph,source,target,schema,2)
