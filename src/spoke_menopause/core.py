"""Graph statistics matching the manuscript's stated methods.

DWPC uses binary, predicate-specific adjacency and degree weighting at both
ends of every edge. It counts walks, including repeated nodes. Gene extraction
is separate and can preserve the historical Q1 parallel-edge multiplicity.
"""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import pickle
import numpy as np
import networkx as nx
from scipy import sparse


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if not len(p):
        return p.copy()
    if np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("P values must be finite and between zero and one")
    order = np.argsort(p, kind="stable")
    ranked = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty(len(p)); result[order] = np.minimum(ranked, 1)
    return result


def load_graph(path):
    """Load a trusted NetworkX pickle; never load an untrusted pickle."""
    with Path(path).open("rb") as f:
        graph = pickle.load(f)
    if not isinstance(graph, nx.Graph):
        raise TypeError("Expected a NetworkX graph")
    return graph


def fingerprint(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return {"filename": Path(path).name, "sha256": h.hexdigest(), "bytes": Path(path).stat().st_size}


class DWPC:
    def __init__(self, graph, damping=-0.4):
        if graph.is_directed():
            raise ValueError("The manuscript DWPC input is an undirected MultiGraph")
        self.graph = graph
        self.index = {n:i for i,n in enumerate(graph)}
        self.damping = damping
        self.adjacency = {}; self.weighted = {}

    def matrix(self, predicate, weighted=True):
        if predicate not in self.adjacency:
            row=[]; col=[]
            for u,v,d in self.graph.edges(data=True):
                if d.get("predicate") == predicate:
                    row.extend((self.index[u],self.index[v])); col.extend((self.index[v],self.index[u]))
            a=sparse.csr_matrix((np.ones(len(row)),(row,col)), shape=(len(self.index),len(self.index)))
            a.sum_duplicates(); a.data[:]=1.0
            self.adjacency[predicate]=a
            degree=np.asarray(a.sum(axis=1)).ravel();degree[degree==0]=1
            d=sparse.diags(degree**self.damping)
            self.weighted[predicate]=d @ a @ d
        return self.weighted[predicate] if weighted else self.adjacency[predicate]

    def score(self, source, target, schema, weighted=True):
        if not schema:
            raise ValueError("A schema must contain at least one predicate")
        if source not in self.index or target not in self.index:
            raise KeyError("Source or target is absent from the graph")
        row=sparse.csr_matrix(([1.0],([0],[self.index[source]])),shape=(1,len(self.index)))
        for predicate in schema:
            row=row @ self.matrix(predicate,weighted)
        value=float(row[0,self.index[target]])
        return value if weighted else int(round(value))


def enumerate_schemas(graph, source_type, target_type, max_length=3):
    transitions=defaultdict(set)
    for u,v,d in graph.edges(data=True):
        a,b=graph.nodes[u]['node_type'],graph.nodes[v]['node_type']
        transitions[a].add((d['predicate'],b))
        if not graph.is_directed(): transitions[b].add((d['predicate'],a))
    result=[]
    def walk(node_type,path):
        if path and node_type==target_type: result.append(tuple(path))
        if len(path)==max_length:return
        for predicate,next_type in sorted(transitions[node_type]):walk(next_type,path+[predicate])
    walk(source_type,[])
    return sorted(set(result),key=lambda x:(len(x),x))


def extract_paths(graph, source, target, schema, cap, parallel_edges=False):
    """Ordered DFS with a per-schema cap, retaining repeated nodes.

    parallel_edges=True reproduces the Q1 v2 extractor's traversal of every
    matching multiedge. False reproduces Q2's one visit per matching neighbor.
    Neighbor insertion order affects capped outputs and is preserved.
    """
    if cap < 1: raise ValueError("cap must be positive")
    paths=[]
    def step(node,i,path):
        if len(paths)>=cap:return
        if i==len(schema):
            if node==target:paths.append(path.copy())
            return
        for neighbor in graph.neighbors(node):
            edges=graph.get_edge_data(node,neighbor)
            values=edges.values() if graph.is_multigraph() else [edges]
            matches=sum(e.get('predicate')==schema[i] for e in values)
            for _ in range(matches if parallel_edges else int(matches>0)):
                step(neighbor,i+1,path+[neighbor])
                if len(paths)>=cap:return
    step(source,0,[source])
    return paths


def gene_counts(graph, source, target, schemas, cap=50, parallel_edges=False):
    occurrences=Counter(); recurrence=defaultdict(set); all_paths=[]
    for i,schema in enumerate(schemas):
        paths=extract_paths(graph,source,target,schema,cap,parallel_edges)
        for path in paths:
            all_paths.append({'schema_index':i,'nodes':path})
            for node in path[1:-1]:
                if graph.nodes[node].get('node_type')=='Gene':
                    occurrences[node]+=1;recurrence[node].add(i)
    return occurrences,recurrence,all_paths


def personalized_pagerank(graph,seeds):
    seeds=list(dict.fromkeys(seeds))
    if not seeds or any(s not in graph for s in seeds):
        raise ValueError("All seed identifiers must be present in the graph")
    return nx.pagerank(graph,alpha=.85,personalization={s:1/len(seeds) for s in seeds},max_iter=500,tol=1e-20)
