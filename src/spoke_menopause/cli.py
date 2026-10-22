"""Explicit offline commands and opt-in graph/API operations."""
import argparse
import csv
import json
import os
from pathlib import Path
from .core import DWPC, enumerate_schemas, fingerprint, gene_counts, load_graph, personalized_pagerank


def repo_root():
    return Path(os.environ.get('SPOKE_RELEASE_ROOT',Path(__file__).resolve().parents[2]))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    for name in ['validate','export-tables','figures']:
        x=sub.add_parser(name);x.add_argument('--root',type=Path,default=repo_root());x.add_argument('--out',type=Path,default=Path('results')/name)
    x=sub.add_parser('dwpc');x.add_argument('--graph',required=True);x.add_argument('--source',required=True);x.add_argument('--target',required=True);x.add_argument('--max-length',type=int,default=3);x.add_argument('--out',type=Path,required=True)
    x=sub.add_parser('pagerank');x.add_argument('--graph',required=True);x.add_argument('--seeds',nargs='+',required=True);x.add_argument('--out',type=Path,required=True)
    x=sub.add_parser('extract-genes');x.add_argument('--graph',required=True);x.add_argument('--source',required=True);x.add_argument('--target',required=True);x.add_argument('--schemas',type=Path,required=True);x.add_argument('--cap',type=int,required=True);x.add_argument('--parallel-edges',action='store_true');x.add_argument('--out',type=Path,required=True)
    x=sub.add_parser('enrichr');x.add_argument('--genes',type=Path,required=True);x.add_argument('--libraries',nargs='+',required=True);x.add_argument('--out',type=Path,required=True)
    x=sub.add_parser('neo4j');x.add_argument('--query-file',type=Path,required=True);x.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.command in ['validate','export-tables']:
        from .paper import validate, export_tables
        if a.command=='validate':validate(a.root,a.out)
        else:export_tables(a.root,a.out)
    elif a.command=='figures':
        from .figures import render_panels
        render_panels(a.root,a.out)
    elif a.command in ['dwpc','pagerank','extract-genes']:
        g=load_graph(a.graph);a.out.parent.mkdir(parents=True,exist_ok=True)
        if a.command=='dwpc':
            engine=DWPC(g);schemas=enumerate_schemas(g,g.nodes[a.source]['node_type'],g.nodes[a.target]['node_type'],a.max_length)
            rows=[{'source':a.source,'target':a.target,'metapath':' → '.join(s),'length':len(s),'dwpc':engine.score(a.source,a.target,s),'path_count':engine.score(a.source,a.target,s,False)} for s in schemas]
        elif a.command=='pagerank':rows=[{'node':n,'pagerank':v} for n,v in personalized_pagerank(g,a.seeds).items()]
        else:
            schemas=json.loads(a.schemas.read_text());counts,rec,paths=gene_counts(g,a.source,a.target,schemas,a.cap,a.parallel_edges)
            rows=[{'gene_id':n,'gene_symbol':g.nodes[n].get('name',n),'total_occurrences':c,'n_metapaths':len(rec[n])} for n,c in counts.items()]
            a.out.with_suffix('.paths.json').write_text(json.dumps(paths))
        if rows:
            with a.out.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        a.out.with_suffix('.provenance.json').write_text(json.dumps({'command':a.command,'graph':fingerprint(a.graph),'parameters':{k:str(v) for k,v in vars(a).items() if k!='graph'}},indent=2))
    elif a.command=='enrichr':
        import requests
        genes=[s.strip() for s in a.genes.read_text().splitlines() if s.strip()]
        a.out.mkdir(parents=True,exist_ok=True)
        response=requests.post('https://maayanlab.cloud/Enrichr/addList',files={'list':(None,'\n'.join(genes)),'description':(None,'SPOKE menopause reproducibility run')},timeout=90);response.raise_for_status();uid=response.json()['userListId']
        for library in a.libraries:
            r=requests.get('https://maayanlab.cloud/Enrichr/enrich',params={'userListId':uid,'backgroundType':library},timeout=90);r.raise_for_status();(a.out/(library+'.json')).write_text(json.dumps(r.json()))
        (a.out/'input.json').write_text(json.dumps({'genes':genes,'libraries':a.libraries,'background':'Enrichr default','userListId':uid},indent=2))
    elif a.command=='neo4j':
        from neo4j import GraphDatabase, READ_ACCESS, Query
        queries=json.loads(a.query_file.read_text());a.out.mkdir(parents=True,exist_ok=True)
        with GraphDatabase.driver(os.environ['NEO4J_URI'],auth=(os.environ['NEO4J_USER'],os.environ['NEO4J_PASSWORD'])) as driver:
            with driver.session(database=os.environ.get('NEO4J_DATABASE'),default_access_mode=READ_ACCESS) as session:
                for name,q in queries.items():
                    rows=[dict(r) for r in session.run(Query(q['cypher'],timeout=300),**q.get('parameters',{}))]
                    (a.out/(name+'.json')).write_text(json.dumps(rows,indent=2,default=str))

if __name__=='__main__':main()
