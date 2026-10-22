"""Portable numerical panels; layout is not the manually assembled manuscript."""
import json
import re
import textwrap
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from .paper import load_tables, records

TEAL='#98bdb7'
DISEASE_COLORS={'Vascular dementia':'#8062aa','Frontotemporal dementia':'#3498db',"Parkinson's disease":'#2ecc71'}


def render_panels(root, out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    tables=load_tables(root); manifest=[]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    def rows(n,b=0): return records(tables,n,b)
    def save(fig,name,source):
        fig.savefig(out/(name+'.png'),dpi=300,bbox_inches='tight',facecolor='white')
        fig.savefig(out/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)
        manifest.append({'panel':name,'source_tables':source})
    def bars(name, labels, values, xlabel, title, source, colors=TEAL, annotations=None, errors=None, legend=None):
        fig,ax=plt.subplots(figsize=(10,max(3.2,len(labels)*.33+1.3)),layout='constrained')
        y=np.arange(len(labels));ax.barh(y,values,color=colors,edgecolor='#6e8791',linewidth=.6)
        ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlabel(xlabel);ax.set_title(title,pad=15)
        ax.set_axisbelow(True);ax.grid(axis='x',color='#e6e9eb')
        if annotations:
            for yi,v,a in zip(y,values,annotations):ax.annotate(a,(v,yi),xytext=(4,0),textcoords='offset points',va='center',fontsize=9)
            ax.set_xlim(0,max(values)*1.18)
        if errors:
            mean,sd=map(np.array,errors);ax.errorbar(mean,y,xerr=np.vstack([np.minimum(mean,sd),sd]),fmt='o',color='black',markersize=3,capsize=2,label='Null mean ± SD (clipped at 0)');ax.legend(loc='lower right',fontsize=9)
        if legend:ax.legend(handles=[Patch(facecolor=c,label=l) for l,c in legend],loc='lower right',fontsize=9)
        save(fig,name,source)
    def enrichment(name,rs,title,source):
        labels=['\n'.join(textwrap.wrap(re.sub(r' \(GO:\d+\)$','',r['Term']),48)) for r in rs]
        colors=[{'KEGG_2021_Human':'#8e7ba4','Reactome_2022':'#66a3ba'}.get(r['Library'],TEAL) for r in rs]
        libraries={'GO_Biological_Process_2023':('GO Biological Process',TEAL),'KEGG_2021_Human':('KEGG','#8e7ba4'),'Reactome_2022':('Reactome','#66a3ba')}
        legend=[libraries[k] for k in dict.fromkeys(r['Library'] for r in rs) if k in libraries]
        bars(name,labels,[-np.log10(r['Adjusted P value']) for r in rs],r'$-\log_{10}$(adjusted P value)',title,source,colors,[f"n={r['Overlapping genes (n)']}" for r in rs],legend=legend if len(legend)>1 else None)
    r=[r for r in rows(4) if r['Shown in Figure 2A']]
    r.sort(key=lambda x:x['DWPC'],reverse=True)
    bars('Figure2A',[r['Metapath'] for r in r],[r['DWPC'] for r in r],'DWPC','Premature menopause–dementia metapaths',[4])
    r=[r for r in rows(5) if r['Shown in Figure 2C']]
    r.sort(key=lambda x:x['Observed recurrence'],reverse=True)
    bars('Figure2C',[x['Gene symbol'] for x in r],[x['Observed recurrence'] for x in r],'Gene occurrences across extracted paths','Intermediate-gene recurrence',[5],errors=([x['Null mean'] for x in r],[x['Null SD'] for x in r]))
    enrichment('Figure2D',[r for r in rows(7) if r['Shown in Figure 2D']],'Biological processes enriched in the 55-gene Q1 set',[7])
    enrichment('Figure3B',[r for r in rows(11) if r['Shown in Figure 3B']],'Biological processes enriched in 27 shared intermediates',[11])
    r=sorted(rows(9),key=lambda x:list(DISEASE_COLORS).index(x['Target disease']))
    bars('Figure3D',[x['Target disease'] for x in r],[x['Cumulative DWPC'] for x in r],'Cumulative DWPC','Metapath connectivity by target disease',[9],[DISEASE_COLORS[x['Target disease']] for x in r],[f"{x['Cumulative DWPC']:.3f}" for x in r])
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    r=[r for r in rows(13,1) if r['PageRank']>0];ax.hist([np.log10(x['PageRank']+1e-10) for x in r],bins=50,color='#d8dadd',edgecolor='white')
    palette=plt.get_cmap('tab10')
    for i,r in enumerate(rows(13)):ax.axvline(np.log10(r['PageRank']+1e-10),color=palette(i),ls='--',lw=2,label=r['Symptom'])
    ax.set_xlabel(r'$\log_{10}$(PageRank + $10^{-10}$)');ax.set_ylabel('Symptom count');ax.legend(fontsize=9);ax.set_title('Dementia proximity of perimenopausal symptoms');save(fig,'Figure4A',[13])
    enrichment('Figure4C',[r for r in rows(15,1) if r['Shown in Figure 4C']],'Enrichment of 15 symptom–dementia route genes',[15])
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained');r=rows(19)
    values=[np.log10(x['PageRank']+1e-10) for x in r];ax.hist(values,bins=60,density=True,color='#cdd6da',edgecolor='white')
    for identifier,label,color in [('DOID:11476','Osteoporosis','#c84952'),('DOID:1307','Dementia','#f28e2b')]:
        value=next(x['PageRank'] for x in r if x['Disease identifier']==identifier);ax.axvline(np.log10(value+1e-10),label=label,color=color,ls='--',lw=2)
    ax.set_xlabel(r'$\log_{10}$(PageRank + $10^{-10}$)');ax.set_ylabel('Density');ax.set_title('Disease proximity to premature menopause');ax.legend();save(fig,'Figure5A',[19])
    r=sorted([r for r in rows(18) if r['Path count']>0],key=lambda x:x['DWPC'],reverse=True)[:11]
    bars('Figure5B',[x['Metapath'] for x in r],[x['DWPC'] for x in r],'DWPC','Osteoporosis–dementia metapaths',[18])
    r=rows(28);bars('Figure6_CGPS',[x['Compound'] for x in r],[x['Distinct cancer diseases (CGPS)'] for x in r],'Distinct connected cancer diseases (CGPS)','Cancer-disease connectivity of selected compounds',[28])
    r=rows(32);bars('Figure6_PPR',[x['Compound'] for x in r],[x['PageRank'] for x in r],'Dementia-seeded PageRank','Proximity of selected breast-cancer therapies',[32])
    for length in [3,4]:
        r=sorted([r for r in rows(36) if r['Number of edges']==length and r['DWPC']>0],key=lambda x:x['DWPC'],reverse=True)[:10]
        bars(f'Figure8_{length}edge',[x['Metapath'] for x in r],[x['DWPC'] for x in r],'DWPC',f'Letrozole–dementia: top {length}-edge metapaths',[36])
    r=rows(39,1);columns=list(r[0])[1:];matrix=np.array([[x[k] for k in columns] for x in r])
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained');ax.imshow(matrix,cmap='Blues',vmin=0,vmax=6,aspect='auto')
    ax.set_xticks(range(len(columns)),columns,rotation=35,ha='right');ax.set_yticks(range(len(r)),[x['Gene symbol'] for x in r])
    labels={0:'',1:'B',2:'I',3:'B / I',5:'D',6:'I / U / D'}
    for i in range(len(r)):
        for j in range(len(columns)):ax.text(j,i,labels.get(int(matrix[i,j]),str(matrix[i,j])),ha='center',va='center',color='white' if matrix[i,j]>=4 else 'black')
    ax.set_title('Hormone ingredient–pharmacogene relationships\nB: binds; I: interacts; U/D: up/downregulates',fontsize=12);save(fig,'Q10_pharmacogenes',[39,40])
    (out/'index.json').write_text(json.dumps({'scope':'Numerical reconstruction from published tables; typography and layout differ from manual manuscript composites. Figure5A uses generic dementia DOID:1307.','panels':manifest},indent=2))
    print(f'Rendered {len(manifest)} numerical panels to {out}')
