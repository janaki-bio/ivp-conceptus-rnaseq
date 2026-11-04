"""report_leading_edge.py — Table 2 of the report: leading-edge genes behind the main gene-set patterns.
Run from the repository root after the R pipeline (needs pandas, numpy). Writes results/tables/leading_edge_drivers.csv"""
import pandas as pd, numpy as np, itertools, os, sys
T=(sys.argv[1] if len(sys.argv)>1 else '.')+'/results/tables/'
def de(f): return pd.read_csv(T+f'DE_{f}_all_genes.csv',dtype={'entrez':str}).drop_duplicates('entrez')
ED,ADJ,EEM=de('ED_IVP_vs_IVV'),de('ED_IVP_vs_IVV_discAdjusted'),de('TE_IVP_vs_IVV')
FULL={'ED':de('fullcohort_ED_full_sex_length'),'EEM':de('fullcohort_TE_full_sex_length')}
MALE={'ED':de('fullcohort_ED_sexmatched_males'),'EEM':de('fullcohort_TE_sexmatched_males')}
sym2e={**dict(zip(EEM.label,EEM.entrez)),**dict(zip(ED.label,ED.entrez))}
markers=set(pd.read_csv(T+'qc_disc_marker_genes.csv',dtype=str).entrez)
gsea={('GO','ED'):'GSEA_GO_BP_ED_IVP_vs_IVV',('GO','EEM'):'GSEA_GO_BP_TE_IVP_vs_IVV',('KEGG','ED'):'GSEA_KEGG_ED_IVP_vs_IVV',('KEGG','EEM'):'GSEA_KEGG_TE_IVP_vs_IVV'}
G={k:pd.read_csv(T+v+'.csv').set_index('Description') for k,v in gsea.items()}
# themes fixed in advance: representative term first, then related terms
themes=[('ED glycosylation','ED',[('GO','glycosylation'),('GO','protein glycosylation'),('GO','glycoprotein biosynthetic process'),('KEGG','N-Glycan biosynthesis')]),
 ('ED lipid metabolism','ED',[('GO','cellular lipid metabolic process'),('GO','lipid biosynthetic process'),('GO','membrane lipid metabolic process'),('GO','protein lipidation'),('KEGG','Glycerophospholipid metabolism')]),
 ('ED development / signalling','ED',[('GO','anatomical structure morphogenesis'),('GO','embryo development'),('GO','cell fate commitment'),('GO','canonical Wnt signaling pathway'),('KEGG','Wnt signaling pathway'),('KEGG','TGF-beta signaling pathway')]),
 ('EEM adhesion / ECM','EEM',[('GO','cell adhesion'),('GO','cell migration'),('GO','wound healing'),('KEGG','ECM-receptor interaction'),('KEGG','Integrin signaling')]),
 ('EEM translation / respiration','EEM',[('GO','cytoplasmic translation'),('GO','ribosome biogenesis'),('GO','cellular respiration'),('GO','oxidative phosphorylation'),('KEGG','Ribosome')]),
 ('Protein folding (ED and EEM)','both',[('GO','protein folding'),('GO','protein refolding'),("GO","'de novo' protein folding"),('GO','chaperone-mediated protein folding')])]
def le(db,tis,term):
    t=G[(db,tis)]
    if term not in t.index: return None
    g=t.loc[term,'core_enrichment'].split('/')
    return set(g) if db=='KEGG' else {sym2e.get(x,x) for x in g}
lab=dict(zip(ED.entrez,ED.label)); lab.update(dict(zip(EEM.entrez,EEM.label)))
rows=[]; drv=[]
for name,tis,terms in themes:
    tl=['ED','EEM'] if tis=='both' else [tis]
    for ti in tl:
        sets={f'{db}:{t}':le(db,ti,t) for db,t in terms}; sets={k:v for k,v in sets.items() if v}
        
        if not sets: continue
        rep=list(sets)[0]; R=sets[rep]
        jac=[len(a&b)/len(a|b) for a,b in itertools.combinations(sets.values(),2)]
        cnt=pd.Series([g for s in sets.values() for g in s]).value_counts()
        P=(ED if ti=='ED' else EEM).set_index('entrez')
        sgn=np.sign(P.loc[[g for g in R if g in P.index],'stat'])
        def agree(tab):
            tb=tab.set_index('entrez'); gg=[g for g in R if g in tb.index and g in P.index]
            return round(100*np.mean(np.sign(tb.loc[gg,'log2FoldChange'])==np.sign(P.loc[gg,'log2FoldChange'])),0)
        rec=cnt[cnt>=max(2,int(np.ceil(len(sets)/2)))].index
        recs=sorted(rec,key=lambda g:-abs(P['stat'].get(g,0)))[:8]
        rows.append(dict(theme=name+('' if tis!='both' else f' – {ti}'),tissue=ti,representative=rep.split(':',1)[1],related_sets=len(sets),
            le_size=len(R),union_size=len(set().union(*sets.values())),median_jaccard=round(np.median(jac),2) if jac else np.nan,
            recurrent=', '.join(lab.get(g,g) for g in recs),
            pct_disc_markers=round(100*len(R&markers)/len(R),0) if ti=='ED' else np.nan,
            same_dir_adj=agree(ADJ) if ti=='ED' else np.nan, same_dir_male=agree(MALE[ti]), same_dir_full=agree(FULL[ti])))
out=pd.DataFrame(rows); out.to_csv(T+'leading_edge_drivers.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',90); print(out.to_string())
