"""report_pathway_figures.py — Figure 5 (GSEA pathway overview) and Figure 10 (within-cohort sensitivity analyses) of the report.
Run from the repository root after the R pipeline (needs pandas, numpy, matplotlib)."""
import os
os.makedirs("report/figs", exist_ok=True)
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ---- Figure 5: pathway overview
T='results/tables/'
cols=[('ED','ED_IVP_vs_IVV'),('ED, disc-\nadjusted','ED_IVP_vs_IVV_discAdjusted'),('EEM','TE_IVP_vs_IVV'),('Interaction\n(EEM − ED)','interaction')]
groups=[
 ('GO: glycosylation and lipids',['glycosylation','protein lipidation','membrane lipid metabolic process','cellular lipid metabolic process'],'GO_BP'),
 ('GO: development and signalling',['anatomical structure morphogenesis','embryo development','cell fate commitment','neuron differentiation','canonical Wnt signaling pathway'],'GO_BP'),
 ('GO: adhesion, migration, tissue repair',['cell adhesion','cell migration','axon guidance','wound healing','blood coagulation'],'GO_BP'),
 ('GO: energy, translation, protein handling',['cellular respiration','oxidative phosphorylation','cytoplasmic translation','ribosome biogenesis','protein folding','endoplasmic reticulum to Golgi vesicle-mediated transport','protein N-linked glycosylation'],'GO_BP'),
 ('KEGG pathways',['Glycerophospholipid metabolism','N-Glycan biosynthesis','Sphingolipid metabolism','TGF-beta signaling pathway','Wnt signaling pathway','Notch signaling pathway','Signaling pathways regulating pluripotency of stem cells','ECM-receptor interaction','Integrin signaling','Ribosome','Protein processing in endoplasmic reticulum','Proteasome','Citrate cycle (TCA cycle)','Spliceosome'],'KEGG'),
]
tabs={(db,tag):pd.read_csv(f'{T}GSEA_{db}_{tag}.csv').set_index('Description') for db in ['GO_BP','KEGG'] for _,tag in cols}
SHORT={'endoplasmic reticulum to Golgi vesicle-mediated transport':'ER-to-Golgi vesicle-mediated transport','Signaling pathways regulating pluripotency of stem cells':'Pathways regulating stem-cell pluripotency','Protein processing in endoplasmic reticulum':'Protein processing in ER'}
rows=[];missing=[]
for g,terms,db in groups:
    rows.append(('H',g,None))
    for t in terms:
        vals=[]
        for _,tag in cols:
            d=tabs[(db,tag)]
            vals.append((d.loc[t,'NES'],d.loc[t,'p.adjust']) if t in d.index and d.loc[t,'p.adjust']<0.05 else (np.nan,np.nan))
        if all(np.isnan(v[0]) for v in vals): missing.append(t)
        rows.append(('T',t,vals))
print('never significant:',missing)
n=len(rows); fig,ax=plt.subplots(figsize=(8.6,0.27*n+1.3))
norm=TwoSlopeNorm(vmin=-3,vcenter=0,vmax=3); cmap=plt.get_cmap('RdBu_r')
ylabels=[];yticks=[]
for i,(k,lab,vals) in enumerate(rows):
    y=n-i
    if k=='H':
        ax.text(-0.62,y,lab,fontweight='bold',fontsize=9,ha='right' if False else 'left',va='center',transform=ax.get_yaxis_transform() if False else ax.transData, clip_on=False)
        continue
    yticks.append(y); ylabels.append(SHORT.get(lab,lab))
    for j,(nes,p) in enumerate(vals):
        if np.isnan(nes): ax.scatter(j,y,s=10,color='#c8c8c8',marker='o',lw=0)
        else:
            sz=40+60*min(-np.log10(p),8)/8*2
            ax.scatter(j,y,s=sz,color=cmap(norm(nes)),edgecolor='#333',lw=0.5)
            ax.text(j+0.2,y,f'{nes:+.1f}',fontsize=6.5,va='center',color='#444')
ax.set_yticks(yticks); ax.set_yticklabels(ylabels,fontsize=8)
ax.set_xticks(range(len(cols))); ax.set_xticklabels([c for c,_ in cols],fontsize=8.5)
ax.xaxis.tick_top(); ax.set_xlim(-0.5,len(cols)-0.3); ax.set_ylim(0.3,n+0.7)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
for j in range(len(cols)): ax.axvline(j,color='#eee',lw=0.8,zorder=0)
# headers placed in axis coords at left
for i,(k,lab,vals) in enumerate(rows):
    if k=='H': ax.annotate(lab,(0,n-i),xycoords=('axes fraction','data'),xytext=(-5,0),textcoords='offset points',ha='right',va='center',fontweight='bold',fontsize=8.5,color='#1f4e79')
for t in ax.texts[:]:
    if t.get_text() in [g for g,_,_ in groups] and t.get_position()[0]==-0.62: t.remove()
sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap); cb=fig.colorbar(sm,ax=ax,fraction=0.025,pad=0.02,shrink=0.35,anchor=(0,0.9))
cb.set_label('GSEA NES\n(IVP vs IVV)',fontsize=8); cb.ax.tick_params(labelsize=7)
plt.tight_layout(); fig.savefig('report/figs/fig5_pathways.png',dpi=200,bbox_inches='tight')

# ---- Figure 8: robustness checks
d=pd.read_csv('results/tables/fullcohort_key_terms_GSEA.csv')
norm=TwoSlopeNorm(vmin=-3,vcenter=0,vmax=3); cmap=plt.get_cmap('RdBu_r')
panels=[('ED',['ED_primary','ED_sexmatched','ED_full','ED_primary_disc','ED_full_disc'],['Primary','Sex-matched\n(males)','Full\ncohort','Primary\n+ disc adj.','Full cohort\n+ disc adj.'],
  ['glycosylation','cellular lipid metabolic process','anatomical structure morphogenesis','neuron differentiation','cell adhesion']),
 ('EEM',['TE_primary','TE_sexmatched','TE_full'],['Primary','Sex-matched\n(males)','Full\ncohort'],
  ['translation','ribosome biogenesis','cellular respiration','protein folding','cell migration','cell adhesion'])]
fig,axs=plt.subplots(2,1,figsize=(7.2,7.4),gridspec_kw={'height_ratios':[5,6]})
for ax,(tis,an,lab,terms) in zip(axs,panels):
    for j,a in enumerate(an):
        for i,t in enumerate(terms):
            r=d[(d.analysis==a)&(d.term==t)]
            y=len(terms)-i
            if r.empty: ax.scatter(j,y,s=10,color='#ccc'); continue
            nes,p=r.NES.iloc[0],r['p.adjust'].iloc[0]; s=p<0.05
            ax.scatter(j,y,s=230,color=cmap(norm(nes)),edgecolor='#222' if s else '#bbb',lw=1.4 if s else 0.8,alpha=1 if s else 0.55)
            ax.text(j,y,f'{nes:+.1f}',ha='center',va='center',fontsize=6.5,color='white' if abs(nes)>1.8 and s else '#222')
    ax.set_yticks(range(len(terms),0,-1)); ax.set_yticklabels(terms,fontsize=8.5)
    ax.set_xticks(range(len(an))); ax.set_xticklabels(lab,fontsize=7.8); ax.xaxis.tick_top()
    ax.set_xlim(-0.6,4.6); ax.set_ylim(0.4,len(terms)+0.6)
    ax.set_title(tis,loc='left',fontweight='bold',color='#1f4e79',fontsize=10,pad=30)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    if tis=='ED': ax.axvline(2.5,color='#999',lw=0.8,ls='--')
sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap); cb=fig.colorbar(sm,ax=axs,fraction=0.02,pad=0.02); cb.set_label('GSEA NES (IVP vs IVV)',fontsize=8); cb.ax.tick_params(labelsize=7)
fig.subplots_adjust(hspace=0.45)
fig.savefig('report/figs/fig10_cohort_sensitivity.png',dpi=200,bbox_inches='tight')
