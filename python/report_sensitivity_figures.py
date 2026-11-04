"""report_sensitivity_figures.py — Figure 11 of the report (embryo-influence heatmap). Run from the repository root
after R/10_sensitivity.R (needs pandas, numpy, matplotlib). Reads results/tables/sens_loeo_terms.csv."""
import os, sys
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
t = pd.read_csv(os.path.join(ROOT, 'results/tables/sens_loeo_terms.csv'))
os.makedirs('report/figs', exist_ok=True)
terms = ['glycosylation', 'cellular lipid metabolic process', 'lipid biosynthetic process', 'anatomical structure morphogenesis',
         'embryo development', 'neuron differentiation', 'canonical Wnt signaling pathway', 'cell adhesion', 'cell migration',
         'translation', 'ribosome biogenesis', 'cellular respiration', 'protein folding']
cols = ['none', 'IVP9', 'IVP10', 'IVP11', 'IVP14', 'VE24', 'VE26', 'VE28', 'VE29', 'IVP9_T1A library']
lab = {'none': 'All 8', 'IVP9_T1A library': 'IVP9_T1A\nonly'}
panels = [('ED', 'ED'), ('ED_adj', 'ED, disc-score-adjusted'), ('EEM', 'EEM'), ('interaction', 'Interaction (EEM − ED)')]
norm = TwoSlopeNorm(vmin=-2.8, vcenter=0, vmax=2.8); cmap = plt.get_cmap('RdBu_r')
fig, axs = plt.subplots(2, 2, figsize=(13, 9.6)); fig.subplots_adjust(hspace=0.38, wspace=0.08)
for k, (ax, (a, title)) in enumerate(zip(axs.flat, panels)):
    x = t[t.analysis == a].pivot_table(index='term', columns='omitted', values='NES')
    x = x.reindex(index=terms, columns=[c for c in cols if c in x.columns])
    p = t[t.analysis == a].pivot_table(index='term', columns='omitted', values='padj').reindex(index=x.index, columns=x.columns)
    im = ax.imshow(x.values, cmap=cmap, norm=norm, aspect='auto')
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            v = x.iat[i, j]
            if np.isnan(v): continue
            ax.text(j, i, f'{v:+.1f}' + ('*' if p.iat[i, j] < 0.05 else ''), ha='center', va='center', fontsize=6.6,
                    color='white' if abs(v) > 1.9 else '#222')
    ax.axvline(0.5, color='black', lw=1.2)
    ax.set_xticks(range(x.shape[1])); ax.set_xticklabels([lab.get(c, c) for c in x.columns], fontsize=7.5, rotation=0)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms if k % 2 == 0 else [''] * len(terms), fontsize=8)
    ax.set_title(title, loc='left', fontsize=10, fontweight='bold', color='#1f4e79', pad=24)
    ax.tick_params(length=0)
    ax.text(0.5 / x.shape[1] + 0.5 * (1 - 1 / x.shape[1]), 1.10, 'embryo omitted', transform=ax.transAxes, ha='center', fontsize=7.5, color='#555')
cb = fig.colorbar(im, ax=axs, fraction=0.02, pad=0.02); cb.set_label('NES (IVP vs IVV)', fontsize=8)
fig.savefig('report/figs/fig11_embryo_influence.png', dpi=200, bbox_inches='tight')
print('done')
