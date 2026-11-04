"""report_signalling_pca.py — (1) box plots of Wnt, BMP/TGF-beta and Notch genes; (2) correlation of principal
components with sample metadata. Run from the repository root (needs pandas, numpy, matplotlib).
Usage: python3 python/report_signalling_pca.py [repo_root] [count_matrix]
Inputs: data/D15_Embryos_RawCountData.txt, metadata/samples.csv, results/tables/*.csv
Outputs: report/figs/fig8_signalling.png, report/figs/figS3_pc_metadata.png, results/tables/pc_metadata_correlation.csv,
         results/tables/signalling_genes.csv"""
import os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math


def _betacf(a, b, x, itmax=300, eps=3e-14):
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1 / (d if abs(d) > 1e-300 else 1e-300); h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c; c = c if abs(c) > 1e-300 else 1e-300; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c; c = c if abs(c) > 1e-300 else 1e-300
        de = d * c; h *= de
        if abs(de - 1) < eps: break
    return h


def _betai(a, b, x):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
    return bt * _betacf(a, b, x) / a if x < (a + 1) / (a + b + 2) else 1 - bt * _betacf(b, a, 1 - x) / b


def pearsonr(x, y):
    """Pearson r and two-sided p-value (t distribution, n-2 df); same result as scipy.stats.pearsonr."""
    x, y = np.asarray(x, float), np.asarray(y, float); n = len(x)
    r = float(np.corrcoef(x, y)[0, 1]); r = max(min(r, 1.0), -1.0)
    df = n - 2
    if abs(r) == 1: return r, 0.0
    t2 = r * r * df / (1 - r * r)
    return r, _betai(df / 2, 0.5, df / (df + t2))

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
COUNTS = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'data/D15_Embryos_RawCountData.txt')
T = os.path.join(ROOT, 'results/tables/'); FIG = os.path.join(ROOT, 'report/figs/')
os.makedirs(FIG, exist_ok=True)
EXCLUDE = ['IVPF10_EC']                       # IVP disc without matching EEM library (primary analysis)
BLUE, RED = '#2a78d6', '#d62728'

# ---------------- data ----------------
cts = pd.read_csv(COUNTS, sep='\t', index_col=0)
cts.columns = [c.replace('_1', '') if c.endswith('_1') else c for c in cts.columns]
cts.index = [i.replace('GeneID:', '') for i in cts.index]
meta = pd.read_csv(os.path.join(ROOT, 'metadata/samples.csv')).set_index('sample').loc[cts.columns]
sex = pd.read_csv(T + 'qc_sex_inference_XIST.csv').set_index('sample')['sex_inferred']
lib = pd.read_csv(T + 'qc_library_sizes.csv').set_index('sample')['reads_M']
disc = pd.read_csv(T + 'qc_disc_score.csv').set_index('sample')['disc_score']
meta = meta.assign(sex=sex, reads_M=lib)

def norm_log(c):
    c = c[(c >= 10).sum(axis=1) >= 4]
    lg = np.log(c.replace(0, np.nan)); geo = lg.mean(axis=1)
    ok = geo.notna() & np.isfinite(geo)
    sf = np.exp((lg[ok].sub(geo[ok], axis=0)).median())
    return np.log2(c / sf + 1)

# ---------------- (1) signalling genes ----------------
sets = {'Wnt': ['WNT3', 'WNT5A', 'WNT10B', 'LEF1', 'NKD1', 'FRAT2'],
        'BMP / TGF-β': ['BMP4', 'ID1', 'ID4', 'SMAD7', 'SMAD9', 'TGFB2'],
        'Notch': ['NOTCH1', 'DLL3', 'HEY1', 'MFNG', 'DTX1']}
de = {k: pd.read_csv(T + f'DE_{f}_all_genes.csv', dtype={'entrez': str}).drop_duplicates('label').set_index('label')
      for k, f in [('ED', 'ED_IVP_vs_IVV'), ('ED_adj', 'ED_IVP_vs_IVV_discAdjusted'), ('EEM', 'TE_IVP_vs_IVV'),
                   ('ED_males', 'fullcohort_ED_sexmatched_males'), ('ED_full', 'fullcohort_ED_full_sex_length'),
                   ('ED_full_disc', 'fullcohort_ED_full_sex_length_disc'),
                   ('EEM_males', 'fullcohort_TE_sexmatched_males'), ('EEM_full', 'fullcohort_TE_full_sex_length')]}
markers = set(pd.read_csv(T + 'qc_disc_marker_genes.csv', dtype=str)['label'])
sl = meta[meta.same_length & ~meta.index.isin(EXCLUDE)]
lx = norm_log(cts[sl.index])
groups = [('ED', 'IVV'), ('ED', 'IVP'), ('TE', 'IVV'), ('TE', 'IVP')]

def q(tab, g):
    v = de[tab]['padj'].get(g, np.nan)
    return 'n.t.' if pd.isna(v) else (f'{v:.0e}' if v < 0.001 else f'{v:.3f}')

rows = []
ncol = 6
fig, axs = plt.subplots(len(sets), ncol, figsize=(15, 3.35 * len(sets)), squeeze=False)
for r, (pw, genes) in enumerate(sets.items()):
    for c in range(ncol):
        ax = axs[r][c]
        if c >= len(genes): ax.axis('off'); continue
        g = genes[c]; e = de['ED'].loc[g, 'entrez']
        y = lx.loc[e] if e in lx.index else None
        data = [y[sl[(sl.region == t) & (sl.type == ty)].index].values for t, ty in groups]
        cols = [BLUE, RED, BLUE, RED]
        bp = ax.boxplot(data, positions=range(4), widths=0.55, patch_artist=True, showfliers=False)
        for k, b in enumerate(bp['boxes']): b.set(facecolor='white', edgecolor=cols[k], lw=1.3)
        for k, m in enumerate(bp['medians']): m.set(color=cols[k], lw=2)
        for k in range(4):
            for part in ('whiskers', 'caps'):
                for w in bp[part][2 * k:2 * k + 2]: w.set(color=cols[k], lw=1)
            ax.scatter(np.full(len(data[k]), k) + np.linspace(-0.1, 0.1, len(data[k])), data[k], s=14, color=cols[k], zorder=3)
        ax.axvline(1.5, color='#ddd', lw=1)
        ax.set_xticks(range(4)); ax.set_xticklabels(['ED\nIVV', 'ED\nIVP', 'EEM\nIVV', 'EEM\nIVP'], fontsize=7.5)
        star = '†' if g in markers else ''
        ax.set_title(f'{g}{star}\nED {q("ED", g)} · adj {q("ED_adj", g)} · EEM {q("EEM", g)}', fontsize=7.6)
        ax.tick_params(labelsize=7); [ax.spines[s].set_visible(False) for s in ('top', 'right')]
        ax.grid(axis='y', color='#eee', lw=0.6)
        if c == 0: ax.set_ylabel(f'{pw}\nlog2 normalised count', fontsize=8.5, fontweight='bold')
        row = dict(pathway=pw, gene=g, disc_marker=g in markers)
        for k in de:
            row[f'{k}_log2FC'] = de[k]['log2FoldChange'].get(g, np.nan); row[f'{k}_padj'] = de[k]['padj'].get(g, np.nan)
        rows.append(row)
plt.tight_layout(h_pad=2.2)
fig.savefig(FIG + 'fig8_signalling.png', dpi=200, bbox_inches='tight'); plt.close(fig)
pd.DataFrame(rows).to_csv(T + 'signalling_genes.csv', index=False)

# ---------------- (2) PC–metadata correlation ----------------
# disc-specific transcript score for every ED library: mean log2 normalised expression of the disc marker genes
_ed = meta[meta.region == 'ED'].index
_x = norm_log(cts[_ed]); _mk = [m for m in pd.read_csv(T + 'qc_disc_marker_genes.csv', dtype=str)['entrez'] if m in _x.index]
meta['disc_score'] = np.nan; meta.loc[_ed, 'disc_score'] = _x.loc[_mk].mean().values

def pc_corr(samples, variables, npc=10):
    x = norm_log(cts[samples.index])
    top = x.var(axis=1).sort_values(ascending=False).index[:500]
    m = x.loc[top].T; m = m - m.mean()
    u, s, vt = np.linalg.svd(m.values, full_matrices=False)
    pcs = u * s; ve = s ** 2 / (s ** 2).sum()
    npc = min(npc, len(samples) - 1)
    R = pd.DataFrame(index=list(variables), columns=range(npc), dtype=float); P = R.copy()
    for name, v in variables.items():
        for k in range(npc):
            if np.std(v) == 0: continue
            R.loc[name, k], P.loc[name, k] = pearsonr(v, pcs[:, k])
    return R, P, ve[:npc]

def b(series, one): return (series == one).astype(float).values
def vars_for(s, disc_score=False):
    v = {'Tissue (EEM)': b(s.region, 'TE'), 'Production (IVP)': b(s.type, 'IVP'), 'Sex (male)': b(s.sex, 'male'),
         'Batch label ("Vitrified")': b(s.cryopreserved, 'Vitrified'), 'Same-length set': s.same_length.astype(float).values,
         'Section (C vs A)': b(s.section, 'C'), 'Library size': s.reads_M.values}
    if disc_score: v['Disc-specific score'] = s.disc_score.values
    return {k: x for k, x in v.items() if np.std(x) > 0}

panels = [('A', 'All libraries (n = 32)', meta),
          ('B', 'Embryonic disc libraries (n = 16)', meta[meta.region == 'ED']),
          ('C', 'Extra-embryonic membrane libraries (n = 16)', meta[meta.region == 'TE'])]
fig, axs = plt.subplots(len(panels), 1, figsize=(8.6, 10.8)); fig.subplots_adjust(hspace=0.55)
out = []
for ax, (lab, title, s) in zip(axs, panels):
    R, P, ve = pc_corr(s, vars_for(s, disc_score=lab == 'B'))
    im = ax.imshow(R.values.astype(float), cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            r, p = R.iat[i, j], P.iat[i, j]
            st = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            ax.text(j, i, f'{r:.2f}{st}', ha='center', va='center', fontsize=6.6, color='white' if abs(r) > 0.6 else '#222')
            out.append(dict(panel=title, variable=R.index[i], PC=j + 1, var_explained=ve[j], r=r, p=p))
    ax.set_xticks(range(R.shape[1])); ax.set_xticklabels([f'PC{j+1}\n{100*ve[j]:.0f}%' for j in range(R.shape[1])], fontsize=7.3)
    ax.set_yticks(range(R.shape[0])); ax.set_yticklabels(R.index, fontsize=8)
    ax.set_title(f'{lab}   {title}', loc='left', fontsize=9.5, fontweight='bold', color='#1f4e79')
    ax.tick_params(length=0)
cb = fig.colorbar(im, ax=axs, fraction=0.025, pad=0.02); cb.set_label('Pearson r (PC score vs variable)', fontsize=8)
fig.savefig(FIG + 'figS3_pc_metadata.png', dpi=200, bbox_inches='tight'); plt.close(fig)
pd.DataFrame(out).to_csv(T + 'pc_metadata_correlation.csv', index=False)
print('done')
