"""report_figures.py — figures for the report (EEM terminology). Run from the repository root."""
import os
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
T, O = "results/tables/", "report/figs/"
os.makedirs(O, exist_ok=True)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
COL = {"IVV": "#2C7BB6", "IVP": "#D7191C"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "axes.titleweight": "bold", "axes.titlesize": 10})
def save(fig, n): fig.savefig(O + n + ".png", dpi=200, bbox_inches="tight"); plt.close(fig)
REG = {"ED": "ED", "TE": "EEM"}

c = pd.read_csv("data/D15_Embryos_RawCountData.txt", sep="\t", index_col=0)
c.columns = [x[:-2] for x in c.columns]; c.index = c.index.str.replace("GeneID:", "")
meta = pd.read_csv("metadata/samples.csv").set_index("sample")
same = meta[meta.same_length].index.tolist()
def norm_log2(s):
    m = c[s]; mm = m[(m > 0).all(axis=1)]; lg = np.log(mm)
    sf = np.exp(lg.sub(lg.mean(axis=1), axis=0).median()); return np.log2(m / sf + 1)

# ---- Fig 1: PCA --------------------------------------------------------------------------
L = norm_log2(same); keep = (c[same] >= 10).sum(axis=1) >= 4; L = L[keep]
top = L.loc[L.var(axis=1).sort_values(ascending=False).index[:500]]
X = top.T.values - top.T.values.mean(axis=0); U, S, Vt = np.linalg.svd(X, full_matrices=False)
pcs = U * S; ve = 100 * S**2 / (S**2).sum()
fig, ax = plt.subplots(figsize=(6.2, 4.8))
for i, s in enumerate(same):
    t, r = meta.loc[s, "type"], meta.loc[s, "region"]; ex = s == "IVPF10_EC"
    ax.scatter(pcs[i, 0], pcs[i, 1], marker="o" if r == "ED" else "^", s=62, color="white" if ex else COL[t],
               edgecolor=COL[t], lw=1.4, zorder=3)
    ax.annotate(s, (pcs[i, 0], pcs[i, 1]), xytext=(4, 3), textcoords="offset points", fontsize=6.5, color=MUTED)
ax.set(xlabel=f"PC1 ({ve[0]:.1f}%)", ylabel=f"PC2 ({ve[1]:.1f}%)")
ax.grid(color=GRID, lw=.5); ax.set_axisbelow(True)
ax.legend([Line2D([], [], marker="o", ls="", color=COL["IVV"]), Line2D([], [], marker="o", ls="", color=COL["IVP"]),
           Line2D([], [], marker="o", ls="", color="grey"), Line2D([], [], marker="^", ls="", color="grey"),
           Line2D([], [], marker="o", ls="", mfc="white", mec=COL["IVP"])],
          ["IVV", "IVP", "ED", "EEM", "IVPF10_EC (excluded)"], fontsize=7.5, loc="best")
save(fig, "fig1_pca")

# ---- Fig 2: volcano + MA, ED and EEM ----------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
for row, reg in enumerate(["ED", "TE"]):
    d = pd.read_csv(T + f"DE_{reg}_IVP_vs_IVV_all_genes.csv").dropna(subset=["padj"])
    cols = d.direction.map({"up": COL["IVP"], "down": COL["IVV"], "ns": "#cfcfcb"})
    ax = axes[row, 0]
    ax.scatter(d.log2FoldChange, -np.log10(d.padj), s=4, c=cols, lw=0, alpha=.75, rasterized=True)
    ax.axhline(-np.log10(0.05), color=MUTED, lw=.8, ls=(0, (3, 3)))
    lab = d[d.significant].nsmallest(12, "padj")
    for _, r in lab.iterrows():
        ax.annotate(r.label, (r.log2FoldChange, -np.log10(r.padj)), fontsize=6, color=INK, xytext=(3, 2), textcoords="offset points")
    nu, nd = (d.direction == "up").sum(), (d.direction == "down").sum()
    ax.set(xlabel="log2 fold change (IVP / IVV)", ylabel="−log10 FDR",
           title=f"{REG[reg]}: {nu} higher, {nd} lower in IVP ({len(d):,} genes tested)")
    ax = axes[row, 1]; dd = d.dropna(subset=["log2FC_shrunk"])
    ax.scatter(np.log10(dd.baseMean + 1), dd.log2FC_shrunk, s=4, c=dd.direction.map({"up": COL["IVP"], "down": COL["IVV"], "ns": "#cfcfcb"}),
               lw=0, alpha=.75, rasterized=True)
    ax.axhline(0, color=MUTED, lw=.8)
    ax.set(xlabel="log10 mean normalised count", ylabel="shrunken log2 fold change", title=f"{REG[reg]}: MA plot")
    for a in axes[row]: a.grid(color=GRID, lw=.5); a.set_axisbelow(True)
fig.tight_layout(); save(fig, "fig2_volcano_MA")

# ---- Fig 3: IVP effect ED vs EEM ------------------------------------------------------------------
b = pd.read_csv(T + "IVPeffect_ED_vs_TE_per_gene.csv"); mag = pd.read_csv(T + "region_effect_magnitude.csv").iloc[0]
cat = np.select([b.sig_ED & b.sig_TE, b.sig_ED, b.sig_TE], ["both", "ED only", "EEM only"], "neither")
cc = {"neither": "#d6d5d0", "ED only": "#2a78d6", "EEM only": "#eb6834", "both": "#0b0b0b"}
fig, ax = plt.subplots(figsize=(5.4, 5.1))
for k in ["neither", "ED only", "EEM only", "both"]:
    m = cat == k
    ax.scatter(b.lfc_ED[m], b.lfc_TE[m], s=4 if k == "neither" else 7, color=cc[k], alpha=.55 if k == "neither" else .8, lw=0,
               label=f"{k} (n={m.sum():,})", rasterized=True)
lim = np.nanpercentile(np.abs(np.r_[b.lfc_ED, b.lfc_TE]), 99.9) * 1.05
ax.plot([-lim, lim], [-lim, lim], color=MUTED, lw=.8, ls=(0, (3, 3)))
ax.axhline(0, color=GRID, lw=.8, zorder=0); ax.axvline(0, color=GRID, lw=.8, zorder=0)
ax.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect="equal", xlabel="IVP vs IVV in ED (shrunken log2FC)", ylabel="IVP vs IVV in EEM (shrunken log2FC)")
ax.text(.02, .98, f"Pearson r = {mag.pearson_r_IVPeffect_ED_vs_TE:.2f}\ninteraction genes (FDR < 0.05): {int(mag.interaction_genes_FDR05)}",
        transform=ax.transAxes, va="top", fontsize=7.5)
ax.legend(title="DEG (FDR < 0.05) in", loc="lower right", fontsize=7, title_fontsize=7, markerscale=2.5)
save(fig, "fig3_ED_vs_EEM")

# ---- Fig 4: ORA dot plot ------------------------------------------------------------------------------------
cl = pd.read_csv(T + "GO_compareCluster_ED_TE_up_down.csv")
cl["Cluster"] = cl.Cluster.str.replace("TE_", "EEM_")
top = cl.sort_values("p.adjust").groupby("Cluster").head(6)
order = [x for x in ["ED_up", "ED_down", "EEM_up", "EEM_down"] if x in top.Cluster.unique()]
terms = list(dict.fromkeys(sum([top[top.Cluster == k].Description.tolist() for k in order], [])))
gr = top.GeneRatio.map(lambda s: eval(s))
fig, ax = plt.subplots(figsize=(7.4, 0.32 * len(terms) + 1.4))
sc = ax.scatter([order.index(k) for k in top.Cluster], [len(terms) - 1 - terms.index(t) for t in top.Description],
                s=gr * 2200, c=-np.log10(top["p.adjust"]), cmap="Blues", vmin=1, edgecolor=INK, lw=.6)
ax.set_xticks(range(len(order)), [k.replace("_", " ") for k in order]); ax.set_xlim(-.6, len(order) - .4)
ax.set_yticks(range(len(terms)), terms[::-1], fontsize=8); ax.grid(color=GRID, lw=.5); ax.set_axisbelow(True)
cb = fig.colorbar(sc, ax=ax, shrink=.6); cb.set_label("−log10 FDR")
for s_, lab in [(0.03, "3%"), (0.06, "6%"), (0.09, "9%")]: ax.scatter([], [], s=s_ * 2200, color="white", edgecolor=INK, label=lab)
ax.legend(title="GeneRatio", loc="upper left", bbox_to_anchor=(1.25, 1), fontsize=7, title_fontsize=7, labelspacing=1.4)
save(fig, "fig4_ORA")

# ---- Fig 5a: disc score ----------------------------------------------------------------------------------
ds = pd.read_csv(T + "qc_disc_score.csv"); ds = ds[ds.region == "ED"]
fig, ax = plt.subplots(figsize=(4.2, 4))
ds = ds.sort_values(["type", "disc_score"]).reset_index(drop=True)
for t in ["IVV", "IVP"]:
    sub = ds[ds.type == t].reset_index(drop=True)
    for k, r in sub.iterrows():
        x = (0 if t == "IVV" else 1) + (0.12 if k % 2 else -0.12); ex = r["sample"] == "IVPF10_EC"
        ax.scatter(x, r.disc_score, s=60, color="white" if ex else COL[t], edgecolor=COL[t], lw=1.4, zorder=3)
        right = k % 2 == 1
        ax.annotate(r["sample"], (x, r.disc_score), xytext=(8 if right else -8, 0), textcoords="offset points",
                    fontsize=7, va="center", ha="left" if right else "right", color=MUTED)
ax.set_xticks([0, 1], ["IVV", "IVP"]); ax.set_xlim(-.9, 1.9); ax.set_ylabel("Disc-specific transcript score")
ax.grid(axis="y", color=GRID, lw=.5); ax.set_axisbelow(True); save(fig, "fig6a_disc_score")

# ---- Fig 5b: GSEA reversal --------------------------------------------------------------------------------
gs = {k: pd.read_csv(T + f"GSEA_GO_BP_{f}.csv").set_index("Description") for k, f in
      [("ED (unadjusted)", "ED_IVP_vs_IVV"), ("ED, adjusted for disc-specific transcript score", "ED_IVP_vs_IVV_discAdjusted"), ("EEM", "TE_IVP_vs_IVV")]}
groups = [("Developmental programmes", ["anatomical structure morphogenesis", "neuron differentiation", "cell adhesion", "canonical Wnt signaling pathway"]),
          ("Metabolism (ED)", ["glycosylation", "cellular lipid metabolic process", "lipid biosynthetic process"])]
ccol = {"ED (unadjusted)": "#2a78d6", "ED, adjusted for disc-specific transcript score": "#eb6834", "EEM": "#199e70"}
mk = {"ED (unadjusted)": "o", "ED, adjusted for disc-specific transcript score": "D", "EEM": "s"}
rows = []
for g, ts in groups: rows.append(("__" + g, None)); rows += [(t, g) for t in ts]
fig, ax = plt.subplots(figsize=(7, 4.2)); y = len(rows); yt, yl = [], []
for t, g in rows:
    y -= 1
    if g is None: ax.text(-3.3, y, t[2:], fontweight="bold", fontsize=8.5, va="center"); continue
    ax.axhline(y, color=GRID, lw=.6, zorder=0); vals = []
    for k, d in gs.items():
        if t in d.index:
            v = d.loc[t, "NES"]; vals.append(v)
            ax.scatter(v, y, marker=mk[k], s=55, color=ccol[k], edgecolor="white", lw=1.2, zorder=3)
    if len(vals) > 1: ax.plot([min(vals), max(vals)], [y, y], color=MUTED, lw=1, zorder=2)
    yt.append(y); yl.append(t)
ax.axvline(0, color=MUTED, lw=1); ax.set_yticks(yt, yl, fontsize=8); ax.set_xlim(-3.3, 3); ax.set_ylim(-.7, len(rows) - .3)
ax.set_xlabel("GSEA NES (IVP vs IVV)   ← lower in IVP | higher in IVP →"); ax.tick_params(axis="y", length=0)
ax.legend([Line2D([], [], marker=mk[k], ls="", color=ccol[k], markersize=7) for k in ccol], list(ccol),
          loc="upper center", bbox_to_anchor=(.4, -.16), ncol=3, fontsize=7.5)
save(fig, "fig6b_gsea_disc_adjustment")

# ---- Fig 6: ED pathway genes vs disc score -------------------------------------------------------------------
ed_s = [s for s in same if meta.loc[s, "region"] == "ED"]; E = norm_log2(ed_s)
disc = pd.read_csv(T + "qc_disc_score.csv").set_index("sample").disc_score
a = pd.read_csv(T + "DE_ED_IVP_vs_IVV_all_genes.csv", dtype={"entrez": str}).set_index("entrez")
b2 = pd.read_csv(T + "DE_ED_IVP_vs_IVV_discAdjusted_all_genes.csv", dtype={"entrez": str}).set_index("entrez")
genes = ["WNT3", "WNT5A", "LEF1", "BMP4", "TGFB2", "DKK3", "SMAD7", "SMAD9"]
fig, axes = plt.subplots(2, 4, figsize=(11, 5.2))
for ax, g in zip(axes.flat, genes):
    eid = a.index[a.label == g][0]; x = disc[ed_s]; yv = E.loc[eid, ed_s]; main = [s for s in ed_s if s != "IVPF10_EC"]
    k, i0 = np.polyfit(x[main], yv[main], 1); xx = np.linspace(x.min() - .1, x.max() + .1, 10)
    ax.plot(xx, k * xx + i0, color=MUTED, lw=1, ls=(0, (3, 2)))
    for s in ed_s:
        t = meta.loc[s, "type"]; ex = s == "IVPF10_EC"
        ax.scatter(x[s], yv[s], s=40, color="white" if ex else COL[t], edgecolor=COL[t], lw=1.3, zorder=3)
    ax.set_title(f"{g}\nFDR {a.loc[eid, 'padj']:.0e} → adjusted {b2.loc[eid, 'padj']:.2g}", fontsize=8.5, fontweight="normal")
    ax.grid(color=GRID, lw=.5); ax.set_axisbelow(True)
for ax in axes[1]: ax.set_xlabel("Disc-specific transcript score")
for ax in axes[:, 0]: ax.set_ylabel("log2 normalised count")
fig.legend([Line2D([], [], marker="o", ls="", color=COL["IVV"]), Line2D([], [], marker="o", ls="", color=COL["IVP"]),
            Line2D([], [], marker="o", ls="", mfc="white", mec=COL["IVP"]), Line2D([], [], color=MUTED, ls=(0, (3, 2)))],
           ["IVV", "IVP", "IVPF10_EC (excluded)", "linear fit"], loc="lower center", ncol=4, bbox_to_anchor=(.5, -.03))
fig.tight_layout(rect=(0, .04, 1, 1)); save(fig, "fig7_ED_genes_vs_disc")

# ---- Fig 7: heat-shock genes both tissues --------------------------------------------------------------
an = [s for s in same if s != "IVPF10_EC"]; N = norm_log2(an)
te = pd.read_csv(T + "DE_TE_IVP_vs_IVV_all_genes.csv", dtype={"entrez": str}).set_index("entrez")
tg = ["HSPA1A", "HSPA6", "DNAJB1", "HERPUD1"]; grp = [("ED", "IVV"), ("ED", "IVP"), ("TE", "IVV"), ("TE", "IVP")]
fig, axes = plt.subplots(1, 4, figsize=(11, 3.3))
for ax, g in zip(axes, tg):
    eid = te.index[te.label == g][0]
    for j, (r, t) in enumerate(grp):
        ss = [s for s in an if meta.loc[s, "region"] == r and meta.loc[s, "type"] == t]; v = N.loc[eid, ss].values
        ax.boxplot(v, positions=[j], widths=.55, showfliers=False, medianprops=dict(color=COL[t], lw=2), boxprops=dict(color=COL[t]),
                   whiskerprops=dict(color=COL[t]), capprops=dict(color=COL[t]))
        ax.scatter(np.full(len(v), j) + np.linspace(-.12, .12, len(v)), v, s=20, color=COL[t], zorder=3, edgecolor="white", lw=.6)
    ax.set_title(f"{g}\nED FDR {a.loc[eid, 'padj']:.1e} · EEM FDR {te.loc[eid, 'padj']:.1e}", fontsize=8, fontweight="normal")
    ax.set_xticks(range(4), ["ED\nIVV", "ED\nIVP", "EEM\nIVV", "EEM\nIVP"]); ax.axvline(1.5, color=GRID)
    ax.grid(axis="y", color=GRID, lw=.5); ax.set_axisbelow(True)
axes[0].set_ylabel("log2 normalised count"); fig.tight_layout(); save(fig, "fig9_heatshock")

# ---- Fig S3: pipeline diagram --------------------------------------------------------------------------------
from matplotlib.patches import FancyBboxPatch
fig, ax = plt.subplots(figsize=(9.5, 5.4)); ax.set_xlim(0, 10); ax.set_ylim(0, 7); ax.axis("off")
def box(x, y, w, h, txt, fc="#eaf1fb"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12", fc=fc, ec="#2a78d6", lw=1.2))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=7.6)
def arr(x1, y1, x2, y2): ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2))
box(3.4, 6.1, 3.2, .7, "Raw counts, same-length conceptuses\n4 IVP vs 4 IVV; ED + EEM per embryo", "#f3f2ee")
box(3.4, 4.95, 3.2, .7, "QC: library size, PCA, correlation,\nsex (XIST), disc-specific transcript score")
arr(5, 6.1, 5, 5.65)
box(0.2, 3.5, 3.0, .85, "DESeq2 per tissue\nIVP vs IVV (~ type)\n± disc-specific transcript score")
box(3.5, 3.5, 3.0, .85, "DESeq2 joint model (paired)\n~ type + type:embryo + type:tissue\n→ tissue × treatment interaction")
box(6.8, 3.5, 3.0, .85, "Classifier sanity check\nleave-one-embryo-out,\nexact permutation test")
for x in (1.7, 5.0, 8.3): arr(5, 4.95, x, 4.35)
box(0.2, 1.9, 3.0, .9, "GO enrichment\nORA (tested-gene background, BH)\nGSEA (all genes ranked)")
box(3.5, 1.9, 3.0, .9, "Gene-level plots\nper-gene IVV vs IVP;\nexpression vs disc-specific score")
arr(1.7, 3.5, 1.7, 2.8); arr(5.0, 3.5, 5.0, 2.8)
box(2.0, 0.35, 6.0, .85, "Interpretation with limitations: sex (all IVP male, all IVV female),\ntwo dams per group, n = 4 per group", "#fbebeb")
arr(1.7, 1.9, 4, 1.2); arr(5, 1.9, 5, 1.2); arr(8.3, 3.5, 6.8, 1.2)
save(fig, "figS1_workflow")
print("ok")
