"""
classifier_figures.py — ROC curves, accuracy and per-embryo predictions for 07_classifier_by_region.R

Run from the repository root after 07 (re-run this whenever 07 is re-run):
    python3 python/classifier_figures.py
Reads  results/tables/classifier_by_region_{predictions,summary,null}.csv
Writes results/figures/classifier_{roc,accuracy,per_embryo,panel}.{png,pdf} and report/figs/figS2_classifier.png
"""
import os
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

T, F = "results/tables/", "results/figures/"
pred = pd.read_csv(T + "classifier_by_region_predictions.csv")
summ = pd.read_csv(T + "classifier_by_region_summary.csv")
null = pd.read_csv(T + "classifier_by_region_null.csv")

ORDER  = ["within_ED", "within_TE", "train_TE_test_ED", "train_ED_test_TE"]
TITLES = {"within_ED": "ED only", "within_TE": "EEM only",
          "train_TE_test_ED": "Train EEM → test ED", "train_ED_test_TE": "Train ED → test EEM"}
TICKS  = {"within_ED": "ED only", "within_TE": "EEM only",
          "train_TE_test_ED": "Train EEM\n→ test ED", "train_ED_test_TE": "Train ED\n→ test EEM"}
NF_MAIN = 100
NF_COL  = {50: "#86b6ef", 100: "#2a78d6", 500: "#104281"}     # one hue, light -> dark
NF_LS   = {50: ":", 100: "-", 500: "--"}
MOD_COL = {"SVM": "#6250d6", "RF": "#199e70"}
TYPE_COL = {"IVV": "#2C7BB6", "IVP": "#D7191C"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titlesize": 10, "axes.titleweight": "bold",
                     "legend.frameon": False})

# P(IVP) for every held-out embryo (RF vote fraction)
pred["p_ivp"] = np.where(pred.truth == "IVP", pred.rf_prob_true, 1 - pred.rf_prob_true)

def roc(y, s):
    """ROC points and AUC (Mann–Whitney, ties = 0.5) for labels y (1 = IVP) and scores s."""
    thr = np.r_[np.inf, np.sort(np.unique(s))[::-1]]
    tpr = [((s >= t) & (y == 1)).sum() / (y == 1).sum() for t in thr]
    fpr = [((s >= t) & (y == 0)).sum() / (y == 0).sum() for t in thr]
    pos, neg = s[y == 1], s[y == 0]
    auc = np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg])
    return np.array(fpr), np.array(tpr), auc

def save(fig, name):
    fig.savefig(F + name + ".png", dpi=220, bbox_inches="tight")
    fig.savefig(F + name + ".pdf", bbox_inches="tight")
    plt.close(fig)

# ---- 1. ROC curves ----------------------------------------------------------------
def draw_roc(axes):
    aucs = []
    for ax, a in zip(axes, ORDER):
        ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls=(0, (3, 3)))
        for nf in [50, 100, 500]:
            d = pred[(pred.analysis == a) & (pred.n_features == nf)]
            fpr, tpr, auc = roc((d.truth == "IVP").astype(int).values, d.p_ivp.values)
            aucs.append((a, nf, auc))
            ax.step(fpr, tpr, where="post", color=NF_COL[nf], ls=NF_LS[nf],
                    lw=2.2 if nf == NF_MAIN else 1.6, label=f"{nf} genes  AUC {auc:.2f}")
        ax.set(xlim=(-0.03, 1.03), ylim=(-0.03, 1.03), aspect="equal",
               xlabel="False positive rate", title=TITLES[a])
        ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.legend(loc="lower right", fontsize=7.5, handlelength=2.2)
    axes[0].set_ylabel("True positive rate (IVP)")
    return pd.DataFrame(aucs, columns=["analysis", "n_features", "auc"])

fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), sharey=True)
auc_tab = draw_roc(axes)
fig.suptitle("ROC: random forest, leave-one-embryo-out (4 IVP vs 4 IVV embryos per curve)",
             x=0.01, ha="left", fontsize=10.5, color=INK)
save(fig, "classifier_roc")
auc_tab.to_csv(T + "classifier_by_region_auc.csv", index=False)

# ---- 2. accuracy ---------------------------------------------------------------------
def draw_accuracy(ax):
    x = np.arange(len(ORDER)); w = 0.42
    for j, (mod, col) in enumerate([("SVM", "accuracy_svm"), ("RF", "accuracy_rf")]):
        pcol = "p_svm" if mod == "SVM" else "p_rf"
        xs = x + (j - 0.5) * w
        main = summ[summ.n_features == NF_MAIN].set_index("analysis").loc[ORDER]
        ax.bar(xs, main[col], width=w - 0.04, color=MOD_COL[mod], label=f"{mod} ({NF_MAIN} genes)",
               edgecolor="white", linewidth=2)
        for nf, mk, dx in [(50, "v", -0.08), (500, "^", 0.08)]:
            other = summ[summ.n_features == nf].set_index("analysis").loc[ORDER]
            ax.scatter(xs + dx, other[col] + 0.035, marker=mk, s=26, color=INK, zorder=3,
                       edgecolor="white", linewidth=0.6)
        for xi, acc, p in zip(xs, main[col], main[pcol]):
            ax.text(xi, acc - 0.03, f"{acc:.2f}\np=" + f"{p:.3f}".lstrip("0"), ha="center", va="top",
                    fontsize=6.2, color="white", fontweight="bold")
    # 95th percentile of the permutation null (accuracy, RF and SVM pooled) per analysis
    for xi, a in zip(x, ORDER):
        nl = null[(null.analysis == a) & (null.n_features == NF_MAIN)]
        q = np.quantile(np.r_[nl.svm, nl.rf], 0.95)
        ax.hlines(q, xi - 0.46, xi + 0.46, color=MUTED, lw=1.4, ls=(0, (4, 2)))
    ax.axhline(0.5, color=MUTED, lw=0.8, ls=":")
    ax.set_xticks(x, [TICKS[a] for a in ORDER])
    ax.set(ylim=(0, 1.12), ylabel="LOEO accuracy (8 embryos)", yticks=np.arange(0, 1.01, 0.25),
           title="Classification accuracy")
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=MOD_COL[m]) for m in MOD_COL] + [
        Line2D([], [], marker="v", ls="", color=INK, label="50 genes"),
        Line2D([], [], marker="^", ls="", color=INK, label="500 genes"),
        Line2D([], [], color=MUTED, ls=(0, (4, 2)), label="null 95th pct"),
        Line2D([], [], color=MUTED, ls=":", label="chance")]
    ax.legend(handles, [f"SVM ({NF_MAIN} genes)", f"RF ({NF_MAIN} genes)", "50 genes", "500 genes",
                        "permutation null, 95th pct", "chance"],
              loc="upper center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, -0.17))

fig, ax = plt.subplots(figsize=(7.5, 4.6))
draw_accuracy(ax)
ax.text(0, -0.42, "p = exact permutation p over 70 embryo-level labellings; minimum attainable = 2/70 = 0.029",
        transform=ax.transAxes, fontsize=7, color=MUTED)
save(fig, "classifier_accuracy")

# ---- 3. per-embryo predicted probability -------------------------------------------
def draw_per_embryo(ax):
    d = pred[pred.n_features == NF_MAIN]
    for i, a in enumerate(ORDER):
        di = d[d.analysis == a].sort_values("p_ivp")
        jitter = np.linspace(-0.3, 0.3, len(di))
        for k, ((_, r), jx) in enumerate(zip(di.iterrows(), jitter)):
            ax.scatter(i + jx, r.p_ivp, s=60, color=TYPE_COL[r.truth], edgecolor="white",
                       linewidth=1.5, zorder=3)
            right = k % 2 == 0
            ax.annotate(r.embryo, (i + jx, r.p_ivp), xytext=(5 if right else -5, 0),
                        textcoords="offset points", fontsize=6, color=MUTED, va="center",
                        ha="left" if right else "right")
    ax.axhline(0.5, color=MUTED, lw=1, ls=(0, (3, 3)))
    ax.text(-0.55, 0.515, "threshold 0.5", fontsize=7, color=MUTED, va="bottom", ha="left")
    ax.set_xticks(range(len(ORDER)), [TICKS[a] for a in ORDER])
    ax.set(ylim=(0, 1), xlim=(-0.6, len(ORDER) - 0.4), ylabel="Predicted P(IVP), random forest",
           title=f"Held-out embryo predictions ({NF_MAIN} genes)")
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.legend([Line2D([], [], marker="o", ls="", color=TYPE_COL[t], markersize=8) for t in TYPE_COL],
              ["true IVV", "true IVP"], loc="upper left", fontsize=8, ncol=2)

fig, ax = plt.subplots(figsize=(8, 4.4))
draw_per_embryo(ax)
save(fig, "classifier_per_embryo")

# ---- 4. combined panel -------------------------------------------------
fig = plt.figure(figsize=(13, 8.6))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.35], hspace=0.42, wspace=0.28)
draw_roc([fig.add_subplot(gs[0, i]) for i in range(4)])
draw_accuracy(fig.add_subplot(gs[1, :2]))
draw_per_embryo(fig.add_subplot(gs[1, 2:]))
fig.text(0.07, 0.965, "All curves overlap at AUC = 1.00: every held-out embryo is ranked correctly in every analysis",
         fontsize=8, color=MUTED)
for lbl, (xx, yy) in zip("ABC", [(0.06, 0.92), (0.06, 0.52), (0.51, 0.52)]):
    fig.text(xx, yy, lbl, fontsize=14, fontweight="bold", color=INK)
save(fig, "classifier_panel")
os.makedirs("report/figs", exist_ok=True)
fig.savefig("report/figs/figS2_classifier.png", dpi=220, bbox_inches="tight")  # Supplementary Figure S2 of the report

print(auc_tab.pivot(index="analysis", columns="n_features", values="auc").loc[ORDER].round(3))
print("figures written to", F)
