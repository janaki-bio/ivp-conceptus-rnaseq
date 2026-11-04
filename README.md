# IVP vs IVV day-15 bovine conceptuses: embryonic disc and extra-embryonic membrane RNA-seq analysis

Reproducible analysis from an MSc research project (Utrecht University). It compares the transcriptome of the
**embryonic disc (ED)** and **extra-embryonic membranes (EEM)** of day-15 bovine conceptuses produced **in vitro (IVP)**
or **in vivo (IVV)**, using only conceptuses of the same length (4 IVP vs 4 IVV; ED and EEM from each embryo).

This repository is not the canonical publication of the dataset. The samples come from the study published as Behrens et al. (2025), *Biol Reprod* 115(2):348–360, https://doi.org/10.1093/biolre/ioaf095.

The full report is in [`report/embryo_transcriptomics_report.pdf`](report/embryo_transcriptomics_report.pdf).

## Main findings (exploratory; group differences, not effects of IVP; see limitations)

- DEGs between the IVP and IVV groups (DESeq2, BH FDR < 0.05): **ED 1,218** (467 higher / 751 lower in IVP); **EEM 721** (423 / 298).
  For **417 genes** the IVP–IVV contrast differs between tissues (paired tissue × group interaction).
- **ED glycosylation and lipid metabolism higher in the IVP group:** the most consistent pattern; it keeps its direction
  when any single conceptus is omitted, after removing X-linked genes and in the full cohort. The disc-score-adjusted
  version depends on one IVP conceptus (IVP14).
- **ED developmental gene sets (morphogenesis, Wnt, BMP/TGF-β, Notch) lower in the IVP group:** sensitive to adjustment
  for an internally derived disc-enrichment score. This is compatible with less disc tissue *or* a different disc-cell
  state; bulk RNA-seq cannot tell them apart.
- **EEM adhesion/extracellular-matrix gene sets higher in the IVP group:** consistent in direction, weaker outside the
  primary comparison. Lower translation/respiration gene sets reverse in other embryos (source unresolved).
- X-linked genes are over-represented among DEGs (8–10% vs 3%), but removing them does not change the gene-set results.
- Inducible heat-shock genes and HEY1/ID1 are lower in the IVP group in both tissues in every analysis (cause unresolved).
- A classifier separates the groups perfectly (supplementary, exploratory); this says nothing about *why* they separate.

## Important limitations

- **Sex confounding:** in the primary comparison every IVP conceptus is inferred male and every IVV conceptus female (XIST).
- **Maternal clustering:** IVP conceptuses were carried by two recipients, IVV conceptuses came from two superovulated
  donors; embryo-to-cow assignments are not available and the models do not account for shared cows.
- **Unresolved sample-sheet label:** the four IVP conceptuses are labelled "Vitrified" (all IVP embryos were transferred
  fresh); its meaning is unresolved and associated variation cannot be separated from the group contrast (see `data/README.md`).
- **Selection:** conceptuses were selected by length category and presence of a disc; measured lengths are not available.
- **Small sample and bulk tissue:** n = 4 per group; ED sections mix disc and membrane tissue, and the disc-specific
  transcript score is an expression-derived proxy, not a measured cell proportion.

## Terminology note

In the code and file names, the non-disc tissue is labelled **`TE`**, for example `region == "TE"` and
`DE_TE_IVP_vs_IVV_*.csv`. In the report it is called **EEM** (extra-embryonic membranes), as in Behrens et al. (2025).
`disc_score` in the code is the "disc-specific transcript score" in the report.

## Repository layout

```text
R/                 analysis pipeline (run in order by run_all.R)
  00_config.R      paths, parameters (FDR, filters, exclusions, seed), helpers
  01_qc.R          library sizes, PCA, correlations, disc-specific transcript score
  02_de_by_region.R  IVP vs IVV per tissue (+ sensitivity analyses)
  03_joint_model.R paired model of all 16 libraries; tissue × treatment interaction
  04_enrichment.R  GO ORA (tested-gene background) and GSEA (all genes ranked)
  05_figures.R     volcano/MA, heatmaps, GO and gene-level plots
  06_classifier_loeo.R      leave-one-embryo-out classifier + exact permutation test
  07_classifier_by_region.R within-tissue and cross-tissue classification
  08_full_cohort_checks.R   male-only and full-cohort sensitivity analyses
  09_xlinked_check.R        contribution of X-linked genes
  10_sensitivity.R          leave-one-embryo-out influence and disc-marker sensitivity
python/            figure scripts for the report (read results/tables)
metadata/samples.csv   one row per library (embryo, type, tissue, section, same-length flag)
data/              put the count matrix here (not committed; see data/README.md)
results/tables/    final result tables (CSV)
results/figures/   pipeline figures (PNG + PDF)
report/            report (PDF + HTML source + figures)
```

## Models

| Analysis | DESeq2 design | Notes |
|---|---|---|
| IVP vs IVV per tissue (primary) | `~ type` | 4 vs 4 embryos; IVV = reference |
| ED, composition-adjusted (sensitivity) | `~ disc_score + type` | disc_score scaled; conservative |
| Tissue × treatment interaction | `~ type + type:embryo_n + type:region` | embryo_n = 1–4 within type; contrast `list("typeIVP.regionTE", "typeIVV.regionTE")`; 6 residual df |

Genes are kept if they have ≥ 10 reads in ≥ 4 libraries. Wald tests are run on raw counts, with BH FDR < 0.05 and
apeglm shrinkage for ranking and plots. Enrichment uses clusterProfiler ORA (background = tested genes) and gseGO
(ranked by the Wald statistic). The classifier uses linear SVM and random forest, leave-one-embryo-out validation,
feature selection inside folds, and an exact permutation test over 70 labellings. Seeds are fixed (`params$seed`).

## How to reproduce

1. Obtain the gene-level count matrix from the corresponding author of Behrens et al. (2025) and place it at `data/D15_Embryos_RawCountData.txt` (see `data/README.md`).
2. Open `ivp_conceptus_rnaseq.Rproj` (R ≥ 4.3; developed with R 4.4.1, DESeq2 1.44.0).
3. Install packages once with `source("install_packages.R")`.
4. Run everything with `source("run_all.R")`. This takes about 20–40 min, mostly classifier permutations and GSEA, and writes `results/sessionInfo.txt`.
5. Optional, for the report figures: run `python3 python/report_figures.py`, `python3 python/report_pathway_figures.py`, `python3 python/report_signalling_pca.py`, `python3 python/report_leading_edge.py`,
   `python3 python/report_sensitivity_figures.py` and `python3 python/classifier_figures.py`
   from the repository root (needs pandas, numpy, matplotlib).
6. Optional, for an HTML summary: `quarto::quarto_render("report.qmd")`.

The package versions used for the results are recorded in `renv.lock`. To install the same versions, run
`install.packages("renv"); renv::restore()` from the repository root (instead of step 3).

## Citation

If you use this code, please cite the repository (see `CITATION.cff`) and Behrens et al. (2025).

## Licence

Code: MIT (see `LICENSE`). Report text and figures: CC BY 4.0. The underlying sequencing data are governed by their
GEO deposition (GSE278243).
