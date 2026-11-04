# =============================================================================
# 07_classifier_by_region.R — predict IVP vs IVV from ONE region, and across regions
#
# Four analyses, all leave-one-embryo-out (LOEO) over the 8 same-length embryos:
#   within_ED         train on ED of 7 embryos  -> predict ED of the held-out embryo
#   within_TE         train on TE of 7 embryos  -> predict TE of the held-out embryo
#   train_TE_test_ED  train on TE of 7 embryos  -> predict ED of the held-out embryo
#   train_ED_test_TE  train on ED of 7 embryos  -> predict TE of the held-out embryo
# The held-out embryo is excluded from training in BOTH regions, so embryo identity
# (genetic background, sex, ...) cannot leak from the training set into the test sample.
#
# Interpretation
#   * Within-region accuracy will very likely be at the ceiling (1.00) for both regions,
#     because IVP and IVV differ at >700 genes in each. A ceiling cannot rank the regions,
#     so the RF mean probability assigned to the true class is also reported as a graded
#     measure.
#   * Cross-region transfer asks whether the IVP signature is SHARED between tissues:
#       high accuracy      -> common IVP signature in ED and TE
#       ~0.5 (chance)      -> region-specific IVP effects
#       well below 0.5     -> signature present but INVERTED between regions
#   * ED predictions may partly reflect lower disc content in IVP sections (see 01_qc.R),
#     not only an IVP expression signature.
#   * IVP is confounded with embryo sex and dam; "IVP" here means IVP + sex + dam.
#
# Statistics: exact permutation test. The whole procedure, including feature selection,
# is repeated for all 70 ways of labelling 4 of the 8 embryos as IVP. Swapping every
# label gives the same accuracy, so the smallest attainable p-value is 2/70 = 0.029.
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(e1071); library(randomForest) })

n_features_grid <- c(50, 100, 500)   # robustness: results should not hinge on one choice

inp  <- readRDS(file.path(paths$objects, "inputs.rds"))
meta <- analysis_samples(inp$samples)          # 16 libraries = 8 embryos x {ED, TE}
cts  <- prefilter(inp$counts[, meta$sample])
# blind VST does not use labels, so computing it once on all samples is leakage-free
V    <- assay(vst(DESeqDataSetFromMatrix(cts, meta, ~ 1), blind = TRUE))

emb_tab  <- meta %>% distinct(embryo, type) %>% arrange(type, embryo)
embryos  <- emb_tab$embryo
true_lab <- setNames(as.character(emb_tab$type), embryos)

# genes x embryos matrix for one region (columns aligned by embryo across regions)
region_mat <- function(reg) {
  m <- meta[meta$region == reg, ]
  M <- V[, m$sample, drop = FALSE]
  colnames(M) <- m$embryo
  M[, embryos]
}
R_raw <- list(ED = region_mat("ED"), TE = region_mat("TE"))

# z-score every gene WITHIN each region. This is label-free and removes the large ED-vs-TE
# baseline difference, which would otherwise make cross-region prediction meaningless.
# (It assumes each region's samples are roughly balanced between IVP and IVV, as they are here: 4 vs 4.)
keep <- matrixStats::rowSds(R_raw$ED) > 0 & matrixStats::rowSds(R_raw$TE) > 0
Z <- lapply(R_raw, function(M) {
  M <- M[keep, , drop = FALSE]
  (M - rowMeans(M)) / matrixStats::rowSds(M)
})

row_welch_t <- function(M, g) {
  m1 <- rowMeans(M[, g, drop = FALSE]); m0 <- rowMeans(M[, !g, drop = FALSE])
  v1 <- matrixStats::rowVars(M[, g, drop = FALSE]); v0 <- matrixStats::rowVars(M[, !g, drop = FALSE])
  (m1 - m0) / sqrt(v1 / sum(g) + v0 / sum(!g) + 1e-8)
}

predict_loeo <- function(src, tgt, labels, nf) {
  y <- factor(labels[embryos], levels = c("IVV", "IVP"))
  out <- tibble(embryo = embryos, truth = y,
                svm = factor(NA, levels = levels(y)), rf = factor(NA, levels = levels(y)),
                rf_prob_true = NA_real_)
  for (i in seq_along(embryos)) {
    tr    <- setdiff(seq_along(embryos), i)                      # held-out embryo excluded
    tt    <- row_welch_t(Z[[src]][, tr], y[tr] == "IVP")         # feature selection: training only
    feats <- names(sort(abs(tt), decreasing = TRUE))[seq_len(nf)]
    xtr   <- t(Z[[src]][feats, tr])
    xte   <- t(Z[[tgt]][feats, i, drop = FALSE])
    fit_svm <- svm(xtr, y[tr], kernel = "linear", cost = 1, scale = TRUE)
    set.seed(params$seed + i)   # same forest for identical inputs -> observed run == its permutation twin
    fit_rf  <- randomForest(xtr, y[tr], ntree = 500)
    out$svm[i] <- predict(fit_svm, xte)
    out$rf[i]  <- predict(fit_rf, xte)
    out$rf_prob_true[i] <- predict(fit_rf, xte, type = "prob")[1, as.character(y[i])]
  }
  out
}

score <- function(p) c(svm     = mean(p$svm == p$truth),
                       rf      = mean(p$rf  == p$truth),
                       rf_prob = mean(p$rf_prob_true))

analyses <- tribble(
  ~analysis,          ~src, ~tgt,
  "within_ED",        "ED", "ED",
  "within_TE",        "TE", "TE",
  "train_TE_test_ED", "TE", "ED",
  "train_ED_test_TE", "ED", "TE"
)
combos      <- combn(embryos, 4, simplify = FALSE)
perm_labels <- lapply(combos, function(s) setNames(ifelse(embryos %in% s, "IVP", "IVV"), embryos))

summary_rows <- list(); pred_rows <- list(); null_rows <- list()
for (k in seq_len(nrow(analyses))) {
  a <- analyses[k, ]
  for (nf in n_features_grid) {
    message(sprintf("%-18s n_features = %d", a$analysis, nf))
    obs_pred <- predict_loeo(a$src, a$tgt, true_lab, nf)
    obs      <- score(obs_pred)
    null     <- t(vapply(perm_labels, function(l) score(predict_loeo(a$src, a$tgt, l, nf)),
                         numeric(3)))

    summary_rows[[length(summary_rows) + 1]] <- tibble(
      analysis = a$analysis, train_region = a$src, test_region = a$tgt, n_features = nf,
      accuracy_svm = unname(obs["svm"]), p_svm = mean(null[, "svm"] >= obs["svm"]),
      accuracy_rf  = unname(obs["rf"]),  p_rf  = mean(null[, "rf"]  >= obs["rf"]),
      rf_mean_prob_true_class = unname(obs["rf_prob"]),
      p_rf_prob = mean(null[, "rf_prob"] >= obs["rf_prob"]),
      null_rf_prob_95pct = unname(quantile(null[, "rf_prob"], 0.95)))
    pred_rows[[length(pred_rows) + 1]] <- mutate(obs_pred, analysis = a$analysis, n_features = nf)
    null_rows[[length(null_rows) + 1]] <- as_tibble(null) %>%
      mutate(analysis = a$analysis, n_features = nf, permutation = seq_len(n()))
  }
}

summ  <- bind_rows(summary_rows)
preds <- bind_rows(pred_rows)
nulls <- bind_rows(null_rows)
write_csv(summ,  file.path(paths$tables, "classifier_by_region_summary.csv"))
write_csv(preds, file.path(paths$tables, "classifier_by_region_predictions.csv"))
write_csv(nulls, file.path(paths$tables, "classifier_by_region_null.csv"))
print(summ, width = Inf)

# ---- figure: graded measure (RF probability of the true class) vs permutation null ----
lvl <- analyses$analysis
p <- ggplot(mutate(nulls, analysis = factor(analysis, lvl)),
            aes(factor(n_features), rf_prob)) +
  geom_violin(fill = "grey85", colour = NA) +
  geom_point(data = mutate(summ, analysis = factor(analysis, lvl)),
             aes(y = rf_mean_prob_true_class), colour = "#D7191C", size = 3) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  facet_wrap(~ analysis, nrow = 1) +
  labs(x = "Number of selected genes", y = "Mean RF probability of the true class",
       title = "Predicting IVP vs IVV from one region, and across regions",
       subtitle = "Red: observed (leave-one-embryo-out). Grey: all 70 embryo-level label permutations.")
save_plot(p, "classifier_by_region", 11, 4.2)

message("07_classifier_by_region.R done")
