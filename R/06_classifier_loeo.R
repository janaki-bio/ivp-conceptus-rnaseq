# =============================================================================
# 06_classifier_loeo.R — honest IVP-vs-IVV classifier (optional validation)
#
# Design choices that avoid information leakage:
#   * Selecting features from DEGs computed on all samples would leak the labels into the features
#     and inflate accuracy (any random labelling would also look predictable).
#   * ED and TE from one embryo are not independent; putting one in training and the other in
#     testing would leak embryo identity.
#
# Here:
#   * leave-one-EMBRYO-out CV (both libraries of the held-out embryo are removed together);
#   * feature selection (top N genes by Welch t on the training fold only) happens INSIDE each fold;
#   * an exact permutation test re-runs the whole procedure for all 70 ways of assigning
#     4 of the 8 embryos to "IVP" — this is the null distribution of accuracy.
# With 8 embryos the result is necessarily coarse; treat it as a sanity check, not a biomarker.
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(e1071); library(randomForest) })

n_features <- 100
inp  <- readRDS(file.path(paths$objects, "inputs.rds"))
meta <- analysis_samples(inp$samples)
cts  <- prefilter(inp$counts[, meta$sample])
# blind VST: unsupervised (does not use labels), so it is safe to compute once on all samples
X    <- assay(vst(DESeqDataSetFromMatrix(cts, meta, ~ 1), blind = TRUE))
X    <- X[rowVars(X) > 0, ]

row_welch_t <- function(M, g) {
  m1 <- rowMeans(M[, g, drop = FALSE]); m0 <- rowMeans(M[, !g, drop = FALSE])
  v1 <- rowVars(M[, g, drop = FALSE]);  v0 <- rowVars(M[, !g, drop = FALSE])
  (m1 - m0) / sqrt(v1 / sum(g) + v0 / sum(!g) + 1e-8)
}

loeo_predict <- function(labels_by_embryo) {
  y <- factor(labels_by_embryo[meta$embryo], levels = c("IVV", "IVP"))
  pred <- tibble(sample = meta$sample, embryo = meta$embryo, truth = y,
                 svm = factor(NA, levels = levels(y)), rf = factor(NA, levels = levels(y)))
  for (emb in unique(meta$embryo)) {
    test  <- meta$embryo == emb
    train <- !test
    tt    <- row_welch_t(X[, train], y[train] == "IVP")
    feats <- names(sort(abs(tt), decreasing = TRUE))[seq_len(n_features)]
    xtr <- t(X[feats, train]); xte <- t(X[feats, test, drop = FALSE])
    fit_svm <- svm(xtr, y[train], kernel = "linear", cost = 1, scale = TRUE)
    set.seed(params$seed + match(emb, unique(meta$embryo)))   # reproducible forest per fold
    fit_rf  <- randomForest(xtr, y[train], ntree = 500)
    pred$svm[test] <- predict(fit_svm, xte)
    pred$rf[test]  <- predict(fit_rf,  xte)
  }
  pred
}

embryos <- meta %>% distinct(embryo, type)
true_labels <- setNames(as.character(embryos$type), embryos$embryo)
pred <- loeo_predict(true_labels)
acc  <- c(svm = mean(pred$svm == pred$truth), rf = mean(pred$rf == pred$truth))
write_csv(pred, file.path(paths$tables, "classifier_loeo_predictions.csv"))

# exact permutation null: all 70 label assignments at the embryo level
combos <- combn(embryos$embryo, 4, simplify = FALSE)
null_acc <- map_dfr(combos, function(ivp_set) {
  lab <- setNames(ifelse(embryos$embryo %in% ivp_set, "IVP", "IVV"), embryos$embryo)
  p   <- loeo_predict(lab)
  tibble(svm = mean(p$svm == p$truth), rf = mean(p$rf == p$truth))
})
perm_p <- c(svm = mean(null_acc$svm >= acc["svm"]), rf = mean(null_acc$rf >= acc["rf"]))

result <- tibble(model = names(acc), loeo_accuracy = acc, permutation_p = perm_p,
                 n_permutations = length(combos), n_features = n_features)
write_csv(result, file.path(paths$tables, "classifier_loeo_summary.csv"))
write_csv(null_acc, file.path(paths$tables, "classifier_loeo_null_accuracies.csv"))
print(result)

p <- null_acc %>% pivot_longer(everything(), names_to = "model", values_to = "accuracy") %>%
  ggplot(aes(accuracy)) + geom_histogram(binwidth = 1 / 16, fill = "grey70", colour = "white") +
  geom_vline(data = result, aes(xintercept = loeo_accuracy), colour = "#D7191C", linewidth = 1) +
  facet_wrap(~ model) +
  labs(x = "Leave-one-embryo-out accuracy", y = "Label permutations",
       title = "Observed accuracy (red) vs all 70 embryo-level label permutations")
save_plot(p, "classifier_loeo_permutation", 8, 4)

message("06_classifier_loeo.R done")
