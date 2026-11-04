# =============================================================================
# 08_full_cohort_checks.R — two pre-stated robustness checks using all embryos
#
# The same-length set (primary analysis) has all-male IVP vs all-female IVV embryos. These checks use the
# full 32-library cohort to test whether the main findings survive. Decide in advance what counts as
# "holds": same direction of effect, and for key genes/terms, FDR < 0.05 where power allows.
#
# Check 1 — SEX-MATCHED CONTRAST: IVP males vs IVV males (VE3, VE5), per tissue.
#           Only males exist in both groups. With 2 IVV males, power is low: judge by DIRECTION
#           agreement with the primary analysis, not by DEG counts.
# Check 3 — FULL-COHORT REPLICATION: all embryos per tissue, design ~ sex + length_group + type
#           (length_group = same-length vs other, a crude proxy; true lengths were not available),
#           ED also with the disc-content score added.
#
# Caveats written into the outputs:
#   * In the ED, every IVP embryo is male, so the sex coefficient is estimated from the 2 IVV males only.
#   * Adjusting for length also removes part of the real IVP effect (IVP embryos are shorter).
#   * IVPF13 is excluded: its ED and EEM pieces disagree on inferred sex (possible mislabel).
#
# Requires: 01, 02 (and qc_sex_inference_XIST.csv in results/tables).
# Outputs:  results/tables/fullcohort_*.csv, results/figures/fullcohort_*.png
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(apeglm); library(clusterProfiler) })

exclude_full <- c("IVPF13_EA", "IVPF13_T1A")
key_genes <- c("WNT10B", "SMAD7", "SMAD9", "WNT3", "WNT5A", "LEF1", "BMP4",
               "HSPA1A", "HSPA6", "DNAJB1", "HERPUD1")
key_terms <- list(
  ED = c("anatomical structure morphogenesis", "neuron differentiation", "cell adhesion",
         "glycosylation", "cellular lipid metabolic process"),
  TE = c("translation", "ribosome biogenesis", "protein folding", "cellular respiration",
         "cell migration", "cell adhesion")
)

inp  <- readRDS(file.path(paths$objects, "inputs.rds"))
disc <- readRDS(file.path(paths$objects, "disc_score.rds"))
sex  <- read_csv(file.path(paths$tables, "qc_sex_inference_XIST.csv"), show_col_types = FALSE) %>%
  select(sample, sex = sex_inferred)

meta_all <- inp$samples %>%
  filter(!sample %in% exclude_full) %>%
  left_join(sex, by = "sample") %>%
  mutate(type         = factor(type, levels = c("IVV", "IVP")),
         sex          = factor(sex, levels = c("female", "male")),
         length_group = factor(ifelse(same_length, "matched", "other"), levels = c("matched", "other"))) %>%
  as.data.frame()
rownames(meta_all) <- meta_all$sample
stopifnot(!anyNA(meta_all$sex))

# disc-content score for ALL ED libraries (same markers and same blind VST as 01_qc.R)
vsd_all <- assay(vst(DESeqDataSetFromMatrix(prefilter(inp$counts), inp$samples, ~ 1), blind = TRUE))
meta_all$disc_score <- colMeans(vsd_all[disc$markers, meta_all$sample])

design_tab <- meta_all %>% count(region, type, sex, length_group)
write_csv(design_tab, file.path(paths$tables, "fullcohort_design.csv"))
print(design_tab)

# ---- generic DESeq2 runner ------------------------------------------------------------
run_model <- function(meta, design, tag) {
  meta <- droplevels(meta)
  cts  <- prefilter(inp$counts[, meta$sample])
  dds  <- DESeqDataSetFromMatrix(cts, colData = meta, design = design)
  dds  <- DESeq(dds, quiet = TRUE)
  res  <- results(dds, name = "type_IVP_vs_IVV", alpha = params$alpha)
  shr  <- lfcShrink(dds, coef = "type_IVP_vs_IVV", type = "apeglm", quiet = TRUE)
  tab  <- tidy_results(res, shr)
  write_de_tables(tab, paste0("fullcohort_", tag))
  message(sprintf("%-32s n = %2d  DEGs %5d (%d up, %d down)", tag, ncol(dds), sum(tab$significant),
                  sum(tab$direction == "up"), sum(tab$direction == "down")))
  tab
}

full <- list()
for (reg in c("ED", "TE")) {
  m <- meta_all[meta_all$region == reg, ]

  # Check 1: sex-matched (males only)
  mm <- m[m$sex == "male", ]
  stopifnot(nlevels(droplevels(mm$type)) == 2, all(table(droplevels(mm$type)) >= 2))
  full[[paste0(reg, "_sexmatched")]] <- run_model(mm, ~ type, paste0(reg, "_sexmatched_males"))

  # Check 3: full cohort, sex + length group
  full[[paste0(reg, "_full")]] <- run_model(m, ~ sex + length_group + type, paste0(reg, "_full_sex_length"))
  if (reg == "ED") {
    m2 <- m; m2$disc_score <- as.numeric(scale(m2$disc_score))
    full[["ED_full_disc"]] <- run_model(m2, ~ sex + length_group + disc_score + type, "ED_full_sex_length_disc")
  }
}

# ---- compare with the primary (same-length) analyses ------------------------------------
primary <- list(ED = readRDS(file.path(paths$objects, "de_ED_IVP_vs_IVV.rds"))$table,
                TE = readRDS(file.path(paths$objects, "de_TE_IVP_vs_IVV.rds"))$table,
                ED_disc = readRDS(file.path(paths$objects, "de_ED_IVP_vs_IVV_discAdjusted.rds"))$table)

compare_to_primary <- function(full_tab, prim_tab, label) {
  j <- inner_join(select(prim_tab, entrez, lfc_p = log2FoldChange, sig_p = significant),
                  select(full_tab, entrez, lfc_f = log2FoldChange, sig_f = significant),
                  by = "entrez")
  p <- filter(j, sig_p, !is.na(lfc_f))
  tibble(analysis                  = label,
         primary_DEGs_compared     = nrow(p),
         same_direction_pct        = round(100 * mean(sign(p$lfc_p) == sign(p$lfc_f)), 1),
         replicated_FDR05_same_dir = sum(p$sig_f & sign(p$lfc_p) == sign(p$lfc_f), na.rm = TRUE),
         spearman_all_genes        = round(cor(j$lfc_p, j$lfc_f, method = "spearman", use = "complete.obs"), 3),
         DEGs_in_check             = sum(full_tab$significant))
}
comparison <- bind_rows(
  compare_to_primary(full$ED_sexmatched, primary$ED,      "ED: sex-matched (IVP males vs IVV males)"),
  compare_to_primary(full$ED_full,       primary$ED,      "ED: full cohort, ~ sex + length"),
  compare_to_primary(full$ED_full_disc,  primary$ED_disc, "ED: full cohort, ~ sex + length + disc (vs disc-adjusted primary)"),
  compare_to_primary(full$TE_sexmatched, primary$TE,      "EEM: sex-matched (IVP males vs IVV males)"),
  compare_to_primary(full$TE_full,       primary$TE,      "EEM: full cohort, ~ sex + length")
)
write_csv(comparison, file.path(paths$tables, "fullcohort_replication_summary.csv"))
print(comparison, width = Inf)

# ---- key genes across analyses ----------------------------------------------------------------
all_tabs <- c(list(ED_primary = primary$ED, ED_primary_disc = primary$ED_disc, TE_primary = primary$TE), full)
key_tab <- imap_dfr(all_tabs, function(t, nm) {
  t %>% filter(label %in% key_genes) %>%
    transmute(analysis = nm, tissue = ifelse(startsWith(nm, "ED"), "ED", "EEM"),
              gene = label, log2FC = log2FoldChange, padj)
})
write_csv(key_tab, file.path(paths$tables, "fullcohort_key_genes.csv"))

lvl <- c("ED_primary", "ED_primary_disc", "ED_sexmatched", "ED_full", "ED_full_disc",
         "TE_primary", "TE_sexmatched", "TE_full")
nice <- c(ED_primary = "Primary\n(same length)", ED_primary_disc = "Primary\n+ disc adj.",
          ED_sexmatched = "Sex-matched\n(males)", ED_full = "Full cohort", ED_full_disc = "Full cohort\n+ disc adj.",
          TE_primary = "Primary\n(same length)", TE_sexmatched = "Sex-matched\n(males)", TE_full = "Full cohort")
p <- key_tab %>%
  mutate(analysis = factor(analysis, lvl), gene = factor(gene, rev(key_genes)),
         sig = !is.na(padj) & padj < params$alpha,
         nlp = pmin(-log10(pmax(padj, 1e-12)), 12)) %>%
  ggplot(aes(analysis, gene)) +
  geom_point(aes(size = nlp, fill = log2FC, colour = sig), shape = 21, stroke = 0.9) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "grey95", high = "#D7191C", midpoint = 0,
                       limits = c(-4, 4), oob = scales::squish, name = "log2FC\n(IVP/IVV)") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"), name = "FDR < 0.05") +
  scale_size(range = c(1, 7), name = "-log10 FDR") +
  scale_x_discrete(labels = nice) +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "Key genes: primary analysis vs full-cohort checks",
       subtitle = "Same colour across columns = same direction of IVP effect") +
  theme(axis.text.x = element_text(size = 8))
save_plot(p, "fullcohort_key_genes", 10, 5.5)

# ---- key GSEA terms across analyses --------------------------------------------------------------
gsea_all_terms <- function(tab) {
  r <- tab %>% filter(!is.na(pvalue), !is.na(stat)) %>% distinct(entrez, .keep_all = TRUE) %>% arrange(desc(stat))
  g <- gseGO(setNames(r$stat, r$entrez), OrgDb = org.Bt.eg.db, keyType = "ENTREZID", ont = "BP",
             minGSSize = 10, maxGSSize = 500, eps = 0, pvalueCutoff = 1, pAdjustMethod = "BH",
             seed = TRUE, verbose = FALSE)
  as.data.frame(g)
}
term_tab <- imap_dfr(all_tabs, function(t, nm) {
  tissue <- ifelse(startsWith(nm, "ED"), "ED", "TE")
  gsea_all_terms(t) %>% filter(Description %in% key_terms[[tissue]]) %>%
    transmute(analysis = nm, tissue = ifelse(tissue == "TE", "EEM", "ED"), term = Description, NES, p.adjust)
})
write_csv(term_tab, file.path(paths$tables, "fullcohort_key_terms_GSEA.csv"))

p <- term_tab %>%
  mutate(analysis = factor(analysis, lvl), sig = p.adjust < params$alpha) %>%
  ggplot(aes(analysis, term)) +
  geom_point(aes(fill = NES, colour = sig), shape = 21, size = 6, stroke = 0.9) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "grey95", high = "#D7191C", midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"), name = "FDR < 0.05") +
  scale_x_discrete(labels = nice) +
  facet_grid(tissue ~ ., scales = "free", space = "free") +
  labs(x = NULL, y = NULL, title = "Key GSEA terms: primary analysis vs full-cohort checks") +
  theme(axis.text.x = element_text(size = 8))
save_plot(p, "fullcohort_key_terms", 10, 6)

message("08_full_cohort_checks.R done")
