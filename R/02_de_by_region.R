# =============================================================================
# 02_de_by_region.R — IVP vs IVV within each region (primary analysis + sensitivity)
#
# Primary:      ED: 4 IVP vs 4 IVV   |   TE: 4 IVP vs 4 IVV      design ~ type
# Sensitivity:  (a) TE without the outlier library IVP9_T1A
#               (b) ED adjusted for disc-content score            design ~ disc_score + type
#               (c) ED including the fresh IVP library IVPF10_EC
#
# NOTE: in this design every IVP embryo is male and every IVV embryo female, and each group comes
# from two dams, so "type" is confounded with embryo sex and dam (see README).
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages(library(apeglm))

inp  <- readRDS(file.path(paths$objects, "inputs.rds"))
disc <- readRDS(file.path(paths$objects, "disc_score.rds"))

run_de <- function(reg, tag, exclude = params$exclude, adjust_disc = FALSE) {
  meta <- analysis_samples(inp$samples, region = reg, exclude = exclude)
  cts  <- prefilter(inp$counts[, meta$sample])

  if (adjust_disc) {
    meta$disc_score <- disc$score$disc_score[match(meta$sample, disc$score$sample)]
    meta$disc_score <- as.numeric(scale(meta$disc_score))
    design <- ~ disc_score + type
  } else {
    design <- ~ type
  }

  dds <- DESeqDataSetFromMatrix(cts, colData = meta, design = design)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = "type_IVP_vs_IVV", alpha = params$alpha)
  shr <- lfcShrink(dds, coef = "type_IVP_vs_IVV", type = "apeglm", quiet = TRUE)

  tab <- tidy_results(res, shr)
  write_de_tables(tab, tag)
  saveRDS(list(dds = dds, res = res, shrunk = shr, table = tab, meta = meta),
          file.path(paths$objects, paste0("de_", tag, ".rds")))
  message(sprintf("%-28s %5d DEGs (%d up, %d down)", tag, sum(tab$significant),
                  sum(tab$direction == "up"), sum(tab$direction == "down")))
  tab
}

tabs <- list(
  ED                 = run_de("ED", "ED_IVP_vs_IVV"),
  TE                 = run_de("TE", "TE_IVP_vs_IVV"),
  TE_no_IVP9         = run_de("TE", "TE_IVP_vs_IVV_noIVP9T1A", exclude = c(params$exclude, "IVP9_T1A")),
  ED_disc_adjusted   = run_de("ED", "ED_IVP_vs_IVV_discAdjusted", adjust_disc = TRUE),
  ED_with_IVPF10     = run_de("ED", "ED_IVP_vs_IVV_withIVPF10", exclude = character(0))
)

summary_tab <- imap_dfr(tabs, de_summary)
write_csv(summary_tab, file.path(paths$tables, "DE_summary.csv"))
print(summary_tab)

# ---- how much of the ED signal tracks disc identity? -------------------------
# Correlate the IVP-vs-IVV fold change in ED with the ED-vs-TE fold change in IVV embryos.
# A strong negative correlation means "down in IVP ED" largely = "disc-specific gene".
ivv   <- analysis_samples(inp$samples) %>% filter(type == "IVV")
rownames(ivv) <- ivv$sample
dds_r <- DESeqDataSetFromMatrix(prefilter(inp$counts[, ivv$sample]), ivv, design = ~ embryo + region)
dds_r <- DESeq(dds_r, quiet = TRUE)
res_r <- results(dds_r, name = "region_TE_vs_ED")
region_lfc <- tibble(entrez = rownames(res_r), ED_vs_TE_lfc_in_IVV = -res_r$log2FoldChange)

comp <- tabs$ED %>% select(entrez, label, log2FoldChange, significant, direction) %>%
  inner_join(region_lfc, by = "entrez")
rho  <- cor(comp$log2FoldChange[comp$significant], comp$ED_vs_TE_lfc_in_IVV[comp$significant],
            method = "spearman", use = "complete.obs")
composition_check <- comp %>% filter(significant) %>% group_by(direction) %>%
  summarise(n = n(),
            median_ED_vs_TE_lfc = median(ED_vs_TE_lfc_in_IVV, na.rm = TRUE),
            pct_disc_enriched = 100 * mean(ED_vs_TE_lfc_in_IVV > 1, na.rm = TRUE)) %>%
  mutate(spearman_rho_all_DEG = rho)
write_csv(composition_check, file.path(paths$tables, "ED_composition_check.csv"))
write_csv(comp, file.path(paths$tables, "ED_IVPeffect_vs_region_identity.csv"))
print(composition_check)

message("02_de_by_region.R done")
