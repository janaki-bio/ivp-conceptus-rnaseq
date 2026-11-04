# =============================================================================
# 03_joint_model.R — one model for all 16 libraries, using the ED/TE pairing
#
# ED and TE come from the SAME embryo, so the libraries are paired. The DESeq2 design
# for "individuals nested within groups" (DESeq2 vignette) is:
#     ~ type + type:embryo_n + type:region
# which gives:
#   typeIVV.regionTE   TE vs ED within IVV embryos (paired)          -> region effect within IVV
#   typeIVP.regionTE   TE vs ED within IVP embryos (paired)          -> region effect within IVP
#   interaction        (TE-ED)_IVP - (TE-ED)_IVV
#                      = (IVP-IVV)_TE - (IVP-IVV)_ED                  -> "is the IVP effect different
#                                                                        in TE than in ED?"
# The interaction tests whether the IVP effect differs between tissues.
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages(library(apeglm))

inp  <- readRDS(file.path(paths$objects, "inputs.rds"))
meta <- analysis_samples(inp$samples)          # 16 libraries (IVPF10_EC excluded)
meta <- meta %>%
  group_by(type) %>%
  mutate(embryo_n = factor(as.integer(factor(embryo)))) %>%   # embryo index 1..4 within each type
  ungroup() %>% as.data.frame()
rownames(meta) <- meta$sample
stopifnot(all(table(meta$type, meta$embryo_n) == 2))         # every embryo contributes ED + TE

cts <- prefilter(inp$counts[, meta$sample])
dds <- DESeqDataSetFromMatrix(cts, colData = meta, design = ~ type + type:embryo_n + type:region)
dds <- DESeq(dds, quiet = TRUE)
print(resultsNames(dds))

res_ivv_region <- results(dds, name = "typeIVV.regionTE", alpha = params$alpha)
res_ivp_region <- results(dds, name = "typeIVP.regionTE", alpha = params$alpha)
res_inter      <- results(dds, contrast = list("typeIVP.regionTE", "typeIVV.regionTE"),
                          alpha = params$alpha)

tabs <- list(
  TE_vs_ED_within_IVV = tidy_results(res_ivv_region),
  TE_vs_ED_within_IVP = tidy_results(res_ivp_region),
  interaction_IVPeffect_TE_minus_ED = tidy_results(res_inter)
)
iwalk(tabs, write_de_tables)
saveRDS(list(dds = dds, results = list(ivv = res_ivv_region, ivp = res_ivp_region, inter = res_inter),
             tables = tabs, meta = meta),
        file.path(paths$objects, "joint_model.rds"))
write_csv(imap_dfr(tabs, de_summary), file.path(paths$tables, "joint_model_summary.csv"))

# ---- Is one region "more affected"? ----------------------------------------
# Compare the per-region IVP effects (02_de_by_region.R) on the SAME genes, using apeglm
# shrunken LFCs so that noisy low-count genes do not dominate.
ed <- readRDS(file.path(paths$objects, "de_ED_IVP_vs_IVV.rds"))$table
te <- readRDS(file.path(paths$objects, "de_TE_IVP_vs_IVV.rds"))$table
both <- inner_join(select(ed, entrez, label, lfc_ED = log2FC_shrunk, sig_ED = significant),
                   select(te, entrez, lfc_TE = log2FC_shrunk, sig_TE = significant),
                   by = "entrez") %>%
  filter(!is.na(lfc_ED), !is.na(lfc_TE))

wt <- wilcox.test(abs(both$lfc_ED), abs(both$lfc_TE), paired = TRUE)
magnitude <- tibble(
  genes_compared               = nrow(both),
  DEG_ED                       = sum(ed$significant),
  DEG_TE                       = sum(te$significant),
  DEG_both_regions             = sum(both$sig_ED & both$sig_TE),
  median_abs_shrunkLFC_ED      = median(abs(both$lfc_ED)),
  median_abs_shrunkLFC_TE      = median(abs(both$lfc_TE)),
  paired_wilcoxon_p            = wt$p.value,
  pearson_r_IVPeffect_ED_vs_TE = cor(both$lfc_ED, both$lfc_TE),
  interaction_genes_FDR05      = sum(tabs$interaction_IVPeffect_TE_minus_ED$significant)
)
write_csv(magnitude, file.path(paths$tables, "region_effect_magnitude.csv"))
write_csv(both, file.path(paths$tables, "IVPeffect_ED_vs_TE_per_gene.csv"))
print(t(magnitude))

p <- ggplot(both, aes(lfc_ED, lfc_TE)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(aes(colour = interaction(sig_ED, sig_TE)), size = 0.7, alpha = 0.6) +
  scale_colour_manual(values = c("FALSE.FALSE" = "grey80", "TRUE.FALSE" = "#1B9E77",
                                 "FALSE.TRUE" = "#7570B3", "TRUE.TRUE" = "#D95F02"),
                      labels = c("neither", "ED only", "TE only", "both"), name = "DEG in") +
  labs(x = "IVP vs IVV, ED (shrunken log2FC)", y = "IVP vs IVV, TE (shrunken log2FC)",
       title = "IVP effect in ED vs TE (per gene)",
       subtitle = sprintf("r = %.2f; interaction DEGs (FDR<%.2f) = %d",
                          magnitude$pearson_r_IVPeffect_ED_vs_TE, params$alpha,
                          magnitude$interaction_genes_FDR05))
save_plot(p, "IVPeffect_ED_vs_TE_scatter", 6.5, 5.5)

p <- both %>% pivot_longer(c(lfc_ED, lfc_TE), names_to = "region", values_to = "lfc") %>%
  mutate(region = sub("lfc_", "", region)) %>%
  ggplot(aes(abs(lfc), colour = region)) + stat_ecdf() +
  coord_cartesian(xlim = c(0, 2)) +
  labs(x = "|shrunken log2FC| IVP vs IVV", y = "Cumulative fraction of genes",
       title = "Magnitude of IVP effect per region")
save_plot(p, "IVPeffect_magnitude_ecdf", 5.5, 4)

message("03_joint_model.R done")
