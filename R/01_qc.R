# =============================================================================
# 01_qc.R — library QC, PCA, sample correlation, disc-content score
# Outputs: results/figures/qc_*, results/tables/qc_*, results/objects/inputs.rds
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(pheatmap); library(ggrepel) })

inp     <- load_inputs()
counts  <- inp$counts
samples <- inp$samples
saveRDS(inp, file.path(paths$objects, "inputs.rds"))

# ---- 1. library size and detected genes (all 32 libraries) ------------------
lib <- tibble(sample   = colnames(counts),
              reads_M  = colSums(counts) / 1e6,
              detected = colSums(counts > 0)) %>%
  left_join(samples, by = "sample")
write_csv(lib, file.path(paths$tables, "qc_library_sizes.csv"))

p <- ggplot(lib, aes(reorder(sample, reads_M), reads_M, fill = type)) +
  geom_col() + coord_flip() + scale_fill_manual(values = type_cols) +
  labs(x = NULL, y = "Assigned reads (millions)", title = "Library sizes (all 32 libraries)")
save_plot(p, "qc_library_sizes", 6, 7)

# ---- 2. PCA: all libraries, then the same-length analysis set --------------
dds_all <- DESeqDataSetFromMatrix(prefilter(counts), samples, design = ~ 1)
vsd_all <- assay(vst(dds_all, blind = TRUE))

plot_pca <- function(df, title, label = TRUE) {
  ve <- attr(df, "varexp")
  p <- ggplot(df, aes(PC1, PC2, colour = type, shape = region)) +
    geom_point(size = 3) +
    scale_colour_manual(values = type_cols) + scale_shape_manual(values = region_shps) +
    labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"), title = title)
  if (label) p <- p + geom_text_repel(aes(label = sample), size = 2.6, show.legend = FALSE)
  p
}
save_plot(plot_pca(pca_df(vsd_all, samples), "All 32 libraries"), "qc_pca_all_libraries", 8, 6)

# same-length set INCLUDING the excluded library, so its position is visible
same <- samples[samples$same_length, ]
save_plot(plot_pca(pca_df(vsd_all[, same$sample], same), "Same-length embryos (17 libraries)"),
          "qc_pca_same_length", 7, 5.5)

for (reg in c("ED", "TE")) {
  m <- same[same$region == reg, ]
  save_plot(plot_pca(pca_df(vsd_all[, m$sample], m), paste0("Same-length embryos, ", reg)),
            paste0("qc_pca_", reg), 6.5, 5)
}

# ---- 3. sample-to-sample correlation (flags outliers, e.g. IVP9_T1A) --------
top <- order(rowVars(vsd_all[, same$sample]), decreasing = TRUE)[1:2000]
cm  <- cor(vsd_all[top, same$sample])
ann <- same[, c("type", "region", "cryopreserved")]
png(file.path(paths$figures, "qc_sample_correlation.png"), width = 1600, height = 1400, res = 200)
pheatmap(cm, annotation_col = ann, annotation_row = ann, display_numbers = TRUE,
         number_format = "%.2f", fontsize_number = 5, main = "Pearson r, top-2000 variable genes (VST)")
dev.off()

# mean correlation of each library with the other libraries of the same region & type,
# using genes that are variable WITHIN that region (otherwise ED-vs-TE genes dominate)
outlier_tab <- map_dfr(c("ED", "TE"), function(reg) {
  s_reg <- same$sample[same$region == reg]
  g     <- order(rowVars(vsd_all[, s_reg]), decreasing = TRUE)[1:2000]
  cr    <- cor(vsd_all[g, s_reg])
  map_dfr(s_reg, function(s) {
    peers <- s_reg[same[s_reg, "type"] == same[s, "type"] & s_reg != s]
    tibble(sample = s, region = reg, mean_r_with_group = mean(cr[s, peers]))
  })
}) %>% arrange(mean_r_with_group)
write_csv(outlier_tab, file.path(paths$tables, "qc_within_group_correlation.csv"))

# ---- 4. disc-content score ---------------------------------------------------
# Disc-enriched genes are defined ONLY from the IVV embryos (paired ED vs TE), so the
# definition is independent of the IVP samples. The score is the mean VST expression of
# those genes; a low score in an ED library means the section carries proportionally less
# embryonic-disc tissue (smaller/less developed disc, or more trophectoderm in the section).
ivv_pairs <- same %>% filter(type == "IVV") %>% select(embryo, region, sample) %>%
  pivot_wider(names_from = region, values_from = sample)
d_ed_te <- vsd_all[, ivv_pairs$ED] - vsd_all[, ivv_pairs$TE]
mean_ivv <- rowMeans(vsd_all[, c(ivv_pairs$ED, ivv_pairs$TE)])
disc_markers <- rownames(vsd_all)[rowMins(d_ed_te) > 1.5 & mean_ivv > 6]
message(length(disc_markers), " disc-enriched marker genes")

disc_score <- tibble(sample = same$sample,
                     disc_score = colMeans(vsd_all[disc_markers, same$sample])) %>%
  left_join(same, by = "sample")
write_csv(annotate_genes(tibble(entrez = disc_markers)), file.path(paths$tables, "qc_disc_marker_genes.csv"))
write_csv(disc_score, file.path(paths$tables, "qc_disc_score.csv"))
saveRDS(list(markers = disc_markers, score = disc_score), file.path(paths$objects, "disc_score.rds"))

p <- ggplot(filter(disc_score, region == "ED"), aes(type, disc_score, colour = type)) +
  geom_point(size = 3, position = position_jitter(width = 0.08, seed = 1)) +
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = type_cols) +
  labs(x = NULL, y = "Mean VST of disc-enriched genes",
       title = "Disc-content score of ED sections",
       subtitle = "Markers defined from IVV ED vs TE pairs only") +
  theme(legend.position = "none")
save_plot(p, "qc_disc_score_ED", 5, 4.5)

message("01_qc.R done")
