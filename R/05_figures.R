# =============================================================================
# 05_figures.R — publication figures
#   volcano + MA (all tested genes), top-DEG heatmaps, GO dot/lollipop plots,
#   gene-level expression plots for key pathways (one panel per gene, one point per embryo)
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(ggrepel); library(pheatmap); library(clusterProfiler); library(enrichplot); library(patchwork)
})

# ---- volcano & MA -------------------------------------------------------------
volcano <- function(tab, title, n_label = 15) {
  tab <- filter(tab, !is.na(padj))
  lab <- tab %>% filter(significant) %>% arrange(padj) %>% head(n_label)
  ggplot(tab, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(colour = direction), size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = c(up = "#D7191C", down = "#2C7BB6", ns = "grey75"),
                        breaks = c("up", "down"), name = NULL) +
    geom_hline(yintercept = -log10(params$alpha), linetype = 2, colour = "grey40") +
    geom_text_repel(data = lab, aes(label = label), size = 2.8, max.overlaps = 30) +
    labs(x = "log2 fold change (IVP / IVV)", y = "-log10 FDR", title = title,
         subtitle = sprintf("%d up, %d down (FDR < %.2f), %d genes tested",
                            sum(tab$direction == "up"), sum(tab$direction == "down"),
                            params$alpha, nrow(tab)))
}
ma_plot <- function(tab, title) {
  ggplot(filter(tab, !is.na(log2FC_shrunk)), aes(log10(baseMean + 1), log2FC_shrunk)) +
    geom_point(aes(colour = direction), size = 0.6, alpha = 0.6) +
    scale_colour_manual(values = c(up = "#D7191C", down = "#2C7BB6", ns = "grey75"),
                        breaks = c("up", "down"), name = NULL) +
    geom_hline(yintercept = 0) +
    labs(x = "log10 mean normalised count", y = "shrunken log2FC (apeglm)", title = title)
}

for (tag in c("ED_IVP_vs_IVV", "TE_IVP_vs_IVV", "ED_IVP_vs_IVV_discAdjusted")) {
  obj <- readRDS(file.path(paths$objects, paste0("de_", tag, ".rds")))
  save_plot(volcano(obj$table, gsub("_", " ", tag)) | ma_plot(obj$table, "MA plot"),
            paste0("volcano_MA_", tag), 12, 5)

  # heatmap of the top 50 DEGs (z-scored VST of that region's samples)
  sig <- obj$table %>% filter(significant) %>% head(50)
  if (nrow(sig) >= 2) {
    vsd <- assay(vst(obj$dds, blind = FALSE))
    z   <- t(scale(t(vsd[sig$entrez, , drop = FALSE])))
    rownames(z) <- make.unique(sig$label)
    ann <- data.frame(type = obj$meta$type, row.names = obj$meta$sample)
    png(file.path(paths$figures, paste0("heatmap_top50_", tag, ".png")),
        width = 1400, height = 2000, res = 200)
    pheatmap(z, annotation_col = ann, annotation_colors = list(type = type_cols),
             fontsize_row = 6, main = paste("Top 50 DEGs,", tag))
    dev.off()
  }
}

# ---- GO plots -------------------------------------------------------------------
lollipop <- function(ego, title, n = 15) {
  df <- as.data.frame(ego)
  if (nrow(df) == 0) return(NULL)
  df <- df %>% arrange(p.adjust) %>% head(n) %>%
    mutate(Description = fct_reorder(Description, -log10(p.adjust)))
  ggplot(df, aes(-log10(p.adjust), Description)) +
    geom_segment(aes(x = 0, xend = -log10(p.adjust), yend = Description), colour = "grey60") +
    geom_point(aes(size = Count, colour = -log10(p.adjust))) +
    scale_colour_gradient(low = "#2C7BB6", high = "#D7191C", name = "-log10 FDR") +
    labs(x = "-log10 FDR (BH)", y = NULL, size = "Genes", title = title)
}

for (tag in c("ED_IVP_vs_IVV", "TE_IVP_vs_IVV", "ED_IVP_vs_IVV_discAdjusted")) {
  enr <- readRDS(file.path(paths$objects, paste0("enrich_", tag, ".rds")))
  for (dir in c("up", "down")) {
    ego <- enr[[paste0("ora_", dir)]]
    p <- if (!is.null(ego)) lollipop(ego, paste("GO BP,", dir, "in IVP,", tag)) else NULL
    if (!is.null(p)) save_plot(p, paste0("GO_lollipop_", dir, "_", tag), 8, 5)
  }
  if (nrow(as.data.frame(enr$gsea)) > 0) {
    p <- dotplot(enr$gsea, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign) +
      ggtitle(paste("GSEA GO BP,", tag))
    save_plot(p, paste0("GSEA_dotplot_", tag), 10, 6)
  }
}

cc <- readRDS(file.path(paths$objects, "enrich_compareCluster.rds"))
if (nrow(as.data.frame(cc)) > 0) {
  save_plot(dotplot(cc, showCategory = 6) + ggtitle("GO BP enrichment of DEG sets"),
            "GO_compareCluster_ED_TE", 9, 8)
}

# ---- gene-level plots for key developmental pathways ---------------------------
# One panel per gene, one point per sample, IVV vs IVP; VST values are already on a
# log2-like scale, so they are plotted as they are (no exp()).
key_terms <- c("GO:0030509" = "BMP signaling pathway",
               "GO:0060070" = "canonical Wnt signaling pathway",
               "GO:0007179" = "TGF-beta receptor signaling pathway",
               "GO:0007219" = "Notch signaling pathway")
go2genes <- suppressMessages(AnnotationDbi::select(org.Bt.eg.db, keys = names(key_terms),
                                                  keytype = "GOALL", columns = "ENTREZID")) %>%
  distinct(GOALL, ENTREZID)

for (tag in c("ED_IVP_vs_IVV", "TE_IVP_vs_IVV")) {
  obj <- readRDS(file.path(paths$objects, paste0("de_", tag, ".rds")))
  vsd <- assay(vst(obj$dds, blind = FALSE))
  for (go in names(key_terms)) {
    genes <- obj$table %>% filter(significant, entrez %in% go2genes$ENTREZID[go2genes$GOALL == go]) %>%
      arrange(padj) %>% head(12)
    if (nrow(genes) == 0) next
    long <- as.data.frame(vsd[genes$entrez, , drop = FALSE]) %>% rownames_to_column("entrez") %>%
      pivot_longer(-entrez, names_to = "sample", values_to = "vst") %>%
      left_join(select(genes, entrez, label, log2FoldChange, padj), by = "entrez") %>%
      left_join(select(obj$meta, sample, type), by = "sample") %>%
      mutate(panel = sprintf("%s\nLFC %.1f, FDR %.1e", label, log2FoldChange, padj))
    p <- ggplot(long, aes(type, vst, colour = type)) +
      geom_boxplot(outlier.shape = NA, width = 0.5) +
      geom_point(position = position_jitter(width = 0.1, seed = 1), size = 1.8) +
      facet_wrap(~ panel, scales = "free_y") +
      scale_colour_manual(values = type_cols) +
      labs(x = NULL, y = "VST expression (log2 scale)",
           title = paste0(key_terms[[go]], " (", go, ") — DEGs, ", tag)) +
      theme(legend.position = "none", strip.text = element_text(size = 7))
    save_plot(p, paste0("genes_", sub(":", "", go), "_", tag), 9, 7)
  }
}

message("05_figures.R done")
