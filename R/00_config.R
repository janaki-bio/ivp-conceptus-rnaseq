# =============================================================================
# 00_config.R — paths, parameters and shared helpers
# Day-15 bovine conceptus RNA-seq: IVP vs IVV, embryonic disc (ED) vs trophectoderm (TE)
# Every other script starts with source("R/00_config.R").
# Working directory must be the repository root (open ivp_conceptus_rnaseq.Rproj).
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Bt.eg.db)
  library(matrixStats)
  library(tidyverse)
})
# Bioconductor packages mask several dplyr verbs (AnnotationDbi::select, S4Vectors::rename, ...).
# Pin the dplyr versions in the global environment so they always win.
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
first     <- dplyr::first
rowVars   <- matrixStats::rowVars
rowMins   <- matrixStats::rowMins

# ---- paths ------------------------------------------------------------------
paths <- list(
  counts  = "data/D15_Embryos_RawCountData.txt",  # gene-level count matrix (32 libraries); see data/README.md
  samples = "metadata/samples.csv",                     # clean sample sheet (one row per library)
  results = "results",
  figures = "results/figures",
  tables  = "results/tables",
  objects = "results/objects"
)
for (p in paths[c("results", "figures", "tables", "objects")]) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

# ---- analysis parameters ----------------------------------------------------
params <- list(
  alpha       = 0.05,          # FDR (BH) threshold for calling DEGs
  lfc_display = 1,             # |log2FC| used ONLY to highlight genes in plots, not for calling DEGs
  min_count   = 10,            # pre-filter: keep genes with >= min_count reads ...
  min_samples = 4,             # ... in >= min_samples libraries (= smallest group size)
  exclude     = c("IVPF10_EC"),# IVP disc library without a matching EEM library; excluded from the
                               # primary analysis
  seed        = 2024
)
set.seed(params$seed)

# ---- helpers ----------------------------------------------------------------

#' Read the count matrix and the sample sheet, aligned and validated.
load_inputs <- function() {
  counts <- read.delim(paths$counts, row.names = 1, check.names = FALSE)
  colnames(counts) <- sub("_1$", "", colnames(counts))        # "IVP10_EC_1" -> "IVP10_EC"
  rownames(counts) <- sub("^GeneID:", "", rownames(counts))   # "GeneID:540784" -> "540784" (Entrez)
  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"

  samples <- read.csv(paths$samples, stringsAsFactors = FALSE)
  stopifnot(setequal(samples$sample, colnames(counts)))
  samples <- samples[match(colnames(counts), samples$sample), ]
  rownames(samples) <- samples$sample
  list(counts = counts, samples = samples)
}

#' Samples used in the primary analysis: same-length embryos minus params$exclude.
analysis_samples <- function(samples, region = NULL, exclude = params$exclude) {
  keep <- samples$same_length & !(samples$sample %in% exclude)
  if (!is.null(region)) keep <- keep & samples$region %in% region
  meta <- samples[keep, , drop = FALSE]
  meta$type   <- factor(meta$type,   levels = c("IVV", "IVP"))   # IVV = reference (control)
  meta$region <- factor(meta$region, levels = c("ED", "TE"))
  meta
}

#' Pre-filter low counts (DESeq2 still does independent filtering on top of this).
prefilter <- function(cts) {
  cts[rowSums(cts >= params$min_count) >= params$min_samples, , drop = FALSE]
}

#' Add gene symbol / name from org.Bt.eg.db (Entrez IDs as keys).
annotate_genes <- function(tab, id_col = "entrez") {
  ids <- as.character(tab[[id_col]])
  safe_map <- function(column) {
    out <- suppressMessages(AnnotationDbi::mapIds(org.Bt.eg.db, keys = ids, keytype = "ENTREZID",
                                                  column = column, multiVals = "first"))
    unname(out[ids])
  }
  tab$symbol   <- safe_map("SYMBOL")
  tab$genename <- safe_map("GENENAME")
  tab$label    <- ifelse(is.na(tab$symbol), paste0("LOC", ids), tab$symbol)
  tab
}

#' DESeq2 results -> tidy, annotated data frame.
tidy_results <- function(res, shrunk = NULL) {
  tab <- as.data.frame(res) %>% rownames_to_column("entrez")
  if (!is.null(shrunk)) {
    tab$log2FC_shrunk <- shrunk$log2FoldChange[match(tab$entrez, rownames(shrunk))]
  }
  tab <- annotate_genes(tab)
  tab$significant <- !is.na(tab$padj) & tab$padj < params$alpha
  tab$direction   <- case_when(!tab$significant ~ "ns",
                               tab$log2FoldChange > 0 ~ "up",
                               TRUE ~ "down")
  tab %>% arrange(padj, pvalue)
}

#' Write the full table plus significant / up / down subsets.
write_de_tables <- function(tab, tag) {
  write_csv(tab, file.path(paths$tables, paste0("DE_", tag, "_all_genes.csv")))
  write_csv(filter(tab, significant), file.path(paths$tables, paste0("DE_", tag, "_significant.csv")))
  write_csv(filter(tab, direction == "up"),   file.path(paths$tables, paste0("DE_", tag, "_up.csv")))
  write_csv(filter(tab, direction == "down"), file.path(paths$tables, paste0("DE_", tag, "_down.csv")))
  invisible(tab)
}

#' One-line DEG summary.
de_summary <- function(tab, tag) {
  tibble(analysis = tag,
         genes_tested = sum(!is.na(tab$padj)),
         DEG = sum(tab$significant),
         up = sum(tab$direction == "up"),
         down = sum(tab$direction == "down"),
         median_abs_lfc_DEG = median(abs(tab$log2FoldChange[tab$significant])),
         median_abs_shrunk_lfc_all = if ("log2FC_shrunk" %in% names(tab))
           median(abs(tab$log2FC_shrunk), na.rm = TRUE) else NA_real_)
}

#' PCA on the top-variance genes of a (VST) matrix.
pca_df <- function(mat, meta, ntop = 500) {
  rv  <- rowVars(mat)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pc  <- prcomp(t(mat[sel, ]))
  ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df  <- cbind(meta[colnames(mat), , drop = FALSE], as.data.frame(pc$x[, 1:min(4, ncol(pc$x))]))
  attr(df, "varexp") <- ve
  df
}

#' Save a ggplot as both PDF (for publication) and PNG (for quick viewing).
save_plot <- function(p, name, width = 7, height = 5) {
  ggsave(file.path(paths$figures, paste0(name, ".pdf")), p, width = width, height = height)
  ggsave(file.path(paths$figures, paste0(name, ".png")), p, width = width, height = height, dpi = 200)
  invisible(p)
}

theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank()))
type_cols   <- c(IVV = "#2C7BB6", IVP = "#D7191C")
region_shps <- c(ED = 16, TE = 17)
