# =============================================================================
# 10_sensitivity.R — embryo-influence and disc-marker sensitivity analyses
#
# 1. Leave-one-embryo-out (LOEO) influence: for each of the 8 primary embryos, drop BOTH of its
#    libraries and refit (a) ED ~ type, (b) ED ~ disc_score + type (original marker set, score
#    re-scaled), (c) EEM ~ type, (d) the paired tissue x group interaction. Record, for a FROZEN list
#    of GO terms and highlighted genes, the direction/size of the group contrast.
# 2. EEM without library IVP9_T1A: same outputs (stability, not DEG counts).
# 3. Disc-marker sensitivity: rebuild the disc-enriched marker set 4 times, each time omitting one
#    of the four primary IVV embryos; recompute the score and refit ED ~ disc_score + type.
#
# The frozen terms/genes were fixed before running this script (main themes of the report).
# Gene-set statistics use fgsea on the frozen GO terms only (gene sets from org.Bt.eg.db, GOALL),
# ranked by the DESeq2 Wald statistic; the reference column ("none") is computed the same way,
# so all columns are directly comparable. Unshrunken log2 fold changes are used throughout.
#
# Requires: 01, 02, 03 (results/objects). Outputs: results/tables/sens_*.csv
# Run time: about 10-20 minutes.
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages(library(fgsea))
set.seed(params$seed)

inp  <- readRDS(file.path(paths$objects, "inputs.rds"))
disc <- readRDS(file.path(paths$objects, "disc_score.rds"))

# ---- frozen terms and genes --------------------------------------------------------------------
frozen_terms <- c(
  "glycosylation", "cellular lipid metabolic process", "lipid biosynthetic process",        # ED metabolism
  "anatomical structure morphogenesis", "embryo development", "neuron differentiation",
  "canonical Wnt signaling pathway",                                                          # ED development
  "cell adhesion", "cell migration",                                                          # EEM adhesion
  "translation", "ribosome biogenesis", "cellular respiration",                               # EEM translation/energy
  "protein folding")                                                                          # shared heat shock
frozen_genes <- c("HSPA1A", "HSPA6", "DNAJB1", "HERPUD1",                                     # heat shock
                  "HEY1", "ID1", "FRAT2", "NKD1", "SMAD9", "SMAD7", "WNT10B",                   # signalling, kept after adj.
                  "WNT3", "WNT5A", "LEF1", "BMP4", "NOTCH1")                                  # signalling, lost after adj.

# GO IDs for the frozen terms, taken from the primary GSEA tables written by 04
go_tabs <- list.files(paths$tables, pattern = "^GSEA_GO_BP_.*\\.csv$", full.names = TRUE)
go_ids  <- map_dfr(go_tabs, ~ read_csv(.x, show_col_types = FALSE) %>% select(ID, Description)) %>%
  distinct() %>% filter(Description %in% frozen_terms)
missing <- setdiff(frozen_terms, go_ids$Description)
if (length(missing)) stop("GO IDs not found for: ", paste(missing, collapse = ", "))
gs <- suppressMessages(AnnotationDbi::select(org.Bt.eg.db, keys = go_ids$ID, keytype = "GOALL",
                                             columns = "ENTREZID"))
gene_sets <- split(gs$ENTREZID, gs$GOALL) %>% map(unique)
names(gene_sets) <- go_ids$Description[match(names(gene_sets), go_ids$ID)]

# ---- helpers ------------------------------------------------------------------------------------
fit <- function(meta, design, coef = "type_IVP_vs_IVV", contrast = NULL, full_mm = NULL) {
  meta <- droplevels(meta)
  cts  <- prefilter(inp$counts[, meta$sample])
  if (is.null(full_mm)) {
    dds <- DESeq(DESeqDataSetFromMatrix(cts, meta, design), quiet = TRUE)
  } else {
    # nested design with 3 vs 4 embryos: pass the model matrix (empty columns removed) directly
    dds <- DESeqDataSetFromMatrix(cts, meta, ~ 1)
    dds <- DESeq(dds, full = full_mm, quiet = TRUE)
  }
  res <- if (is.null(contrast)) results(dds, name = coef, alpha = params$alpha)
         else results(dds, contrast = contrast, alpha = params$alpha)
  annotate_genes(as.data.frame(res) %>% rownames_to_column("entrez"))
}

gsea_frozen <- function(tab) {
  r <- tab %>% filter(!is.na(stat)) %>% distinct(entrez, .keep_all = TRUE)
  st <- setNames(r$stat, r$entrez)
  set.seed(params$seed)
  fg <- fgsea(gene_sets, st, minSize = 10, maxSize = 2000, eps = 0)
  tibble(term = fg$pathway, NES = fg$NES, padj = fg$padj, size = fg$size)
}

gene_rows <- function(tab) {
  tab %>% filter(label %in% frozen_genes) %>% distinct(label, .keep_all = TRUE) %>%
    transmute(gene = label, log2FC = log2FoldChange, padj)
}

stability <- function(tab, ref) {
  j <- inner_join(select(ref, entrez, lfc_r = log2FoldChange, q_r = padj),
                  select(tab, entrez, lfc = log2FoldChange, q = padj), by = "entrez")
  p <- filter(j, !is.na(q_r), q_r < params$alpha, !is.na(lfc))
  tibble(DEGs = sum(tab$padj < params$alpha, na.rm = TRUE),
         primary_DEGs_compared = nrow(p),
         primary_DEGs_same_direction_pct = round(100 * mean(sign(p$lfc) == sign(p$lfc_r)), 1),
         primary_DEGs_still_FDR05 = sum(!is.na(p$q) & p$q < params$alpha),
         spearman_lfc_all_genes = round(cor(j$lfc, j$lfc_r, method = "spearman", use = "complete.obs"), 3),
         pearson_lfc_all_genes  = round(cor(j$lfc, j$lfc_r, use = "complete.obs"), 3))
}

interaction_fit <- function(meta) {
  meta <- meta %>% group_by(type) %>% mutate(embryo_n = factor(as.integer(factor(embryo)))) %>%
    ungroup() %>% as.data.frame()
  rownames(meta) <- meta$sample
  meta$type <- factor(meta$type, levels = c("IVV", "IVP")); meta$region <- factor(meta$region, levels = c("ED", "TE"))
  mm <- model.matrix(~ type + type:embryo_n + type:region, meta)
  mm <- mm[, colSums(mm != 0) > 0, drop = FALSE]                     # drop empty embryo columns (3 vs 4)
  colnames(mm) <- make.names(colnames(mm))                            # "typeIVP:regionTE" -> "typeIVP.regionTE"
  stopifnot(qr(mm)$rank == ncol(mm), all(c("typeIVP.regionTE", "typeIVV.regionTE") %in% colnames(mm)))
  fit(meta, ~ type + type:embryo_n + type:region, contrast = list("typeIVP.regionTE", "typeIVV.regionTE"),
      full_mm = mm)
}

prim <- analysis_samples(inp$samples)                                # 16 libraries
prim$disc_score <- disc$score$disc_score[match(prim$sample, disc$score$sample)]
embryos <- unique(prim$embryo)

run_set <- function(meta_all, omitted) {
  ed <- meta_all[meta_all$region == "ED", ]; te <- meta_all[meta_all$region == "TE", ]
  ed_adj <- ed; ed_adj$disc_score <- as.numeric(scale(ed_adj$disc_score))
  list(ED     = fit(ed, ~ type),
       ED_adj = fit(ed_adj, ~ disc_score + type),
       EEM    = fit(te, ~ type),
       interaction = interaction_fit(meta_all))
}

# ---- 1. leave-one-embryo-out ----------------------------------------------------------------------
message("Reference fits (all 8 embryos)")
ref <- run_set(prim, "none")
terms_out <- list(); genes_out <- list(); stab_out <- list()
collect <- function(res, omitted) {
  for (a in names(res)) {
    terms_out[[length(terms_out) + 1]] <<- gsea_frozen(res[[a]]) %>% mutate(analysis = a, omitted = omitted)
    genes_out[[length(genes_out) + 1]] <<- gene_rows(res[[a]]) %>% mutate(analysis = a, omitted = omitted)
    stab_out[[length(stab_out) + 1]]   <<- stability(res[[a]], ref[[a]]) %>% mutate(analysis = a, omitted = omitted)
  }
}
collect(ref, "none")
for (e in embryos) {
  message("Omitting embryo ", e)
  collect(run_set(prim[prim$embryo != e, ], e), e)
}

# ---- 2. EEM without IVP9_T1A ------------------------------------------------------------------------
te9 <- prim[prim$region == "TE" & prim$sample != "IVP9_T1A", ]
r9  <- fit(te9, ~ type)
terms_out[[length(terms_out) + 1]] <- gsea_frozen(r9) %>% mutate(analysis = "EEM", omitted = "IVP9_T1A library")
genes_out[[length(genes_out) + 1]] <- gene_rows(r9) %>% mutate(analysis = "EEM", omitted = "IVP9_T1A library")
stab_out[[length(stab_out) + 1]]   <- stability(r9, ref$EEM) %>% mutate(analysis = "EEM", omitted = "IVP9_T1A library")

write_csv(bind_rows(terms_out), file.path(paths$tables, "sens_loeo_terms.csv"))
write_csv(bind_rows(genes_out), file.path(paths$tables, "sens_loeo_genes.csv"))
write_csv(bind_rows(stab_out),  file.path(paths$tables, "sens_loeo_stability.csv"))

# ---- 3. disc-marker sensitivity (leave one IVV embryo out of the marker definition) ---------------
vsd_all <- assay(vst(DESeqDataSetFromMatrix(prefilter(inp$counts), inp$samples, ~ 1), blind = TRUE))
ivv_pairs <- prim %>% filter(type == "IVV") %>% select(embryo, region, sample) %>%
  pivot_wider(names_from = region, values_from = sample)
ed_prim <- prim[prim$region == "ED", ]
marker_out <- list(); marker_terms <- list()
for (k in c("none", ivv_pairs$embryo)) {
  keep <- if (k == "none") ivv_pairs else ivv_pairs[ivv_pairs$embryo != k, ]
  d    <- vsd_all[, keep$ED, drop = FALSE] - vsd_all[, keep$TE, drop = FALSE]
  mu   <- rowMeans(vsd_all[, c(keep$ED, keep$TE)])
  mk   <- rownames(vsd_all)[rowMins(d) > 1.5 & mu > 6]
  sc   <- colMeans(vsd_all[mk, ed_prim$sample])
  m    <- ed_prim; m$disc_score <- as.numeric(scale(sc))
  r    <- fit(m, ~ disc_score + type)
  ivp  <- sc[m$type == "IVP"]; ivv <- sc[m$type == "IVV"]
  marker_out[[length(marker_out) + 1]] <- tibble(
    IVV_embryo_left_out = k, n_markers = length(mk),
    overlap_with_original = length(intersect(mk, disc$markers)),
    score_r_with_original = round(cor(sc, disc$score$disc_score[match(ed_prim$sample, disc$score$sample)]), 3),
    mean_score_IVP = round(mean(ivp), 2), mean_score_IVV = round(mean(ivv), 2),
    IVP_lower_on_average = mean(ivp) < mean(ivv),
    n_IVP_below_lowest_IVV = sum(ivp < min(ivv)),
    adjusted_DEGs = sum(r$padj < params$alpha, na.rm = TRUE),
    adjusted_DEGs_up = sum(r$padj < params$alpha & r$log2FoldChange > 0, na.rm = TRUE),
    adjusted_DEGs_down = sum(r$padj < params$alpha & r$log2FoldChange < 0, na.rm = TRUE))
  marker_terms[[length(marker_terms) + 1]] <- gsea_frozen(r) %>% mutate(IVV_embryo_left_out = k)
}
write_csv(bind_rows(marker_out),   file.path(paths$tables, "sens_marker_loo_summary.csv"))
write_csv(bind_rows(marker_terms), file.path(paths$tables, "sens_marker_loo_terms.csv"))
print(bind_rows(marker_out), width = Inf)

message("10_sensitivity.R done")
