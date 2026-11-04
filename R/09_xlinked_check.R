# =============================================================================
# 09_xlinked_check.R — how much do X-linked genes contribute to the IVP vs IVV results?
#
# Why: in the same-length set all IVP embryos are male and all IVV embryos female, so
# X-linked genes (incomplete X-inactivation / dosage) can appear "IVP-regulated" purely because of sex.
#
# 1. Count X-linked genes among tested genes and DEGs per tissue; Fisher test for X enrichment.
# 2. Supplementary table of X-linked DEGs.
# 3. Re-run ORA and GSEA without X-linked genes and compare the key terms with the primary results.
#
# Requires: 02 (de_*.rds objects). Outputs: results/tables/xlinked_*.csv, results/figures/xlinked_key_terms.png
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages(library(clusterProfiler))

tags <- c(ED = "ED_IVP_vs_IVV", EEM = "TE_IVP_vs_IVV", ED_adj = "ED_IVP_vs_IVV_discAdjusted")
key_terms <- list(
  ED     = c("glycosylation", "cellular lipid metabolic process", "lipid biosynthetic process",
             "anatomical structure morphogenesis", "neuron differentiation", "cell adhesion",
             "canonical Wnt signaling pathway"),
  EEM    = c("translation", "ribosome biogenesis", "protein folding", "cellular respiration",
             "cell migration", "cell adhesion"),
  ED_adj = c("glycosylation", "cellular lipid metabolic process", "anatomical structure morphogenesis",
             "neuron differentiation", "cell adhesion")
)

add_chr <- function(tab) {
  chr <- suppressMessages(AnnotationDbi::mapIds(org.Bt.eg.db, keys = tab$entrez, keytype = "ENTREZID",
                                                column = "CHR", multiVals = "first"))
  tab$chr <- unname(chr[tab$entrez])
  tab
}

tabs <- map(tags, ~ add_chr(readRDS(file.path(paths$objects, paste0("de_", .x, ".rds")))$table))

# ---- 1. counts and Fisher test ------------------------------------------------------------
summary_tab <- imap_dfr(tabs, function(t, nm) {
  t <- filter(t, !is.na(padj))
  isX <- !is.na(t$chr) & t$chr == "X"
  ft  <- fisher.test(table(factor(t$significant, c(TRUE, FALSE)), factor(isX, c(TRUE, FALSE))))
  tibble(analysis = nm,
         genes_tested = nrow(t), X_tested = sum(isX), pct_X_tested = round(100 * mean(isX), 2),
         DEG = sum(t$significant), X_DEG = sum(isX & t$significant),
         pct_X_DEG = round(100 * mean(isX[t$significant]), 2),
         X_up = sum(isX & t$direction == "up"), X_down = sum(isX & t$direction == "down"),
         chr_unknown_pct = round(100 * mean(is.na(t$chr)), 1),
         fisher_OR_X_among_DEG = unname(ft$estimate), fisher_p = ft$p.value)
})
write_csv(summary_tab, file.path(paths$tables, "xlinked_summary.csv"))
print(summary_tab, width = Inf)

# ---- 2. supplementary table of X-linked DEGs ------------------------------------------------
xdeg <- imap_dfr(tabs[c("ED", "EEM")], function(t, nm) {
  t %>% filter(significant, chr == "X") %>%
    transmute(tissue = nm, entrez, symbol = label, genename, log2FoldChange = round(log2FoldChange, 2),
              padj = signif(padj, 3), direction)
}) %>% arrange(tissue, padj)
write_csv(xdeg, file.path(paths$tables, "xlinked_DEGs.csv"))

# ---- 3. enrichment with and without X-linked genes ---------------------------------------------
gsea_terms <- function(t, terms) {
  r <- t %>% filter(!is.na(pvalue), !is.na(stat)) %>% distinct(entrez, .keep_all = TRUE) %>% arrange(desc(stat))
  g <- gseGO(setNames(r$stat, r$entrez), OrgDb = org.Bt.eg.db, keyType = "ENTREZID", ont = "BP",
             minGSSize = 10, maxGSSize = 500, eps = 0, pvalueCutoff = 1, pAdjustMethod = "BH",
             seed = TRUE, verbose = FALSE)
  as.data.frame(g) %>% filter(Description %in% terms) %>% select(Description, NES, p.adjust)
}
ora_terms <- function(t, terms) {
  uni <- t$entrez[!is.na(t$padj)]
  map_dfr(c("up", "down"), function(d) {
    e <- enrichGO(t$entrez[t$direction == d], universe = uni, OrgDb = org.Bt.eg.db, keyType = "ENTREZID",
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                  minGSSize = 10, maxGSSize = 500)
    if (is.null(e)) return(tibble(Description = character(), direction = character(),
                                  GeneRatio = character(), ora_padj = numeric()))
    as.data.frame(e) %>% filter(Description %in% terms) %>%
      transmute(Description, direction = d, GeneRatio, ora_padj = p.adjust)
  })
}

comparison <- imap_dfr(tabs, function(t, nm) {
  noX <- filter(t, is.na(chr) | chr != "X")
  bind_rows(
    gsea_terms(t,   key_terms[[nm]]) %>% mutate(genes = "all genes"),
    gsea_terms(noX, key_terms[[nm]]) %>% mutate(genes = "X-linked removed")
  ) %>% mutate(analysis = nm) %>%
    left_join(bind_rows(ora_terms(t,   key_terms[[nm]]) %>% mutate(genes = "all genes"),
                        ora_terms(noX, key_terms[[nm]]) %>% mutate(genes = "X-linked removed")),
              by = c("Description", "genes"), relationship = "many-to-many")
})
write_csv(comparison, file.path(paths$tables, "xlinked_key_terms_comparison.csv"))

wide <- comparison %>% distinct(analysis, Description, genes, NES, p.adjust) %>%
  pivot_wider(names_from = genes, values_from = c(NES, p.adjust))
print(wide, n = Inf, width = Inf)

p <- comparison %>% distinct(analysis, Description, genes, NES, p.adjust) %>%
  mutate(sig = p.adjust < params$alpha,
         analysis = recode(analysis, ED = "ED", EEM = "EEM", ED_adj = "ED (disc-score adjusted)")) %>%
  ggplot(aes(NES, Description, shape = genes, colour = genes)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(aes(alpha = sig), size = 3) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), name = "FDR < 0.05") +
  scale_colour_manual(values = c("all genes" = "#2a78d6", "X-linked removed" = "#eb6834"), name = NULL) +
  scale_shape_manual(values = c("all genes" = 16, "X-linked removed" = 17), name = NULL) +
  facet_grid(analysis ~ ., scales = "free_y", space = "free_y") +
  labs(x = "GSEA NES (IVP vs IVV)", y = NULL, title = "Key GO terms with and without X-linked genes")
save_plot(p, "xlinked_key_terms", 8.5, 6.5)

message("09_xlinked_check.R done")
