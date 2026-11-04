# =============================================================================
# 04_enrichment.R — GO enrichment done the standard way
#
# ORA  (over-representation; what DAVID does):
#       enrichGO on up- and down-regulated DEGs separately,
#       universe = all genes actually tested in that comparison (NOT the whole genome),
#       BH-adjusted p-values, redundant terms collapsed with simplify().
# GSEA (rank-based):
#       gseGO on ALL tested genes ranked by the DESeq2 Wald statistic.
#       (GSEA requires the full ranked list; it is not run on pre-filtered DEG lists.)
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(clusterProfiler); library(enrichplot) })

ora_one <- function(genes, universe) {
  if (length(genes) < 5) return(NULL)
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Bt.eg.db,
                  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  minGSSize = 10, maxGSSize = 500, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(ego)
  clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

run_enrichment <- function(tag, do_ora = TRUE) {
  tab <- readRDS(file.path(paths$objects, paste0("de_", tag, ".rds")))$table
  out <- list()

  if (do_ora) {
    universe <- tab$entrez[!is.na(tab$padj)]
    out$ora_up   <- ora_one(tab$entrez[tab$direction == "up"],   universe)
    out$ora_down <- ora_one(tab$entrez[tab$direction == "down"], universe)
    for (nm in c("ora_up", "ora_down")) {
      if (!is.null(out[[nm]])) {
        write_csv(as.data.frame(out[[nm]]),
                  file.path(paths$tables, paste0("GO_", nm, "_", tag, ".csv")))
      }
    }
  }

  ranked <- tab %>% filter(!is.na(pvalue), !is.na(stat)) %>%
    distinct(entrez, .keep_all = TRUE) %>% arrange(desc(stat))
  gene_list <- setNames(ranked$stat, ranked$entrez)
  out$gsea <- gseGO(geneList = gene_list, OrgDb = org.Bt.eg.db, keyType = "ENTREZID",
                    ont = "BP", minGSSize = 10, maxGSSize = 500, eps = 0,
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, verbose = FALSE)
  out$gsea <- setReadable(out$gsea, OrgDb = org.Bt.eg.db, keyType = "ENTREZID")
  write_csv(as.data.frame(out$gsea), file.path(paths$tables, paste0("GSEA_GO_BP_", tag, ".csv")))

  # KEGG needs internet access to rest.kegg.jp at run time; skipped gracefully if unavailable
  out$gsea_kegg <- tryCatch(
    gseKEGG(geneList = gene_list, organism = "bta", keyType = "ncbi-geneid",
            minGSSize = 10, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
            pAdjustMethod = "BH", seed = TRUE, verbose = FALSE),
    error = function(e) { message("KEGG skipped for ", tag, ": ", conditionMessage(e)); NULL })
  if (!is.null(out$gsea_kegg)) {
    write_csv(as.data.frame(out$gsea_kegg), file.path(paths$tables, paste0("GSEA_KEGG_", tag, ".csv")))
  }

  saveRDS(out, file.path(paths$objects, paste0("enrich_", tag, ".rds")))
  message(sprintf("%-28s ORA up %s | ORA down %s | GSEA %d terms", tag,
                  if (is.null(out$ora_up)) "NA" else nrow(as.data.frame(out$ora_up)),
                  if (is.null(out$ora_down)) "NA" else nrow(as.data.frame(out$ora_down)),
                  nrow(as.data.frame(out$gsea))))
  out
}

enr <- list(
  ED          = run_enrichment("ED_IVP_vs_IVV"),
  TE          = run_enrichment("TE_IVP_vs_IVV"),
  ED_adjusted = run_enrichment("ED_IVP_vs_IVV_discAdjusted")
)

# Interaction (IVP effect TE vs ED): GSEA only, from the joint model
inter <- readRDS(file.path(paths$objects, "joint_model.rds"))$tables$interaction_IVPeffect_TE_minus_ED
saveRDS(list(table = inter), file.path(paths$objects, "de_interaction.rds"))
enr$interaction <- run_enrichment("interaction", do_ora = FALSE)

# ---- side-by-side comparison of the four DEG sets -----------------------------
ed <- readRDS(file.path(paths$objects, "de_ED_IVP_vs_IVV.rds"))$table
te <- readRDS(file.path(paths$objects, "de_TE_IVP_vs_IVV.rds"))$table
clusters <- list(ED_up   = ed$entrez[ed$direction == "up"],
                 ED_down = ed$entrez[ed$direction == "down"],
                 TE_up   = te$entrez[te$direction == "up"],
                 TE_down = te$entrez[te$direction == "down"])
universe_union <- union(ed$entrez[!is.na(ed$padj)], te$entrez[!is.na(te$padj)])
cc <- compareCluster(geneClusters = clusters, fun = "enrichGO", universe = universe_union,
                     OrgDb = org.Bt.eg.db, keyType = "ENTREZID", ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
saveRDS(cc, file.path(paths$objects, "enrich_compareCluster.rds"))
write_csv(as.data.frame(cc), file.path(paths$tables, "GO_compareCluster_ED_TE_up_down.csv"))

message("04_enrichment.R done")
