# Install everything the pipeline needs (R >= 4.3, Bioconductor >= 3.18).
cran <- c("tidyverse", "matrixStats", "pheatmap", "ggrepel", "patchwork",
          "e1071", "randomForest", "quarto", "knitr", "BiocManager")
bioc <- c("DESeq2", "apeglm", "AnnotationDbi", "org.Bt.eg.db",
          "clusterProfiler", "enrichplot", "fgsea", "DOSE")

missing_cran <- setdiff(cran, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran)
missing_bioc <- setdiff(bioc, rownames(installed.packages()))
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
