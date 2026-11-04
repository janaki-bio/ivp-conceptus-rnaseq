# =============================================================================
# run_all.R — run the whole pipeline from raw counts to figures
# Open ivp_conceptus_rnaseq.Rproj in RStudio (so the working directory is this folder), then:
#   source("install_packages.R")   # first time only
#   source("run_all.R")
#   quarto::quarto_render("report.qmd")   # or click Render in RStudio
# =============================================================================
steps <- c("R/01_qc.R",
           "R/02_de_by_region.R",
           "R/03_joint_model.R",
           "R/04_enrichment.R",
           "R/05_figures.R",
           "R/06_classifier_loeo.R",   # optional; takes a few minutes (70 permutations)
           "R/07_classifier_by_region.R",  # optional; ~5-15 min
           "R/08_full_cohort_checks.R",     # robustness checks on all embryos
           "R/09_xlinked_check.R",          # X-linked gene contribution
           "R/10_sensitivity.R")            # embryo-influence and disc-marker sensitivity

for (s in steps) {
  message("\n==== ", s, " ====")
  t0 <- Sys.time()
  source(s, local = new.env())
  message(sprintf("     finished in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("\nAll done. Tables: results/tables  Figures: results/figures")
