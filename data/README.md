# Input data

The pipeline starts from a gene-level read-count matrix:

`data/D15_Embryos_RawCountData.txt`: tab-separated; first column = gene IDs as `GeneID:<Entrez ID>`; one column per
library named `<library>_1` (32 libraries; see `metadata/samples.csv`). Reads were aligned to ARS-UCD1.3 and counted
per gene.

The matrix is **not distributed** with this repository. Raw sequencing data are available from GEO, accession
**GSE278243** (Behrens et al., 2025). To obtain the gene-level count matrix used here, please contact the corresponding
author of Behrens et al. (2025).

The analysis uses the 17 same-length libraries flagged `same_length = TRUE` in `metadata/samples.csv` (IVPF10_EC is
excluded from the primary analysis). `08_full_cohort_checks.R` uses all libraries except IVPF13, whose two tissue
libraries disagree on inferred sex.

**Note on `metadata/samples.csv`:** all IVP embryos were transferred fresh (Behrens et al., 2025). The `cryopreserved`
column reproduces a label from the project sample sheet ("Vitrified" for IVP9, IVP10, IVP11 and IVP14; "Fresh" for all
other embryos). It is treated as a **processing-batch label**, not as cryopreservation. It is not used as a model
covariate (in the same-length set it is identical to production system), but it is shown in the PC–metadata
correlation (`python/report_signalling_pca.py`; Supplementary Figure S3 of the report), where it is associated with
several principal components.

**IVPF13:** its ED and EEM libraries disagree on inferred sex (ED male, EEM female), which indicates a probable sample
mix-up; it is excluded from the full-cohort checks.
