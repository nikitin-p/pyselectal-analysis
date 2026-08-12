#!/usr/bin/env Rscript
# Stage 16 step 0: install the two packages Stage 16 needs that the system R
# library does not already carry, and verify the rest are present.
#
# clusterProfiler and org.Sc.sgd.db go into the system library, alongside
# seqArchR (Stage 15). The r_cage conda environment is not involved.
#
# Usage: Rscript scripts/09_go_enrichment/00_install.R

NEEDED  <- c("clusterProfiler", "org.Sc.sgd.db")
ASSUMED <- c("GenomicRanges", "AnnotationDbi", "enrichplot", "GOSemSim",
             "optparse", "ggplot2")

have <- function(p) requireNamespace(p, quietly = TRUE)

missing_assumed <- ASSUMED[!vapply(ASSUMED, have, logical(1))]
if (length(missing_assumed)) {
    stop("Packages expected to be present are missing: ",
         paste(missing_assumed, collapse = ", "),
         ". Install them before running Stage 16.")
}

if (!have("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org")

to_install <- NEEDED[!vapply(NEEDED, have, logical(1))]
if (length(to_install)) {
    message("Installing: ", paste(to_install, collapse = ", "))
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
} else {
    message("clusterProfiler and org.Sc.sgd.db are already installed.")
}

still_missing <- NEEDED[!vapply(NEEDED, have, logical(1))]
if (length(still_missing)) {
    stop("Installation failed for: ", paste(still_missing, collapse = ", "))
}

suppressMessages({
    library(clusterProfiler)
    library(org.Sc.sgd.db)
})
n_orf <- length(AnnotationDbi::keys(org.Sc.sgd.db, keytype = "ORF"))
message(sprintf("org.Sc.sgd.db carries %d ORF keys.", n_orf))
if (n_orf < 5000) stop("org.Sc.sgd.db returned too few ORF keys; the install is broken.")

message("Stage 16 dependencies are ready.")
sessionInfo()
