#!/usr/bin/env Rscript
# Building a TxDb from the annotation GTF, across Bioconductor versions.
#
# makeTxDbFromGFF() lived in GenomicFeatures until Bioconductor moved it into
# the txdbmaker package, where GenomicFeatures >= 1.61.1 makes the old entry
# point formally defunct rather than merely deprecated. Four scripts build a
# TxDb this way and each of them must run both on an environment that predates
# the move and on one that follows it, so the choice of entry point is made
# here once instead of at every call site.
#
# The caller is responsible for having loaded GenomicFeatures.

make_txdb_from_gff <- function(gtf_file, format = "gtf", ...) {
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdbmaker::makeTxDbFromGFF(gtf_file, format = format, ...)
  } else {
    makeTxDbFromGFF(gtf_file, format = format, ...)
  }
}
