#!/usr/bin/env Rscript
# Stage 16 step 3: GO Biological Process enrichment by initiator dinucleotide.
#
# The question concerns the non-canonical RR initiators AA and AG, but all
# sixteen dinucleotides are tested against the same universe. The canonical YR
# initiators (CA, CG, TA, TG) therefore come along at no cost and act as the
# built-in contrast, and no dinucleotide is special-cased in the code.
#
# A gene is assigned to the initiator of its dominant TSS, which is the assigned
# TSS with the highest mean_total. Using every TSS would place well-expressed
# genes in several classes at once and blur exactly the contrast being tested.
#
# Input : results/go/tss_gene_assignment.tsv, results/go/universe.txt
# Output: results/go/go_initiator.tsv, results/go/go_initiator_set_sizes.tsv

suppressMessages(library(optparse))

DINUCS <- as.vector(t(outer(c("A", "C", "G", "T"), c("A", "C", "G", "T"), paste0)))

# a : the assignment table.
#
# Returns the named list of 16 character vectors, as before, with a "stats"
# attribute carrying data.frame(class, n_tss, median_total). Each gene
# contributes exactly one dominant TSS by construction, so a per-class TSS
# count would just restate n_genes and is reported as NA rather than as if it
# carried independent information; median_total (the dominant TSS's own
# pooled read depth) is reported, so it is comparable to the endtype family's
# depth-selection diagnostic even though the initiator family is not
# depth-selected by a read floor. Kept as an attribute rather than changing
# the return shape, so callers that only want the sets (including the unit
# tests) are unaffected.
build_initiator_sets <- function(a) {
    a <- a[!is.na(a$initiator) & a$initiator %in% DINUCS, ]
    # Dominant TSS per gene: highest mean_total, ties broken by cluster_id so
    # the choice is deterministic.
    ord <- order(a$gene, -a$mean_total, a$cluster_id)
    a   <- a[ord, ]
    dom <- a[!duplicated(a$gene), ]

    out   <- setNames(vector("list", length(DINUCS)), DINUCS)
    stats <- vector("list", length(DINUCS))
    names(stats) <- DINUCS
    has_total <- "n_total" %in% colnames(dom)
    for (dn in DINUCS) {
        sel <- dom$initiator == dn
        out[[dn]] <- unique(dom$gene[sel])
        stats[[dn]] <- data.frame(
            class        = dn,
            n_tss        = NA_integer_,
            median_total = if (has_total && sum(sel)) median(dom$n_total[sel]) else NA_real_,
            stringsAsFactors = FALSE)
        message(sprintf("  %s  %5d genes", dn, length(out[[dn]])))
    }
    attr(out, "stats") <- do.call(rbind, stats)
    out
}

main <- function() {
    opt <- parse_args(OptionParser(option_list = list(
        make_option("--outdir",   type = "character"),
        make_option("--min_set",  type = "integer", default = 10L),
        make_option("--pvalue",   type = "double",  default = 0.05),
        make_option("--qvalue",   type = "double",  default = 0.2),
        make_option("--simplify", type = "double",  default = 0.7),
        make_option("--orgdb",    type = "character", default = "org.Sc.sgd.db"),
        make_option("--keytype",  type = "character", default = "ORF"),
        make_option("--force",    action = "store_true", default = FALSE)
    )))

    res_out <- file.path(opt$outdir, "go_initiator.tsv")
    if (file.exists(res_out) && !opt$force) {
        message("[skip] outputs exist; use --force to regenerate")
        quit(save = "no")
    }

    .f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    source(file.path(dirname(sub("^--file=", "", .f[1])), "go_helpers.R"))

    a <- read.delim(file.path(opt$outdir, "tss_gene_assignment.tsv"),
                    stringsAsFactors = FALSE)
    universe <- readLines(file.path(opt$outdir, "universe.txt"))
    message(sprintf("Loaded %d assignments over %d genes in the universe",
                    nrow(a), length(universe)))

    message("Building foreground sets:")
    sets <- build_initiator_sets(a)

    message("Running enrichment:")
    res <- run_enrichment(sets, universe, family = "initiator",
                          orgdb = opt$orgdb, keytype = opt$keytype,
                          p_cut = opt$pvalue, q_cut = opt$qvalue,
                          simplify_cut = opt$simplify, min_set = opt$min_set,
                          set_stats = attr(sets, "stats"))

    write_enrichment(res, opt$outdir, "go_initiator")
    sessionInfo()
}

if (sys.nframe() == 0L) main()
