#!/usr/bin/env Rscript
# Stage 16 step 2: GO Biological Process enrichment by 5' end type.
#
# A TSS joins type T's foreground when its pooled pct_T reaches the
# --quantile quantile of pct_T across all retained TSSs and it carries at least
# --min_type_reads reads of that type. A gene joins the set when any of its
# assigned TSSs qualifies.
#
# A single fixed percentage cutoff was rejected: the global shares are so uneven
# (1Sg is about 74 % of reads, 1Sa well under 1 %) that any one threshold either
# admits nearly every TSS or none. The sets are deliberately not mutually
# exclusive, because the question is about preference, not partitioning.
#
# The 'other' category is excluded, since it collapses several rare types and an
# enrichment within it could not be attributed to one 5' end chemistry.
#
# Input : results/go/tss_gene_assignment.tsv, results/go/universe.txt
# Output: results/go/go_endtype.tsv, results/go/go_endtype_set_sizes.tsv

suppressMessages(library(optparse))

END_TYPES <- c("1Sg", "2Sgg", "3Sggg", "1Sa", "1St", "1Sc", "M")

# a              : the assignment table
# types          : end types to build sets for
# quantile_cut   : quantile of pct_T defining the foreground
# min_type_reads : minimum pooled reads of type T at a qualifying TSS
#
# Returns the named list of character vectors of gene names, as before, with a
# "stats" attribute carrying data.frame(class, n_tss, median_total). n_tss
# counts TSS clusters, not genes, since several clusters can point at the same
# gene; median_total is the median pooled n_total (depth across all types) of
# those same clusters. go_min_type_reads is an absolute read count applied to
# types whose global share spans two orders of magnitude, so for the rarer
# types reaching it requires unusually deep TSSs; n_tss/median_total exist so
# that this depth selection is visible in the artifact rather than only in
# prose (see the report's "How this was computed" / "What it shows" for
# sec-go). Kept as an attribute rather than changing the return shape, so
# callers that only want the sets (including the unit tests) are unaffected.
build_endtype_sets <- function(a, types, quantile_cut, min_type_reads) {
    out   <- setNames(vector("list", length(types)), types)
    stats <- vector("list", length(types))
    names(stats) <- types
    has_total <- "n_total" %in% colnames(a)
    for (ty in types) {
        pct_col <- paste0("pct_", ty)
        n_col   <- paste0("n_",   ty)
        if (!all(c(pct_col, n_col) %in% colnames(a))) {
            stop("Assignment table lacks columns for type ", ty)
        }
        pct <- a[[pct_col]]
        cut <- as.numeric(quantile(pct, probs = quantile_cut, na.rm = TRUE))
        # Guard the degenerate case where the quantile is 0 and would otherwise
        # admit every TSS that has none of this type at all.
        sel <- pct >= cut & pct > 0 & a[[n_col]] >= min_type_reads
        out[[ty]] <- unique(a$gene[sel])
        stats[[ty]] <- data.frame(
            class        = ty,
            n_tss        = sum(sel),
            median_total = if (has_total && sum(sel)) median(a$n_total[sel]) else NA_real_,
            stringsAsFactors = FALSE)
        message(sprintf("  %-6s cutoff pct >= %6.2f  ->  %5d TSS, %5d genes",
                        ty, cut, sum(sel), length(out[[ty]])))
    }
    attr(out, "stats") <- do.call(rbind, stats)
    out
}

main <- function() {
    opt <- parse_args(OptionParser(option_list = list(
        make_option("--outdir",         type = "character"),
        make_option("--quantile",       type = "double",  default = 0.90),
        make_option("--min_type_reads", type = "integer", default = 2L),
        make_option("--min_set",        type = "integer", default = 10L),
        make_option("--pvalue",         type = "double",  default = 0.05),
        make_option("--qvalue",         type = "double",  default = 0.2),
        make_option("--simplify",       type = "double",  default = 0.7),
        make_option("--orgdb",          type = "character", default = "org.Sc.sgd.db"),
        make_option("--keytype",        type = "character", default = "ORF"),
        make_option("--force",          action = "store_true", default = FALSE)
    )))

    res_out <- file.path(opt$outdir, "go_endtype.tsv")
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
    sets <- build_endtype_sets(a, END_TYPES, opt$quantile, opt$min_type_reads)

    message("Running enrichment:")
    res <- run_enrichment(sets, universe, family = "endtype",
                          orgdb = opt$orgdb, keytype = opt$keytype,
                          p_cut = opt$pvalue, q_cut = opt$qvalue,
                          simplify_cut = opt$simplify, min_set = opt$min_set,
                          set_stats = attr(sets, "stats"))

    write_enrichment(res, opt$outdir, "go_endtype")
    sessionInfo()
}

if (sys.nframe() == 0L) main()
