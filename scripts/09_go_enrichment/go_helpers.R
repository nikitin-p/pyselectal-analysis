# Shared GO enrichment machinery for Stage 16.
#
# Both families — end types and initiator dinucleotides — must be tested
# identically, so the statistics live here and neither enrichment script
# reimplements them. Source it with:
#
#   .f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
#   source(file.path(dirname(sub("^--file=", "", .f[1])), "go_helpers.R"))
#
# Benjamini-Hochberg correction is applied by enrichGO across GO terms within a
# single test. Families are corrected separately and never pooled, because they
# answer separate questions.

suppressMessages({
    library(clusterProfiler)
    library(AnnotationDbi)
})

# Turn enrichGO's "12/345" ratio strings into a numeric fold enrichment.
.ratio <- function(x) {
    parts <- strsplit(as.character(x), "/", fixed = TRUE)
    vapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
}

.empty_results <- function() {
    data.frame(family = character(), class = character(), ID = character(),
               Description = character(), GeneRatio = character(),
               BgRatio = character(), fold_enrichment = numeric(),
               pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(),
               Count = integer(), geneID = character(),
               stringsAsFactors = FALSE)
}

# Pull one class's value out of the caller-supplied set_stats table, or NA if
# set_stats is absent or has no row for that class.
.set_stat <- function(set_stats, cls, col) {
    if (is.null(set_stats) || !nrow(set_stats)) return(NA)
    v <- set_stats[[col]][set_stats$class == cls]
    if (length(v)) v[1] else NA
}

# sets      : named list of character vectors of gene identifiers
# universe  : character vector of gene identifiers
# family    : label written into the family column ("endtype" or "initiator")
# set_stats : optional data.frame(class, n_tss, median_total) describing how
#             depth-selected each foreground is. Not every family has a
#             natural n_tss (see the calling scripts); columns absent from
#             set_stats, or set_stats itself, are reported as NA rather than
#             fabricated. This is purely descriptive metadata carried into the
#             set-sizes table — it plays no part in the enrichment statistics.
# Returns list(results = data.frame, sizes = data.frame).
run_enrichment <- function(sets, universe, family, orgdb, keytype,
                           p_cut, q_cut, simplify_cut, min_set,
                           set_stats = NULL) {

    stopifnot(is.list(sets), !is.null(names(sets)), length(universe) > 0)

    suppressMessages(library(orgdb, character.only = TRUE))
    db    <- get(orgdb)
    valid <- AnnotationDbi::keys(db, keytype = keytype)

    universe_in <- unique(universe)
    universe    <- intersect(universe_in, valid)
    n_unmapped  <- length(universe_in) - length(universe)

    # A large unmapped fraction means the gene identifiers do not belong to this
    # annotation package, which would invalidate every result below, so it is
    # reported rather than silently dropped.
    message(sprintf("Universe: %d genes, %d mapped to %s (%d unmapped, %.1f %%)",
                    length(universe_in), length(universe), orgdb,
                    n_unmapped, 100 * n_unmapped / length(universe_in)))

    if (!length(universe)) {
        stop("No gene in the universe maps to ", orgdb, " under keyType = ",
             keytype, ". This indicates a reference or identifier mismatch.")
    }
    if (n_unmapped > 0.5 * length(universe_in)) {
        stop(sprintf("More than half the universe (%d of %d) does not map to %s. ",
                     n_unmapped, length(universe_in), orgdb),
             "This is a reference or identifier mismatch, not a filtering artefact.")
    }

    all_res  <- list()
    all_size <- list()

    for (cls in names(sets)) {
        genes  <- intersect(unique(sets[[cls]]), universe)
        n_all  <- length(unique(sets[[cls]]))
        status <- NA_character_
        n_term <- 0L

        if (length(genes) < min_set) {
            status <- "skipped_small"
        } else {
            eg <- enrichGO(gene          = genes,
                           universe      = universe,
                           OrgDb         = db,
                           keyType       = keytype,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = p_cut,
                           qvalueCutoff  = q_cut,
                           readable      = FALSE)

            if (is.null(eg) || nrow(as.data.frame(eg)) == 0) {
                status <- "no_terms"
            } else {
                eg <- clusterProfiler::simplify(eg, cutoff = simplify_cut,
                                                by = "p.adjust",
                                                select_fun = min)
                df <- as.data.frame(eg)
                if (nrow(df) == 0) {
                    status <- "no_terms"
                } else {
                    status <- "ok"
                    n_term <- nrow(df)
                    df$fold_enrichment <- .ratio(df$GeneRatio) / .ratio(df$BgRatio)
                    df$family <- family
                    df$class  <- cls
                    all_res[[cls]] <- df[, c("family", "class", "ID", "Description",
                                             "GeneRatio", "BgRatio",
                                             "fold_enrichment", "pvalue",
                                             "p.adjust", "qvalue", "Count",
                                             "geneID")]
                }
            }
        }

        message(sprintf("  %-10s %5d genes (%5d mapped) -> %s (%d terms)",
                        cls, n_all, length(genes), status, n_term))

        all_size[[cls]] <- data.frame(family = family, class = cls,
                                      n_genes = n_all, n_mapped = length(genes),
                                      n_tss = .set_stat(set_stats, cls, "n_tss"),
                                      median_total = .set_stat(set_stats, cls, "median_total"),
                                      status = status, n_terms = n_term,
                                      stringsAsFactors = FALSE)
    }

    results <- if (length(all_res)) do.call(rbind, all_res) else .empty_results()
    rownames(results) <- NULL

    list(results = results, sizes = do.call(rbind, all_size))
}

# Write both tables. A family that produced nothing still writes a header-only
# results file, so that "tested and found nothing" never looks like "not run".
write_enrichment <- function(res, outdir, prefix) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write.table(res$results, file.path(outdir, paste0(prefix, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(res$sizes, file.path(outdir, paste0(prefix, "_set_sizes.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("Wrote %s.tsv (%d rows) and %s_set_sizes.tsv (%d rows)",
                    prefix, nrow(res$results), prefix, nrow(res$sizes)))
}
