#!/usr/bin/env Rscript
# Stage 16 step 4: figures for the two enrichment families.
#
# Two views per family. The dotplot shows each class's own strongest terms,
# which answers "what is this class about". The matrix shows the union of those
# terms across all classes, which answers the harder question of whether a term
# is specific to one class or shared by all of them; a term that lights up
# everywhere is a property of the expressed transcriptome, not of the class.
#
# Input : results/go/go_endtype.tsv, results/go/go_initiator.tsv
# Output: results/figures/17_go_*.pdf

suppressMessages({
    library(optparse)
    library(ggplot2)
})

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

CANON     <- "#e63946"                    # canonical (YR) initiator colour, reserved project-wide
HIGHLIGHT_COLOR <- "#2a9d8f"               # AA/AG highlight; distinct from CANON, not a claim of canonicity
YR        <- c("CA", "CG", "TA", "TG")
HIGHLIGHT <- c("AA", "AG")                # the initiators question 2 asks about

opt <- parse_args(OptionParser(option_list = list(
    make_option("--outdir",  type = "character"),
    make_option("--figdir",  type = "character"),
    make_option("--top_n",   type = "integer", default = 8L),
    make_option("--force",   action = "store_true", default = FALSE)
)))

dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

out_files <- c(
    file.path(opt$figdir, "17_go_endtype_dotplot.pdf"),
    file.path(opt$figdir, "17_go_endtype_matrix.pdf"),
    file.path(opt$figdir, "17_go_initiator_dotplot.pdf"),
    file.path(opt$figdir, "17_go_initiator_matrix.pdf"))
if (all(file.exists(out_files)) && !opt$force) {
    message("[skip] outputs exist; use --force to regenerate")
    quit(save = "no")
}

wrap_terms <- function(x, width = 45) {
    vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
           character(1), USE.NAMES = FALSE)
}

# Empty results still produce a figure carrying the reason, so that a missing
# PDF always means "the stage did not run" and never "the stage found nothing".
empty_panel <- function(path, msg) {
    p <- ggplot() +
        annotate("text", x = 0, y = 0, label = msg, size = 5) +
        theme_void()
    ggsave(path, p, width = 8, height = 3)
    message("Wrote ", path, " (empty: ", msg, ")")
}

plot_family <- function(family, res_path, dot_path, mat_path,
                        class_levels, yr_classes = character(0),
                        highlight_classes = character(0)) {

    res <- read.delim(res_path, stringsAsFactors = FALSE)
    if (!nrow(res)) {
        empty_panel(dot_path, paste0("No GO BP term passed the cutoffs for any ",
                                     family, " class."))
        empty_panel(mat_path, paste0("No GO BP term passed the cutoffs for any ",
                                     family, " class."))
        return(invisible(NULL))
    }

    present <- class_levels[class_levels %in% res$class]
    res$class <- factor(res$class, levels = present)

    lab_cols <- ifelse(present %in% yr_classes, CANON,
                 ifelse(present %in% highlight_classes, HIGHLIGHT_COLOR, "black"))

    # ---- dotplot: each class's own top terms -----------------------------
    res <- res[order(res$class, res$p.adjust), ]
    top <- do.call(rbind, lapply(split(res, res$class), function(d) {
        head(d, opt$top_n)
    }))
    top$term <- wrap_terms(top$Description)
    term_order <- unique(top$term[order(top$class, -top$fold_enrichment)])
    top$term <- factor(top$term, levels = rev(term_order))

    p <- ggplot(top, aes(x = class, y = term)) +
        geom_point(aes(size = fold_enrichment, colour = -log10(p.adjust))) +
        scale_colour_viridis_c(option = "C", name = expression(-log[10]~FDR)) +
        scale_size_continuous(name = "fold enrichment", range = c(1.5, 7)) +
        labs(x = NULL, y = NULL,
             title = paste0("GO Biological Process by ", family),
             subtitle = paste0("top ", opt$top_n,
                               " terms per class after redundancy collapse")) +
        theme_bw(base_size = 9) +
        theme(axis.text.x = element_text(colour = lab_cols, face = "bold"),
              panel.grid.minor = element_blank())

    h <- max(4, 0.28 * length(levels(top$term)) + 2)
    ggsave(dot_path, p, width = 9, height = h, limitsize = FALSE)
    message("Wrote ", dot_path)

    # ---- matrix: union of those terms across every class -----------------
    keep <- unique(top$Description)
    m <- res[res$Description %in% keep, ]
    grid <- expand.grid(Description = keep, class = present,
                        stringsAsFactors = FALSE)
    grid <- merge(grid, m[, c("Description", "class", "p.adjust",
                              "fold_enrichment")],
                  by = c("Description", "class"), all.x = TRUE)
    grid$term <- wrap_terms(grid$Description)
    grid$term <- factor(grid$term, levels = rev(unique(wrap_terms(keep))))
    grid$class <- factor(grid$class, levels = present)
    grid$score <- -log10(grid$p.adjust)

    q <- ggplot(grid, aes(x = class, y = term, fill = score)) +
        geom_tile(colour = "grey90") +
        scale_fill_viridis_c(option = "C", na.value = "grey96",
                             name = expression(-log[10]~FDR)) +
        labs(x = NULL, y = NULL,
             title = paste0("Term specificity across ", family, " classes"),
             subtitle = "grey means the term did not pass the cutoffs for that class") +
        theme_bw(base_size = 9) +
        theme(axis.text.x = element_text(colour = lab_cols, face = "bold"),
              panel.grid = element_blank())

    ggsave(mat_path, q, width = 9, height = h, limitsize = FALSE)
    message("Wrote ", mat_path)
}

DINUCS <- as.vector(t(outer(c("A", "C", "G", "T"), c("A", "C", "G", "T"), paste0)))

plot_family("end type",
            file.path(opt$outdir, "go_endtype.tsv"),
            file.path(opt$figdir, "17_go_endtype_dotplot.pdf"),
            file.path(opt$figdir, "17_go_endtype_matrix.pdf"),
            END_TYPES_REAL, character(0), character(0))

plot_family("initiator",
            file.path(opt$outdir, "go_initiator.tsv"),
            file.path(opt$figdir, "17_go_initiator_dotplot.pdf"),
            file.path(opt$figdir, "17_go_initiator_matrix.pdf"),
            DINUCS, YR, HIGHLIGHT)

sessionInfo()
