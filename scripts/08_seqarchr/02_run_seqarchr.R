#!/usr/bin/env Rscript
# Stage 08 step 2: de novo sequence-architecture discovery with seqArchR.
#
# Runs unsupervised NMF-based clustering on the promoter windows. End-type
# labels are deliberately NOT used here — they are only cross-tabulated
# afterwards, so the architecture clusters remain an independent partition.
#
# Input : results/seqarchr/promoter_windows.fa
# Output: results/seqarchr/seqarchr_result.rds
#         results/seqarchr/seq_cluster_assignment.tsv

suppressMessages({
    library(optparse)
    library(seqArchR)
    library(Biostrings)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--fasta",      type = "character"),
    make_option("--outdir",     type = "character"),
    make_option("--upstream",   type = "integer", default = 150L),
    make_option("--downstream", type = "integer", default = 50L),
    make_option("--total_itr",  type = "integer", default = 4L),
    make_option("--chunk_size", type = "integer", default = 500L),
    make_option("--cv_folds",   type = "integer", default = 5L),
    make_option("--n_cores",    type = "integer", default = 4L),
    make_option("--n_arch",     type = "integer", default = 8L),
    make_option("--seed",       type = "integer", default = 42L),
    make_option("--force",      action = "store_true", default = FALSE),
    make_option("--recollate",  action = "store_true", default = FALSE,
                help = "Reuse the saved factorisation and only redo collation")
)))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
rds_out <- file.path(opt$outdir, "seqarchr_result.rds")
tsv_out <- file.path(opt$outdir, "seq_cluster_assignment.tsv")

if (file.exists(rds_out) && file.exists(tsv_out) && !opt$force && !opt$recollate) {
    message("[skip] outputs exist; use --force to regenerate")
    quit(save = "no")
}

set.seed(opt$seed)

raw <- prepare_data_from_FASTA(opt$fasta, raw_seq = TRUE)
message(sprintf("Loaded %d sequences of width %d", length(raw), unique(width(raw))))

# The factorisation is the expensive part and is independent of the collation
# cut, so --recollate reuses it and only re-derives the architectures.
if (opt$recollate && file.exists(rds_out)) {
    message("[recollate] reusing saved factorisation: ", rds_out)
    res <- readRDS(rds_out)
} else {

ohe <- get_one_hot_encoded_seqs(raw, sinuc_or_dinuc = "dinuc")
message(sprintf("One-hot dinucleotide matrix: %d features x %d sequences",
                nrow(ohe), ncol(ohe)))

cfg <- set_config(
    parallelize = opt$n_cores > 1,
    n_cores     = opt$n_cores,
    chunk_size  = opt$chunk_size,
    cv_folds    = opt$cv_folds,
    bound       = 1e-8,
    flags       = list(debug = FALSE, time = TRUE, verbose = TRUE, plot = FALSE)
)

# Collate clusters on every iteration except the last, which is the seqArchR
# recommendation: intermediate collation merges redundant architectures while
# the final pass is left uncollated so genuinely distinct ones are not fused.
ocollation <- c(rep(TRUE, opt$total_itr - 1L), FALSE)

message(sprintf("Running seqArchR: %d iterations, chunk_size=%d, cores=%d",
                opt$total_itr, opt$chunk_size, opt$n_cores))

res <- seqArchR(
    config         = cfg,
    seqs_ohe_mat   = ohe,
    seqs_raw       = raw,
    seqs_pos       = seq(-opt$upstream, opt$downstream - 1L),
    total_itr      = opt$total_itr,
    set_ocollation = ocollation
)

saveRDS(res, rds_out)
}

# ---- flatten final-iteration cluster assignment ---------------------------
# The final seqArchR iteration is deliberately left uncollated, so chunks in
# which NMF selected K = 1 survive as separate clusters of exactly chunk_size.
# Those are an artefact of chunking rather than distinct architectures, so the
# raw clusters are collated here by hierarchical clustering of their basis
# vectors, which is seqArchR's own post-hoc merge step.
final_labels <- res$seqsClustLabels[[opt$total_itr]]
raw_list     <- seqArchR::get_seqs_clust_list(final_labels)
message(sprintf("Raw clusters before collation: %d (%d of them exactly chunk_size)",
                length(raw_list), sum(lengths(raw_list) == opt$chunk_size)))

# collate_seqArchR_result() returns a dendrogram over the raw clusters' basis
# vectors. It is cut at n_arch, which is chosen from the dendrogram's own height
# profile and never from the end-type labels: selecting the cut that maximises
# the label association would make the association test circular.
hc <- seqArchR::collate_seqArchR_result(
    result       = res,
    iter         = opt$total_itr,
    clust_method = "hc",
    aggl_method  = "ward.D",
    dist_method  = "euclid",
    return_order = TRUE,
    flags        = list(debugFlag = FALSE, verboseFlag = FALSE)
)
k_use <- min(opt$n_arch, length(raw_list))
grp   <- cutree(hc, k = k_use)
clust_list <- lapply(split(seq_along(raw_list), grp),
                     function(i) unlist(raw_list[i], use.names = FALSE))
message(sprintf("Collated %d raw clusters into %d architectures",
                length(raw_list), length(clust_list)))

saveRDS(hc, file.path(opt$outdir, "arch_dendrogram.rds"))
write.table(data.frame(k = seq_along(rev(hc$height)), merge_height = rev(hc$height)),
            file.path(opt$outdir, "arch_dendrogram_heights.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

assign_df <- do.call(rbind, lapply(seq_along(clust_list), function(i) {
    idx <- clust_list[[i]]
    data.frame(cluster_id = names(raw)[idx],
               arch_cluster = sprintf("A%02d", i),
               stringsAsFactors = FALSE)
}))
write.table(assign_df, tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("Found %d architecture clusters", length(clust_list)))
print(table(assign_df$arch_cluster))

sessionInfo()
