#!/usr/bin/env Rscript
# Stage 22 step 9: do long clips start with an untemplated (cap-derived) G?
#
# Reads the case-encoded long-clip table of step 4 (longclip_case.tsv.gz, ms:Z
# with UPPER = match, lower = mismatch to the genome laid directly upstream of
# the aligned part) and locates every aligned 5' end as
#   promoter   within 100 nt upstream of an annotated gene start (gene BED),
#   intron     inside an annotated intron (nuclear and chrM separately),
#   elsewhere  in neither; positions in both are dropped.
#   A  % of reads whose 5'-most base is a G that mismatches the genome, for
#      unplaced and placed (spliced / templated) reads per region, with a
#      Fisher test of promoter against nuclear intron among unplaced reads;
#   B  among unplaced reads, % whose first clip base mismatches the genome,
#      split by whether that base is G, against the 75 % chance level.
#
# Outputs: <outdir>/lc_capG_untemplated_G.tsv (A), lc_capG_first_base.tsv (B);
#          <figdir>/29_longclip_capG.{pdf,png}
#
# Usage (from project root):
#   Rscript scripts/14_longclip/09_cap_g.R \
#       [--case results/longclip/bam/longclip_case.tsv.gz] \
#       [--genes data/sacCer3.genes.bed] [--gtf data/sacCer3.ensGene.gtf] \
#       [--outdir results/longclip] [--figdir results/figures] [--force]
suppressPackageStartupMessages({library(data.table);library(GenomicRanges);library(GenomicFeatures)
  library(ggplot2);library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
case_f  <- get_arg("--case",   "results/longclip/bam/longclip_case.tsv.gz")
genes_f <- get_arg("--genes",  "data/sacCer3.genes.bed")
gtf_f   <- get_arg("--gtf",    "data/sacCer3.ensGene.gtf")
outdir  <- get_arg("--outdir", "results/longclip")
figdir  <- get_arg("--figdir", "results/figures")
force   <- "--force" %in% args
fig_f <- file.path(figdir, "29_longclip_capG.pdf")
if (file.exists(fig_f) && !force) { message("[skip] ", fig_f, " exists (use --force)"); quit(save = "no") }
dir.create(outdir, showWarnings = FALSE, recursive = TRUE); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "utils", "txdb.R"))
cs <- fread(case_f)
g <- fread(genes_f, col.names=c("chr","s0","e","gene","sc","strand"))
prom <- with(g, GRanges(chr, IRanges(ifelse(strand=="+", pmax(1,s0+1-100), e+1),
                                     ifelse(strand=="+", s0, e+100)), strand=strand))
txdb <- suppressWarnings(make_txdb_from_gff(gtf_f, format="gtf"))
intr <- unique(unlist(intronsByTranscript(txdb)))
q <- GRanges(cs$chr, IRanges(cs$pos5, width=1), strand=cs$strand)
inP <- overlapsAny(q, prom); inI <- overlapsAny(q, intr)
REG <- c("promoter\n(100 nt upstream)", "intron\n(nuclear)", "intron\n(chrM)", "elsewhere")
cs[, region := fcase(inI & !inP & chr!="chrM", REG[2], inI & !inP & chr=="chrM", REG[3],
                     inP & !inI, REG[1], !inP & !inI, REG[4], default=NA_character_)]
cs <- cs[!is.na(region)][, region := factor(region, REG)]
cs[, b1 := substr(ms,1,1)][, `:=`(G1 = toupper(b1)=="G", mmG = b1=="g", mm1 = b1 %in% c("a","c","g","t"))]
cs[, grp := factor(ifelse(class=="unplaced","unplaced","placed (spliced / templated)"),
                   c("unplaced","placed (spliced / templated)"))]
wilson <- function(x, n, z=1.96) { p <- x/n; d <- 1+z^2/n
  c0 <- (p+z^2/(2*n))/d; h <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/d; list(lo=100*(c0-h), hi=100*(c0+h)) }
a <- cs[, .(n=.N, x=sum(mmG)), by=.(grp, region)][, pct := 100*x/n][, c("lo","hi") := wilson(x,n)]
u <- cs[class=="unplaced" & region %in% REG[1:2]]
ft <- fisher.test(table(u$region, u$mmG)[REG[1:2], ])
th <- theme_classic(base_size=9) + theme(plot.title=element_text(face="bold", size=10),
  plot.subtitle=element_text(size=8, colour="grey30"), legend.position="top", legend.key.size=unit(0.35,"cm"))
pd <- position_dodge(0.8)
pA <- ggplot(a, aes(region, pct, fill=grp)) +
  geom_col(position=pd, width=0.75, colour="white", linewidth=0.4) +
  geom_errorbar(aes(ymin=lo, ymax=hi), position=pd, width=0.2, linewidth=0.35) +
  geom_text(aes(y=hi+1.5, label=sprintf("%.0f%%\nn = %s", pct, format(n, big.mark=",", trim=TRUE))), position=pd, vjust=0,
            size=2.3, lineheight=0.85, colour="grey20") +
  scale_fill_manual(values=c("unplaced"="#dc2626","placed (spliced / templated)"="#1d4ed8"), name=NULL) +
  scale_y_continuous(limits=c(0,112), breaks=seq(0,100,25), expand=c(0,0)) +
  labs(x="region of the aligned 5' end", y="% reads whose 5'-most base is a\nG that mismatches the genome",
       title="A  Long clips starting with an untemplated G",
       subtitle=sprintf("unplaced, promoter vs nuclear intron: Fisher OR = %.2f (95%% CI %.2f-%.2f), p = %.2g; bars: 95%% Wilson CI",
                        ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value)) + th
b <- cs[class=="unplaced", .(n=.N, x=sum(mm1)), by=.(region, first=ifelse(G1,"G","A / C / T"))]
b[, pct := 100*x/n][, c("lo","hi") := wilson(x,n)][, first := factor(first, c("G","A / C / T"))]
pB <- ggplot(b, aes(region, pct, colour=first)) +
  geom_hline(yintercept=75, linetype=2, linewidth=0.35, colour="grey50") +
  geom_errorbar(aes(ymin=lo, ymax=hi), position=position_dodge(0.4), width=0.15, linewidth=0.4) +
  geom_point(position=position_dodge(0.4), size=2.4) +
  annotate("text", x=4.45, y=73, label="chance (75%)", size=2.3, hjust=1, colour="grey40") +
  scale_colour_manual(values=c("G"="#d97706","A / C / T"="#7c3aed"), name="first clip base") +
  scale_y_continuous(limits=c(40,100)) +
  labs(x="region of the aligned 5' end", y="% whose first base mismatches\nthe adjacent upstream genome",
       title="B  Unplaced clips: a G-initial clip is more often untemplated",
       subtitle="the clip is laid directly upstream of the aligned part; error bars: 95% Wilson CI") + th
ggsave(fig_f, pA / pB, width=7.5, height=8)
ggsave(sub("\\.pdf$", ".png", fig_f), pA / pB, width=7.5, height=8, dpi=150)
flat <- function(d) copy(d)[, region := gsub("\n", " ", region)]
fwrite(flat(a), file.path(outdir, "lc_capG_untemplated_G.tsv"), sep="\t")
fwrite(flat(b), file.path(outdir, "lc_capG_first_base.tsv"), sep="\t")
message(sprintf("unplaced, promoter vs nuclear intron, untemplated G: Fisher OR = %.2f (%.2f-%.2f), p = %.2g",
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
print(a); print(b)
message("Written: ", fig_f)
sessionInfo()
