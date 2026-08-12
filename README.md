# pyselectal-analysis

This is a computational framework for analyzing CAGE and NET-CAGE datasets to characterize 5′-end diversity at transcription start sites (TSS): soft-clip patterns (1Sg / 2Sgg / 3Sggg / M / other), initiator dinucleotide usage, promoter sequence architecture, and their dynamics across conditions and subcellular fractions.

---

## Datasets

| Dataset | Conditions | Genome | Status |
| --- | --- | --- | --- |
| Yeast CAGE (Börlin et al. 2019, FEMS Yeast Res.) | anaerobic, ethanol, glucose, nitrogen × 2 reps | sacCer3 | in progress |
| LUHMES CAGE + NET-CAGE | 3 differentiation time points | hg38 | awaiting data |
| NET-CAGE paper cell lines | 2–3 lines, CAGE + NET-CAGE | hg38 | awaiting data |
| TLDR-seq | — | hg38 | awaiting data |

## Reference genomes

| Alias | Source | BSgenome package |
| --- | --- | --- |
| `sacCer3` | UCSC Golden Path | `BSgenome.Scerevisiae.UCSC.sacCer3` |
| `hg38` | GENCODE GRCh38 primary assembly | `BSgenome.Hsapiens.UCSC.hg38` |
| sacCer3 | UCSC S. cerevisiae reference | `data/sacCer3.fa` |

---

## Conda environments

| Env | Purpose |
| --- | --- |
| `star` | STAR aligner |
| `tools` | samtools and other CLI tools |
| `pysam` | pyselectal (BAM filtering and counting) |
| `r_cage` | R + Bioconductor (most analysis scripts) |
| `seqarchr` | scikit-learn backend used by seqArchR through reticulate (Stage 15) |

Two stages deliberately run against the **system** R library rather than the `r_cage`
environment. Stage 15 needs `seqArchR`, `optparse`, and `getopt`, and it reaches its
factorisation backend through `RETICULATE_PYTHON`, which the wrapper points at the
`seqarchr` conda environment (`conda_env_seqarchr` in `params.yaml`). Stage 16 needs
`clusterProfiler` and `org.Sc.sgd.db`. Stage 14 needs only Quarto and a plain R
installation.

---

## Quick start

All parameters are stored in `config/params.yaml`; sample metadata is stored in `config/samples.tsv`. There are no hardcoded paths outside those files.

```bash
# Chained pipeline: stages 1–9, then 17 and 18, then the report (stage 14)
bash run_analysis.sh

# Skip mapping if BAMs already exist in results/bam/
bash run_analysis.sh --skip-mapping

# Force re-run even if outputs exist
bash run_analysis.sh --skip-mapping --force
```

The `--force` flag can be added to any individual stage script as well.

The remaining stages are run on their own, because they are either optional, expensive,
or dependent on a manual choice of parameters:

```bash
bash scripts/05_tss_end_types/run_motif.sh            # stages 9a and 9b
bash scripts/05_tss_end_types/run_motif_by_class.sh   # stage 9d
bash scripts/06_cager/run.sh                          # stage 11
bash scripts/06_expression/run.sh                     # stage 12
bash scripts/08_seqarchr/run.sh                       # stage 15 (~25 min)
bash scripts/09_go_enrichment/run.sh                  # stage 16
```

**Re-run hazard.** The `annotation_bed_yeast` parameter is now filled in, so a forced
re-run of `scripts/05_tss_end_types/run.sh` would newly populate the `annotation`,
`dist_nearest`, and `nearest_gene` columns of `results/tss/tss_clustered.tsv`, which are
currently blank and which stages 14, 15, and 16 all read. Running `bash run_analysis.sh
--force` therefore changes those columns underneath stages 14–16 without re-running them,
and the three stages should be re-run afterwards.

---

## Pipeline stages

### Stage 1 — STAR mapping + MAPQ filter

**Status:** The analysis is complete for the yeast dataset.  
**Input:** Trimmed and UMI-processed FASTQ files are expected in `data/`.  
**Output:** The filtered alignments are written to `results/bam/{sample}_uniq.bam`.

STAR is run with `--alignEndsType Local` so that soft-clips are preserved; only uniquely mapped reads (MAPQ = 255) are retained via `samtools view -q 255`. Mapping is run sequentially per sample and is multi-threaded within each sample (`threads: 10` in `params.yaml`). For sacCer3, the index was built with `--genomeSAindexNbases 11`.

Key STAR flags:

```text
--alignEndsType Local
--outSAMmultNmax 1
--outSAMattributes NH HI AS NM
--limitBAMsortRAM 3000000000   # required for some samples
```

---

### Stage 2 — Subsampling to equal depth

**Script:** `scripts/02_subsample/01_subsample.sh`  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/bam/*_uniq.bam`  
**Output:** The subsampled BAM files are written to `results/bam_subsampled/{sample}_sub.bam`.

The script performs the following steps:

1. The read count of every library is obtained with `samtools view -c -F 256` (counting reads, not alignments);
2. The Q1, Q3, and IQR are computed across all libraries;
3. Any library with depth below the Tukey lower fence (Q1 − 1.5×IQR) is **excluded** — excluded libraries are logged and dropped from all downstream analysis;
4. The subsampling target is set as `floor(subsample_factor × min(retained depths))` (default factor 0.95);
5. Each retained library is subsampled with `samtools view -s SEED.FRAC` (seed 42).

For the yeast dataset, all 8 libraries were retained and the target was 2,190,922 reads.

---

### Stage 3 — 5′-end type counting

**Script:** `scripts/03_five_prime/01_count.sh`  
**Tool:** pyselectal v3.0 (`--count` mode)  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** The subsampled BAM files listed in the `bam` column of `config/samples.tsv`.  
**Output:** The per-sample count tables are written to `results/five_prime/{sample}_counts.tsv`.

pyselectal inspects the CIGAR string and the first soft-clipped base of each read and outputs `type\tcount` rows sorted by descending count. Types include `1Sg`, `2Sggg`, `2Sggc`, `36Mctcac`, etc.

Key parameters from `params.yaml`:

- `collapse_threshold: 0.1` — types accounting for less than 0.1 % of reads are merged into `other`;
- `mapped_prefix: 5` — the number of matched bases shown for M-type labels.

---

### Stage 4 — Per-sample histograms, replicate heatmap, chi-square

**Script:** `scripts/03_five_prime/02_plot.R`  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/five_prime/`, `config/samples.tsv`, `config/params.yaml`  
**Output:** The figures and tables are written to `results/figures/`: `03_per_sample_histogram.pdf`, `03_replicate_heatmap.pdf`, `03_chisq_pairwise.tsv`.

Raw pyselectal types are collapsed into four broad categories:

| Category | Rule |
| --- | --- |
| `1Sg` | exactly one soft-clipped G |
| `2Sgg` | two soft-clipped bases, first is G |
| `M` | no 5′ soft-clip (mapped start) |
| `other` | all remaining soft-clips |

The script produces the following outputs:

- **Histogram** (`03_per_sample_histogram.pdf`): the category frequency is shown per sample in a faceted bar chart;
- **Heatmap** (`03_replicate_heatmap.pdf`): a samples × categories matrix is clustered on both axes. Values are Z-scores of log10(count + 1), which stabilizes variance across categories with different count magnitudes. The distance metric is Euclidean (pheatmap default) and the linkage is controlled by `replicate_hclust_method` in `params.yaml` (default: `complete`). Consistent replicates are expected to cluster together;
- **Chi-square** (`03_chisq_pairwise.tsv`): reads are pooled per condition and a pairwise chi-square test is applied to the count matrix, with BH-corrected p-values.

---

### Stage 5 — Soft-clip nucleotide composition

**Script:** `scripts/03_five_prime/03_softclip_composition.R`  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/five_prime/`, `config/samples.tsv`, `config/params.yaml`  
**Output:** The figures and tables are written to `results/figures/`: `03_softclip_nuc_composition.pdf`, `03_softclip_nuc_by_sample.pdf`, `03_softclip_nuc_composition.tsv`.

The nucleotide composition of non-1Sg soft-clipped bases is compared against literature RT non-templated addition preferences (Zhu et al. 2001, Biotechniques). The reference frequencies can be overridden via `rt_nontemplated_prefs` in `params.yaml`.

---

### Stage 6 — Per-type BAM selection + initiator dinucleotide proportions

**Scripts:** `scripts/04_dinuc/00_select_bam.sh`, `scripts/04_dinuc/run.sh` (calls `01_dinuc_proportion.R`)  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/bam_subsampled/`  
**Output:** The type-filtered BAMs are written to `results/bam_selected/{sample}_1Sg.bam` and `{sample}_2Sgg.bam`; the dinucleotide tables are written to `results/dinuc_1Sg/` and `results/dinuc_2Sgg/`.

pyselectal `--select` is used to filter each BAM to reads matching the 1Sg or 2Sgg pattern. The script `01_dinuc_proportion.R` then extracts the 5′ CTSS position of each read, queries a 2-bp window centred on the CTSS (strand-aware, CTSS ±1 bp), and computes raw-count-weighted dinucleotide proportions.

The genome is specified via the `genome` column in `samples.tsv` as either a BSgenome alias (`sacCer3`, `hg38`) or a path to a FASTA file for custom assemblies.

Per-sample outputs are `{sample}_dinuc_proportion.tsv` and `{sample}_dinuc_proportion.pdf`; the combined table is `all_dinuc_proportions.tsv`.

---

### Stage 7 — YR vs YC enrichment across conditions

**Script:** `scripts/04_dinuc/02_yr_enrichment.R`  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/dinuc_1Sg/all_dinuc_proportions.tsv`, `results/dinuc_2Sgg/all_dinuc_proportions.tsv`  
**Output:** The figures and tables are written to `results/figures/`: `07_yr_*.pdf`, `07_yr_chisq_pairwise.tsv`.

Dinucleotides are classified into three categories:

- **YR** — pyrimidine at −1 and purine at +1: CA, CG, TA, TG;
- **YC** — pyrimidine at −1 and cytosine at +1: CC, TC;
- **other** — all remaining dinucleotides.

The script produces the following outputs:

- **`07_yr_per_sample.pdf`** — a stacked bar chart showing YR / YC / other proportions per sample, faceted by end type (1Sg / 2Sgg);
- **`07_yr_by_condition.pdf`** — the mean ± SD YR proportion by condition, grouped by end type;
- **`07_dinuc_heatmap_{1Sg,2Sgg}.pdf`** — a tile heatmap of the mean proportion per dinucleotide × condition;
- **`07_dinuc_zscore_heatmap_{1Sg,2Sgg}.pdf`** — a hierarchical clustering heatmap where rows are initiator dinucleotides and columns are samples. The values are Z-scores of log10(count + 1) within each sample (column-wise scaling stabilizes variance and removes inter-sample depth differences). The distance metric is set by `dinuc_dist_method` (default: `euclidean`) and the linkage by `dinuc_hclust_method` (default: `ward.D2`);
- **`07_yr_chisq_pairwise.tsv`** — pairwise chi-square results for YR vs non-YR across conditions, BH-corrected, computed separately for each end type.

---

### Stage 8 — TSS × end-type matrix + initiator annotation (scripted)

**Scripts:** `scripts/05_tss_end_types/run.sh`  
**Output:** The outputs are written to `results/tss/` (`typed_ctss.tsv.gz`, `tss_matrix.tsv`) and `results/figures/08_*.pdf`.

Reads are classified by 5′-end type directly from the CIGAR string and sequence in R (strand-aware; no additional BAM files are required). CTSS are aggregated and clustered with CAGEr `distclu` (maxDist 20 bp); the dominant TSS per cluster is annotated with the initiator dinucleotide via `promoters() + getSeq()` (always 5′→3′). The figures include a heatmap sorted by pct_1Sg, an initiator boxplot, and a stacked bar chart.

---

### Stage 9 — TSS hierarchical clustering + BED annotation (scripted)

**Scripts:** `scripts/05_tss_end_types/04_cluster_annotate.R`, `05_plot_clusters.R`  
**Output:** The outputs are written to `results/tss/` (`tss_clustered.tsv`, `tss_hclust.rds`) and `results/figures/09_*.pdf`.

TSS are clustered hierarchically by their mean end-type vector using Euclidean distance and Ward.D2 linkage (the default number of clusters is 6, controlled by `--k_clusters`). Strand-aware BED annotation assigns each TSS to one of three classes — sense_genic, antisense_genic, or intergenic — together with the distance to the nearest gene. The annotation BED path is set by `annotation_bed_yeast` in `params.yaml`.

---

### Stage 9a — Sequence context of purity-matched start sites

**Scripts:** `scripts/05_tss_end_types/run_motif.sh` (calls `06_motif_analysis.R`)  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/tss/tss_matrix.tsv`, `data/sacCer3.ensGene.gtf`  
**Output:** The sequences and tables are written to `results/motifs/` and the figures to `results/figures/11_*.pdf`.

Start sites are first restricted to those carrying a YR initiator and lying inside an
annotated promoter window, so the comparison is made within a homogeneous pool. Only
then are the sites ranked by end-class purity and the top decile of each class retained,
which gives three matched sets (`1Sg`-dominated, `M`-dominated, and other-soft-clip
dominated). The genomic sequence from 50 bp upstream to 10 bp downstream is extracted
strand-aware, and the three sets are compared as sequence logos, as per-position
nucleotide frequencies, and as pairwise frequency differences.

Because the YR filter is applied before extraction, the pyrimidine at −1 and the purine
at 0 are imposed rather than observed and carry no information; likewise, the analysis
cannot ask whether `M` positions are promoters, because every set is promoter-localised
by construction.

---

### Stage 9b — Motif comparison by occurrence-weighted end type

**Scripts:** `scripts/05_tss_end_types/run_motif.sh` (calls `07_motif_by_occurrence.R`)  
**Status:** The analysis is complete for the yeast dataset.  
**Output:** `results/motifs/motif_occurrence_summary.tsv` and `results/figures/12_*.pdf`.

The Stage 9a comparison is repeated across the finer end-type vocabulary (`1Sa`, `1Sc`,
`1St`, `2Sgg`, `3Sggg`, and `M`), each contrasted against `1Sg` as the reference class,
so that a difference specific to one clipped base can be separated from a difference
shared by all clipped classes.

---

### Stage 9c — The G ladder

**Script:** `scripts/05_tss_end_types/08_g_ladder.R`  
**Status:** The analysis is complete for the yeast dataset.  
**Output:** The tables are written to `results/g_ladder/` and the figures to `results/figures/13_*.pdf`.

If the reverse transcriptase added each extra G independently with a fixed probability
*p*, the classes would form a geometric ladder in which the `2Sgg`/`1Sg` ratio and the
`3Sggg`/`2Sgg` ratio were both equal to *p*. Every sufficiently deep start position is
therefore characterised by its counts of `1Sg`, `2Sgg`, and `3Sggg`, the step ratio is
computed at each rung, and a dispersion statistic quantifies how far the position-level
variation exceeds binomial sampling under a single global rate. The whole calculation is
repeated after removing the outlier positions identified by their extreme extra-G rate,
so that the conclusion does not rest on a handful of aberrant loci. The outputs also
break the rate down by initiator dinucleotide and by sample, and fit a GLM whose
coefficients are written to `glm_sequence_coefficients.tsv`.

---

### Stage 9d — Dinucleotide class as a function of distance from the start site

**Scripts:** `scripts/05_tss_end_types/09_dinuc_class_by_position.R`, `run_motif_by_class.sh` (calls `10_motif_by_dinuc_class.R`)  
**Status:** The analysis is complete for the yeast dataset.  
**Output:** The tables are written to `results/dinuc_position/` and the figures to `results/figures/14_*.pdf` and `15_*.pdf`.

The initiator analysis of Stage 6 looks at a single dinucleotide, the −1/+1 pair. This
analysis slides that window across the region around the start site and asks, at each
offset, how the dinucleotide ending at that offset distributes over the four
purine/pyrimidine classes YR, YY, RY, and RR. Coordinates run from −4 to −1 and from +1
to +4 with no zero position, so the dinucleotide labelled +1 is exactly the classical
initiator. The observed proportions are compared against a background computed from the
dinucleotide composition of the whole sacCer3 genome and reported as a log2
observed-to-expected ratio. Dinucleotides that are impossible for a given end class at a
given offset, such as those beginning with the G that a clipped G would mask, are
excluded from both the observation and the background so that the comparison stays
internally consistent.

The companion script `10_motif_by_dinuc_class.R` rebuilds the sequence logos of Stage 9a
separately within each initiator dinucleotide class, using an asymmetric −35 to +20
window, because the yeast core promoter extends further upstream of the start site than
downstream of it. `3Sggg` is excluded from the default type list, since its reads do not
survive a four-way split into readable logos.

---

### Stage 10 — Genomic localization by end type

**Scripts:** `scripts/06_localization/run.sh`, `01_localization.R`  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/bam_selected/`, `results/bam_subsampled/`, `data/sacCer3.ensGene.gtf`  
**Output:** The outputs are written to `results/localization/` and `results/figures/10_*.pdf`.

The script compares the genomic localization of different 5′ end types (1Sg, 2Sgg, 3Sggg, M) across annotation categories to assess whether multi-G variants represent genuine TSS signals.

**Two annotation modes are available:**

- **Fast mode** (`--fast`): Uses direct TxDb-based annotation with `findOverlaps()`. Low memory, suitable for local runs. Categories: Promoter > 5'UTR > 3'UTR > CDS > Intron > Intergenic.
- **Full mode** (default): Uses ChIPseeker for annotation and CAGEfightR for bidirectional enhancer detection. Requires >32 GB RAM. Categories: Promoter > Enhancer > 5'UTR > 3'UTR > Exon > Intron > Downstream > Intergenic.

The reference scripts in `cageseq2_bin/` demonstrate state-of-the-art CAGE analysis with CAGEr and CAGEfightR (including proper enhancer calling via bidirectional balance score).

**Key findings (yeast dataset):**

| End Type | Promoter | CDS | Intergenic |
| --- | --- | --- | --- |
| 1Sg | 97.2% | 0.8% | 2.0% |
| 2Sgg | 97.1% | 1.0% | 1.9% |
| 3Sggg | 97.6% | 0.8% | 1.6% |
| M | 69.9% | 18.6% | 11.3% |

2Sgg and 3Sggg localize almost identically to 1Sg (~97% at promoters), supporting that they are genuine TSS signals with additional non-templated G additions. M-type shows dramatically different distribution (only 70% promoter).

---

### Stage 11 — CAGEr promoter clusters + gene-type enrichment (scripted)

**Scripts:** `scripts/06_cager/run.sh`  
**Output:** `results/cager/`

The scripts perform a fresh CAGEr run: `filterLowExpCTSS` (TPM ≥ 1), `distclu` (maxDist 20, keepSingletonsAbove 5), IQW calculation, and gene-type enrichment analysis among coding and noncoding subclasses.

---

### Stage 12 — Expression dynamics (scripted)

**Scripts:** `scripts/06_expression/run.sh`  
**Output:** `results/expression/`

Expression dynamics of TSS clusters with ≥ `expr_min_pct_1Sg` % 1Sg reads are tracked across all conditions.

---

### Stage 13 (planned)

The distribution of reads in repetitive vs unique regions is characterized; TLDR-seq observed-vs-expected 5′ ends are compared.

---

### Stage 14 — HTML analysis report

**Scripts:** `scripts/07_report/run.sh`  
**Status:** The report is complete for the yeast dataset.  
**Input:** The summary tables written by stages 3–18, namely `results/five_prime/`, `results/dinuc_all/`, `results/dinuc_1Sg/`, `results/tss/`, `results/localization/`, `results/motifs/`, `results/g_ladder/`, `results/dinuc_position/`, `results/seqarchr/`, `results/go/`, `results/dominance/`, `results/a8/`, and `results/figures/*.tsv`.  
**Output:** The rendered document is written to `reports/cage_5prime_report.html`.

The source document is `reports/cage_5prime_report.qmd`. It renders to a single self-contained HTML file with an embedded table of contents, and it covers the whole yeast analysis in fifteen sections. Each section states in prose how the quantity was computed, presents the figure, and then gives the biological reading; the document closes with a ranked synthesis and an explicit statement of what the analysis cannot settle.

The report regenerates every figure at render time from the pipeline's own summary tables rather than embedding the PDF files produced by the individual stage scripts. This choice has three consequences worth noting:

- the rendered HTML is fully self-contained and can be shared as a single file;
- the report cannot drift away from the analysis that it describes, because it reads the same tables;
- no BAM file is opened, so the stage runs in well under a minute and can be re-run freely.

The stage requires Quarto on `PATH` and an R installation with `data.table`, `ggplot2`, `dplyr`, `tidyr`, `scales`, `patchwork`, `pheatmap`, `Biostrings`, `ggseqlogo`, `ggdendro`, and `knitr`. It does not require the `r_cage` conda environment.

```bash
bash scripts/07_report/run.sh           # skips if the HTML already exists
bash scripts/07_report/run.sh --force   # re-renders unconditionally
```

The report is rendered last, so any stage that is re-run should be followed by a forced
re-render.

---

### Stage 15 — De novo promoter sequence architectures (seqArchR)

**Scripts:** `scripts/08_seqarchr/run.sh` (calls `01_prepare_seqs.R`, `02_run_seqarchr.R`, `03_analyze_clusters.R`)  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/tss/tss_clustered.tsv` and the genome FASTA.  
**Output:** The tables and objects are written to `results/seqarchr/` and the figures to `results/figures/16_arch_*.pdf`.

This stage discovers promoter sequence architectures around the dominant TSSs without
supervision, and then asks whether those architectures are associated with the 5′ end
type. seqArchR clusters the windows by non-negative matrix factorisation on a one-hot
dinucleotide encoding and never sees the end-type labels, so the association it produces
is independent rather than circular. The window is asymmetric, running 150 bp upstream
and 50 bp downstream, because the yeast TATA box sits roughly 40–120 bp before the start
site.

**Collation is required rather than optional.** The final seqArchR iteration is
deliberately left uncollated, so every chunk in which the factorisation selected K = 1
survives as a separate cluster of exactly `seqarchr_chunk_size` sequences; the raw yeast
run produced 47 such clusters, of which 23 were chunks rather than architectures.
`02_run_seqarchr.R` therefore cuts the collation dendrogram from
`collate_seqArchR_result()` at `seqarchr_n_arch`. The cut point is chosen from where the
dendrogram's merge heights flatten, recorded in `arch_dendrogram_heights.tsv`, and never
from the end-type labels, because selecting the cut that maximised the association would
make the association test circular. Cramér's V stays in the 0.30–0.50 range for cuts
between 3 and 12, so no conclusion depends on the choice.

```bash
bash scripts/08_seqarchr/run.sh              # full run, roughly 25 minutes
bash scripts/08_seqarchr/run.sh --recollate  # reuse the saved factorisation, redo collation only
```

**Key findings (yeast dataset):** eight architectures are recovered, the chi-square
p-value is effectively zero, and Cramér's V is 0.483, with the `M` fraction ranging from
6 % to 89 % against a global 27 %. That `M` fraction correlates with the frequency of a
genomic G at the first templated base (Spearman 0.83), which independently recovers the
masking confound argued from the labels elsewhere in the report. At least one
architecture is strongly `M`-enriched without being correspondingly G-rich, and that
residual is the candidate for genuinely different biology.

---

### Stage 16 — GO Biological Process enrichment

**Scripts:** `scripts/09_go_enrichment/run.sh` (calls `01_assign_genes.R`, `02_enrich_endtype.R`, `03_enrich_initiator.R`, `04_plot.R`)  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/tss/tss_clustered.tsv`, `data/sacCer3.genes.bed`  
**Output:** The tables are written to `results/go/` and the figures to `results/figures/17_go_*.pdf`.

This stage tests whether the promoters favouring a given 5′ end type, or carrying a given
initiator dinucleotide, are enriched for particular GO Biological Process terms. Each TSS
cluster is assigned to the gene it plausibly drives, defined as the nearest annotated 5′
end downstream on the same strand within `go_upstream_window` bp upstream or
`go_downstream_window` bp into the gene. Both families are then tested with
`clusterProfiler::enrichGO` against a shared universe of genes that have an assigned and
expressed TSS, so a hit means enrichment relative to the detected transcriptome rather
than to the genome.

The end-type foregrounds use a per-type quantile (`go_quantile`) rather than one fixed
percentage, because the global shares of the types are too uneven for any single cutoff
to be fair to all of them; the initiator foregrounds assign each gene to the initiator of
its dominant TSS. Foreground sets smaller than `go_min_set_size` are skipped rather than
tested.

---

### Stage 17 — Single-type-dominated start sites

**Scripts:** `scripts/10_dominant_endtype/run.sh` (calls `01_build_fine_ctss.R`, `02_dominance.R`, `03_expression.R`, `04_context.R`)  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/bam_subsampled/`, `results/tss/tss_matrix.tsv`  
**Output:** The tables are written to `results/dominance/` and the figures to `results/figures/18_dom_*.pdf`.  
**Tests:** `bash test/run_dominance_test.sh`

This stage asks whether any individual start site predominantly emits one specific
unusual 5′ end type, and how productive such sites are. Every read in the subsampled BAMs
is re-classified by the literal sequence of its 5′ soft clip, which resolves the
collapsed `other` class of Stage 8 into an open-ended vocabulary (`1Sg`, `2Sga`, `3Sgtg`,
and so on); clips longer than `dom_max_clip` collapse to `longclip` and clips carrying an
ambiguous base collapse to `Nclip`. The fine-typed positions are re-clustered with the
Stage 8 distance, and the script **asserts** that the resulting cluster set is identical
to `results/tss/tss_matrix.tsv`, which is what licenses comparing this stage's classes
against the rest of the pipeline.

A cluster is called dominated only after four filters applied in sequence:

1. the pooled depth reaches `dom_min_depth`;
2. the top type's fraction reaches the dominance threshold;
3. that fraction is significantly above the type's genome-wide background under a
   one-sided binomial test with Benjamini-Hochberg correction;
4. the same type dominates in at least `dom_min_samples` of the samples that carry at
   least `dom_min_reads_per_sample` reads.

A **supplementary aggregate test** pools every type other than `1Sg` and `M` and applies
the same filters to the pooled share, because the per-type cascade takes an argmax over
the fine types and therefore cannot see a site whose unusual reads are split across
sibling clips.

**Key findings (yeast dataset):** 6,546 of 27,383 clusters are dominated at threshold
0.7, of which 5,675 are `1Sg` and 823 are `M`, and only 82 of the 10,455 clusters deep
enough to be classifiable are predominantly unusual even when the unusual types are
pooled. The templated-artefact control is structurally near-vacuous, because STAR
`--local` extends a templated 5′ base into the alignment instead of clipping it.

---

### Stage 18 — A−8 saturation of ⟨5′ end type, dinucleotide⟩ pairs

**Scripts:** `scripts/11_a8_saturation/run.sh` (calls `01_build_windows.R`, `02_a8_test.R`, `03_plot.R`)  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/dominance/fine_typed_ctss.tsv.gz`, `results/tss/tss_matrix.tsv`, the genome FASTA.  
**Output:** The tables are written to `results/a8/` and the figures to `results/figures/19_a8_*.pdf`.  
**Tests:** `bash test/a8/run_a8_test.sh`

This stage asks whether a given combination of 5′ end type and initiator dinucleotide is
saturated with genuine start sites, using adenine at position −8 as an operational marker
of real initiation. The grid holds 112 pairs, of which 24 are structurally empty and 88
are testable in principle; in the primary configuration 52 reach the floor of
`a8_min_positions` distinct positions and are tested. The observed A−8 frequency of each
pair is compared against that pair's own upstream promoter background by a two-sided
binomial test, so the comparison is self-normalising against pairs that simply sit in
more A-rich promoters.

Several choices are frozen before any pairwise test is run, and they matter for reading
the result:

- the background window is fixed at −20 … −11, chosen from the pooled per-position A
  profile, which is flat at 32.3–33.5 % there against 39.8 % at −8, while −9 and −7 are
  shoulders of the peak at 35.8 % and 35.6 % and −2 … −7 are all elevated;
- position −1 is excluded from every background set, because it **is** the first base of
  the initiator dinucleotide and is therefore fixed within a pair, so including it would
  bias the four `A`-initiated dinucleotides toward spurious depletion and the rest toward
  spurious enrichment;
- the test counts distinct CTSS rather than reads, because CAGE depth is skewed enough
  that read-weighting amounts to pseudo-replication, and a cluster bootstrap over distclu
  clusters puts a confidence interval on each effect size, because start sites within one
  cluster share most of their window sequence;
- positions below Stage 8's depth floor receive no cluster of their own and are attached
  to the nearest cluster within `a8_max_dist` rather than made into singleton blocks,
  because they were excluded by depth and not by isolation, and 39.5 % of them do lie
  within that distance of a clustered position;
- the 24 structurally empty pairs, which local alignment forbids, are **asserted** empty,
  which doubles as a check on strand handling;
- p-values are computed in log space and BH-corrected in log space, because the exact
  binomial underflows to zero at these depths.

Three sensitivity configurations each vary one fixed choice, and their agreement with the
primary configuration is reported in `19_a8_sensitivity.pdf`.

**Key findings (yeast dataset):** no tested pair is depleted for A−8, the effect sizes
running from +0.03 to +1.27 in log₂ units, so the marker never argues against a pair. 49
of the 52 tested pairs are significant and the same 49 are bootstrap-robust. Every
clipped class sits far above `M`, with median log₂ fold changes of 1.04 for `1Sc`, 0.89
for `3Sggg`, 0.70 for `1Sa`, 0.69 for `2Sgg`, 0.61 for `1Sg`, and 0.52 for `1St`, against
0.20 for `M`. The three pairs that fail significance are all `M` on a G- or A-initiated
dinucleotide, which is precisely where the G-masking argument of stages 15 and 17
predicts hidden capped ends and where this marker nonetheless finds essentially nothing.

---

## Repo layout

```text
config/
  samples.tsv           all samples (bam column points to subsampled BAMs)
  samples_1Sg.tsv       generated by 00_select_bam.sh — 1Sg-filtered BAMs
  samples_2Sgg.tsv       generated by 00_select_bam.sh — 2Sgg-filtered BAMs
  params.yaml           all pipeline parameters and conda env names
data/                   raw FASTQ (per-dataset subdirs or symlinks)
results/
  bam/                  STAR output, MAPQ-filtered
  bam_subsampled/       depth-equalised BAMs (*_sub.bam)
  bam_selected/         type-filtered BAMs (*_1Sg.bam, *_2Sgg.bam, *_3Sggg.bam)
  five_prime/           pyselectal count TSVs per sample
  dinuc/                dinucleotide proportions (all reads)
  dinuc_1Sg/            dinucleotide proportions (1Sg reads only)
  dinuc_2Sgg/           dinucleotide proportions (2Sgg reads only)
  dinuc_3Sggg/          dinucleotide proportions (3Sggg reads only)
  dinuc_all/            dinucleotide proportions merged across all samples
  dinuc_position/       dinucleotide class as a function of offset (Stage 9d)
  tss/                  typed_ctss.tsv.gz, tss_matrix.tsv, tss_clustered.tsv, tss_hclust.rds
  motifs/               purity-matched TSS sequences and PWM tables (Stages 9a, 9b, 9d)
  g_ladder/             G-ladder step ratios, dispersion, GLM coefficients (Stage 9c)
  localization/         genomic localization by end type (category_counts.tsv, category_summary.tsv)
  seqarchr/             promoter windows, factorisation, architectures (Stage 15)
  go/                   TSS-to-gene assignment, universe, enrichment tables (Stage 16)
  dominance/            fine-typed CTSS and dominance calls (Stage 17)
  a8/                   pair windows, A-8 tests, position profile (Stage 18)
  cager/                CAGEr promoter clusters
  cageseq/              CAGEr SE objects
  annotation/           ChIPseeker output
  expression/           (Stage 12) TSS cluster expression dynamics
  figures/              all PDFs and TSVs from analysis scripts
reports/                cage_5prime_report.qmd and the rendered HTML report
scripts/
  00_setup/             test BAM generation
  01_map/               STAR mapping wrapper
  02_subsample/         depth equalisation with low-coverage exclusion
  03_five_prime/        5′-end counting, plotting, soft-clip composition
  04_dinuc/             BAM selection, dinucleotide profiling, YR enrichment
  05_tss_end_types/     CIGAR-based read classification, distclu clustering, annotation,
                        motif analysis, G ladder, dinucleotide class by position
  06_localization/      genomic localization by end type (ChIPseeker/TxDb annotation)
  06_cager/             (scripted) CAGEr promoter clusters, gene-type enrichment
  06_expression/        (scripted) TSS cluster expression dynamics
  07_report/            Quarto HTML report rendering
  08_seqarchr/          de novo promoter sequence architectures
  09_go_enrichment/     GO Biological Process enrichment by end type and initiator
  10_dominant_endtype/  single-type-dominated start sites
  11_a8_saturation/     A-8 saturation of end type and dinucleotide pairs
  99_additional/        (planned) repeats, TLDR-seq
  utils/                log.sh, parse_yaml.sh
logs/                   per-script timestamped logs
envs/                   R session info, environment exports
test/                   synthetic BAMs and reference for pipeline testing
```

## Key parameters (`config/params.yaml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `collapse_threshold` | 0.1 | pyselectal: collapse types < 0.1 % into `other` |
| `mapped_prefix` | 5 | pyselectal: bases shown for M-type labels |
| `subsample_factor` | 0.95 | target = 0.95 × min(retained library depths) |
| `subsample_seed` | 42 | random seed for samtools subsampling |
| `replicate_hclust_method` | complete | linkage for Stage 4 replicate heatmap |
| `dinuc_hclust_method` | ward.D2 | linkage for Stage 7 dinucleotide z-score heatmap |
| `dinuc_dist_method` | euclidean | distance for Stage 7 dinucleotide clustering |
| `chisq_min_count` | 5 | minimum expected count per chi-square cell |
| `expr_min_pct_1Sg` | 70 | % 1Sg threshold selecting high-1Sg clusters in Stage 12 |
| `seqarchr_upstream` / `seqarchr_downstream` | 150 / 50 | asymmetric promoter window for Stage 15 |
| `seqarchr_chunk_size` | 500 | sequences per NMF chunk; an uncollated chunk has exactly this size |
| `seqarchr_n_arch` | 8 | collation dendrogram cut, chosen from merge heights only |
| `conda_env_seqarchr` | seqarchr | environment supplying the scikit-learn backend |
| `go_quantile` | 0.90 | per-type quantile of pct_T defining each end-type foreground |
| `go_upstream_window` / `go_downstream_window` | 500 / 100 | TSS-to-gene assignment window |
| `go_min_set_size` | 10 | foregrounds smaller than this are skipped rather than tested |
| `dom_min_depth` | 25 | minimum pooled reads for a cluster to be considered in Stage 17 |
| `dom_thresholds` | [0.7, 0.8] | dominance fractions, both evaluated |
| `dom_max_clip` | 4 | clips longer than this collapse to `longclip` |
| `dom_min_total` / `dom_max_dist` | 5 / 20 | must equal Stage 8's `--min_total` and `--max_dist` |
| `a8_obs_coord` | −8 | the position under test in Stage 18 |
| `a8_bg_plateau` | −20 … −11 | frozen background window for the binomial test |
| `a8_min_positions` | 200 | distinct qualifying CTSS needed for a pair to be tested |
| `a8_n_boot` | 1000 | cluster bootstrap resamples, primary configuration only |
| `a8_min_total` / `a8_max_dist` | 5 / 20 | must equal Stage 8's `--min_total` and `--max_dist` |

## Tests

```bash
bash test/run_dominance_test.sh   # stage 17
bash test/a8/run_a8_test.sh       # stage 18
```
