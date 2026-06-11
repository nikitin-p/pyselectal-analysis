# pyselectal-analysis

This is a computational framework for analyzing CAGE and NET-CAGE datasets to characterize 5′-end diversity at transcription start sites (TSS): soft-clip patterns (1Sg / 2Sg / M / other), initiator dinucleotide usage, and their dynamics across conditions and subcellular fractions.

---

## Datasets

| Dataset | Conditions | Genome | Status |
| --- | --- | --- | --- |
| Yeast CAGE (Börlin et al. 2019, FEMS Yeast Res.) | anaerobic, ethanol, glucose, nitrogen × 2 reps | CENPK113-7D (custom) | in progress |
| LUHMES CAGE + NET-CAGE | 3 differentiation time points | hg38 | awaiting data |
| NET-CAGE paper cell lines | 2–3 lines, CAGE + NET-CAGE | hg38 | awaiting data |
| TLDR-seq | — | hg38 | awaiting data |

## Reference genomes

| Alias | Source | BSgenome package |
| --- | --- | --- |
| `sacCer3` | UCSC Golden Path | `BSgenome.Scerevisiae.UCSC.sacCer3` |
| `hg38` | GENCODE GRCh38 primary assembly | `BSgenome.Hsapiens.UCSC.hg38` |
| Custom yeast | CENPK113-7D assembly, Börlin et al. 2019 | FASTA path in `samples.tsv` |

---

## Conda environments

| Env | Purpose |
| --- | --- |
| `star` | STAR aligner |
| `tools` | samtools and other CLI tools |
| `pysam` | pyselectal (BAM filtering and counting) |
| `r_cage` | R + Bioconductor (all analysis scripts) |

---

## Quick start

All parameters are stored in `config/params.yaml`; sample metadata is stored in `config/samples.tsv`. There are no hardcoded paths outside those files.

```bash
# Run from project root
bash scripts/02_subsample/01_subsample.sh
bash scripts/03_five_prime/01_count.sh
conda run -n r_cage Rscript scripts/03_five_prime/02_plot.R
conda run -n r_cage Rscript scripts/03_five_prime/03_softclip_composition.R
bash scripts/04_dinuc/run.sh
```

The flag `--force` can be added to any bash script to rerun it even if the outputs already exist.

---

## Pipeline stages

### Stage 1 — STAR mapping + MAPQ filter

**Status:** The analysis is complete for the yeast dataset.  
**Input:** Trimmed and UMI-processed FASTQ files are expected in `data/`.  
**Output:** The filtered alignments are written to `results/bam/{sample}_uniq.bam`.

STAR is run with `--alignEndsType Local` so that soft-clips are preserved; only uniquely mapped reads (MAPQ = 255) are retained via `samtools view -q 255`. Mapping is run sequentially per sample and is multi-threaded within each sample (`threads: 10` in `params.yaml`). For the CENPK113-7D custom assembly, the index was built with `--genomeSAindexNbases 11`.

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

pyselectal inspects the CIGAR string and the first soft-clipped base of each read and outputs `type\tcount` rows sorted by descending count. Types include `1Sg`, `2Sgg`, `2Sgc`, `36Mctcac`, etc.

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
| `2Sg` | two soft-clipped bases, first is G |
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
**Output:** The type-filtered BAMs are written to `results/bam_selected/{sample}_1Sg.bam` and `{sample}_2Sg.bam`; the dinucleotide tables are written to `results/dinuc_1Sg/` and `results/dinuc_2Sg/`.

pyselectal `--select` is used to filter each BAM to reads matching the 1Sg or 2Sg pattern. The script `01_dinuc_proportion.R` then extracts the 5′ CTSS position of each read, queries a 2-bp window centred on the CTSS (strand-aware, CTSS ±1 bp), and computes raw-count-weighted dinucleotide proportions.

The genome is specified via the `genome` column in `samples.tsv` as either a BSgenome alias (`sacCer3`, `hg38`) or a path to a FASTA file for custom assemblies.

Per-sample outputs are `{sample}_dinuc_proportion.tsv` and `{sample}_dinuc_proportion.pdf`; the combined table is `all_dinuc_proportions.tsv`.

---

### Stage 7 — YR vs YC enrichment across conditions

**Script:** `scripts/04_dinuc/02_yr_enrichment.R`  
**Status:** The analysis is complete for the yeast dataset.  
**Input:** `results/dinuc_1Sg/all_dinuc_proportions.tsv`, `results/dinuc_2Sg/all_dinuc_proportions.tsv`  
**Output:** The figures and tables are written to `results/figures/`: `07_yr_*.pdf`, `07_yr_chisq_pairwise.tsv`.

Dinucleotides are classified into three categories:

- **YR** — pyrimidine at −1 and purine at +1: CA, CG, TA, TG;
- **YC** — pyrimidine at −1 and cytosine at +1: CC, TC;
- **other** — all remaining dinucleotides.

The script produces the following outputs:

- **`07_yr_per_sample.pdf`** — a stacked bar chart showing YR / YC / other proportions per sample, faceted by end type (1Sg / 2Sg);
- **`07_yr_by_condition.pdf`** — the mean ± SD YR proportion by condition, grouped by end type;
- **`07_dinuc_heatmap_{1Sg,2Sg}.pdf`** — a tile heatmap of the mean proportion per dinucleotide × condition;
- **`07_dinuc_zscore_heatmap_{1Sg,2Sg}.pdf`** — a hierarchical clustering heatmap where rows are initiator dinucleotides and columns are samples. The values are Z-scores of log10(count + 1) within each sample (column-wise scaling stabilizes variance and removes inter-sample depth differences). The distance metric is set by `dinuc_dist_method` (default: `euclidean`) and the linkage by `dinuc_hclust_method` (default: `ward.D2`);
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

### Stage 10 — CAGEr promoter clusters + gene-type enrichment (scripted)

**Scripts:** `scripts/06_cager/run.sh`  
**Output:** `results/cager/`

The scripts perform a fresh CAGEr run: `filterLowExpCTSS` (TPM ≥ 1), `distclu` (maxDist 20, keepSingletonsAbove 5), IQW calculation, and gene-type enrichment analysis among coding and noncoding subclasses.

---

### Stage 11 — Expression dynamics (scripted)

**Scripts:** `scripts/06_expression/run.sh`  
**Output:** `results/expression/`

Expression dynamics of TSS clusters with ≥ `expr_min_pct_1Sg` % 1Sg reads are tracked across all conditions.

---

### Stage 12 (planned)

The distribution of reads in repetitive vs unique regions is characterized; TLDR-seq observed-vs-expected 5′ ends are compared.

---

## Repo layout

```text
config/
  samples.tsv           all samples (bam column points to subsampled BAMs)
  samples_1Sg.tsv       generated by 00_select_bam.sh — 1Sg-filtered BAMs
  samples_2Sg.tsv       generated by 00_select_bam.sh — 2Sg-filtered BAMs
  params.yaml           all pipeline parameters and conda env names
data/                   raw FASTQ (per-dataset subdirs or symlinks)
results/
  bam/                  STAR output, MAPQ-filtered
  bam_subsampled/       depth-equalised BAMs (*_sub.bam)
  bam_selected/         type-filtered BAMs (*_1Sg.bam, *_2Sg.bam)
  five_prime/           pyselectal count TSVs per sample
  dinuc/                dinucleotide proportions (all reads)
  dinuc_1Sg/            dinucleotide proportions (1Sg reads only)
  dinuc_2Sg/            dinucleotide proportions (2Sg reads only)
  dinuc_all/            dinucleotide proportions merged across all samples
  tss/                  typed_ctss.tsv.gz, tss_matrix.tsv, tss_clustered.tsv, tss_hclust.rds
  cager/                CAGEr promoter clusters
  cageseq/              CAGEr SE objects
  annotation/           ChIPseeker output
  expression/           (Stage 11) TSS cluster expression dynamics
  figures/              all PDFs and TSVs from analysis scripts
scripts/
  00_setup/             test BAM generation
  01_map/               STAR mapping wrapper
  02_subsample/         depth equalisation with low-coverage exclusion
  03_five_prime/        5′-end counting, plotting, soft-clip composition
  04_dinuc/             BAM selection, dinucleotide profiling, YR enrichment
  05_tss_end_types/     CIGAR-based read classification, distclu clustering, annotation, figures
  06_cager/             (scripted) CAGEr promoter clusters, gene-type enrichment
  06_expression/        (planned) TSS cluster expression dynamics
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
