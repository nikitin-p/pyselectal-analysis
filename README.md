# pyselectal-analysis

A computational framework for analyzing CAGE and NET-CAGE datasets to study 5′-end diversity at transcription start sites (TSS), integrating read normalization, sequence bias profiling, initiator motif analysis, and clustering of TSS across cell types and conditions.

## Overview

This project characterizes TSS initiator usage and soft-clip patterns (1Sg / 2Sg / other) from CAGE, NET-CAGE, and TLDR-seq data, and tracks their dynamics across conditions.

### Datasets

- LUHMES CAGE + NET-CAGE across 3 differentiation time points
- CAGE + NET-CAGE from 2–3 cell lines of the original NET-CAGE paper (2 M urea)
- Yeast CAGE across multiple conditions (Lu & Lin, FEMS Yeast Res. 2019)
- TLDR-seq (comparison of observed vs expected 5′ ends)

### Genomes

hg38, sacCer3.

## Pipeline stages

1. STAR mapping (`--local`, SE) + MAPQ 255 filter for uniquely mapped reads
2. Subsampling to equal depth across retained BAMs (0.9 × min)
3. 5′-end extraction and classification (1Sg / 2Sg / other); per-sample histograms and replicate reproducibility
4. Nucleotide composition of non-1Sg soft-clipped bases
5. Genomic dinucleotide profiling around top 5′ ends
6. YR vs YC enrichment and statistical comparison across conditions
7. Per-initiator TSS × end-type distributions and visualization
8. Hierarchical clustering of TSS + ChIPseeker annotation
9. CAGEr promoter clusters and gene-type enrichment
10. Expression dynamics of interesting TSS clusters across conditions
11. Repetitive vs unique regions; TLDR-seq observed-vs-expected

## Repo layout

```text
config/      sample sheets, parameter files
data/        raw FASTQ (per-dataset subdirs or symlinks)
results/     bam/, bam_subsampled/, five_prime/, dinuc/, tss/,
             cager/, annotation/, figures/
scripts/     numbered per stage (01_map, 02_subsample, ...)
logs/        per-script timestamped logs
envs/        R session info, environment exports
```
