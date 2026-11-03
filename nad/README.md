# Research branch: are some 5′ end types NAD caps?

## Question

This branch asks whether any of the non-canonical 5′ end types seen in the yeast
CAGE data, and above all `1Sa`, mark RNAs that carry a nicotinamide adenine
dinucleotide (NAD) cap instead of an m7G cap. The reference dataset is the NAD
captureSeq study of Zhang et al. 2020 (Nat. Commun. 11:5508,
doi 10.1038/s41467-020-19326-3), which is deposited as GEO GSE146368, BioProject
PRJNA610259 and SRA SRP251550.

## Why the question is reasonable

RNA polymerase uses NAD as a non-canonical initiating nucleotide in place of ATP,
so every NAD-RNA begins with the adenosine of NAD at genomic position +1. Beyond
that adenosine sits a pyrophosphate and a nicotinamide riboside, which carries no
Watson–Crick base. That riboside does have a free 2′,3′ cis-diol, so the periodate
oxidation used in CAGE cap trapping can in principle biotinylate an NAD cap just
as it biotinylates the ribose of m7G. A CAGE read from an NAD-RNA should therefore
start at a genomic A, and whatever reverse transcriptase does opposite the
riboside decides how that read looks:

- if reverse transcriptase adds a non-templated T opposite the riboside, the read carries a soft-clipped A and is typed `1Sa`;
- if reverse transcriptase stops at the adenosine, the read starts on a mapped A and is typed `M`;
- if reverse transcriptase adds a C, as it does opposite m7G, the read is typed `1Sg` and is indistinguishable from an ordinary capped read that starts at an A.

Only the first outcome would make an end type an NAD marker, so the analysis has
to test the end-type composition at NAD 5′ ends rather than assume `1Sa`. Other
adenosine-containing metabolite caps, such as FAD, would produce the same read
structure, so a positive result means "a non-m7G adenosine-derived cap" rather
than NAD specifically.

## What the paper established

The findings that bear on this branch are listed below:

- NAD-RNAs in budding yeast are abundant (about 1,460 enriched species in wild type and 3,000–4,000 in the decapping mutants), short (mostly 20–170 nt), mostly 3′-truncated, and map mainly to mRNA 5′ ends;
- the NAD-modified fraction of any nuclear transcript is low, at 1–5 % by qRT-PCR;
- NAD is incorporated by RNA polymerase II at transcription initiation, often from a TSS other than the canonical m7G TSS of the same gene, and efficient incorporation needs a `YAAG` core motif in which the first A is +1 and the G at +3 matters most;
- three decapping enzymes remove NAD caps: Rai1 (nuclear and acting first), Dxo1 and Npy1 (cytosolic), so the `dxo1Δ rai1Δ npy1Δ` triple mutant shows the least distorted map of where NAD is incorporated, and the authors derived the motif from it;
- mitochondrial transcripts were not analysed in the paper, but the authors cite Bird et al. 2018, which reported that up to 60 % of some yeast mitochondrial mRNAs carry NAD, with mtRNAP preferring an `RA` start.

The low nuclear modification ratio sets the expected effect size. At an NAD-favoured
nuclear TSS no more than a few percent of CAGE reads can come from NAD-RNAs, which
is the same order as the overall `1Sa` share in our yeast data (0.45 % of typed reads).

## Experimental design of GSE146368

The series holds 96 samples, all of which are single-end NextSeq 500 runs of
BY4742 (an S288C derivative) grown in YPD to OD600 ≈ 0.8. The paper states 75 nt
reads, but the libraries of every strain carrying `dxo1Δ` or `rai1Δ`, together with
all 24 transcriptome libraries, were deposited at 85–86 nt. Every library is a
biological triplicate. In NAD captureSeq, "S" is the fully treated sample (+ADPRC)
and "N" is the negative control without ADPRC, so enrichment is always S versus N.

| Group | Strains | Samples | Runs | Reads (M) | FASTQ (GB) |
|---|---|---|---|---|---|
| Transcriptome (rRNA-depleted, random-sheared, NEBNext directional) | all 8 | 24 | 24 | 848 | 38.9 |
| NAD captureSeq, unfragmented (UF) | all 8 | 48 | 48 | 1,433 | 54.6 |
| NAD captureSeq, small fragmented (SF, 20–170 nt) | WT, `npy1Δ` | 12 | 12 | 302 | 10.2 |
| NAD captureSeq, large fragmented (LF, 170–350 nt) | WT, `npy1Δ` | 12 | 12 | 260 | 8.8 |

The eight strains are WT, `npy1Δ`, `dxo1Δ`, `rai1Δ`, `dxo1Δ npy1Δ`, `dxo1Δ rai1Δ`,
`rai1Δ npy1Δ` and `dxo1Δ rai1Δ npy1Δ`. The runs are SRR11236378 to SRR11236473,
and the UF libraries per strain are as follows.

| UF strain | S runs | N runs | Reads S+N (M) | Read length |
|---|---|---|---|---|
| WT | SRR11236402–404 | SRR11236405–407 | 60 | 75 |
| `npy1Δ` | SRR11236408–410 | SRR11236411–413 | 103 | 75 |
| `dxo1Δ` | SRR11236414–416 | SRR11236417–419 | 120 | 86 |
| `rai1Δ` | SRR11236420–422 | SRR11236423–425 | 152 | 85 |
| `dxo1Δ npy1Δ` | SRR11236426–428 | SRR11236429–431 | 343 | 86 |
| `dxo1Δ rai1Δ` | SRR11236432–434 | SRR11236435–437 | 218 | 86 |
| `rai1Δ npy1Δ` | SRR11236438–440 | SRR11236441–443 | 45 | 85 |
| `dxo1Δ rai1Δ npy1Δ` | SRR11236444–446 | SRR11236447–449 | 392 | 86 |

The GSM-to-SRR mapping for every run is reproduced by the ENA file report
`https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA610259&result=read_run&fields=run_accession,experiment_accession,sample_title,fastq_ftp&format=tsv`.

## Which datasets are useful

The selection follows from what each library can say about NAD 5′ ends at single-
nucleotide resolution:

| Priority | Datasets | Runs | Reason |
|---|---|---|---|
| 1 | UF `dxo1Δ rai1Δ npy1Δ`, S and N | 6 | This set is the least decapped and the deepest (392 M reads); the authors' YAAG motif comes from it. |
| 1 | UF WT, S and N | 6 | This set shows the physiological NAD landscape that a WT CAGE library can actually sample. |
| 2 | UF of the six single and double mutants, S and N | 36 | The gradient from WT to the triple mutant tests whether a candidate CAGE signature tracks NAD occupancy and supports replication of the priority-1 calls. |
| 3 | SF WT and `npy1Δ`, S and N | 12 | These libraries were sheared before capture and PNK-phosphorylated, so most fragments start internally, enrichment is weak (maximum fold change below 7) and true 5′ ends are diluted. |
| Skip | LF WT and `npy1Δ` | 12 | These libraries have the same drawback as SF with lower enrichment still. |
| Skip | Transcriptome, all strains | 24 | Random-primed, sheared libraries carry no 5′-end information; the processed count tables already give expression for any normalisation. |

The minimum useful set is therefore 12 runs (about 450 M reads and 16 GB of FASTQ),
and the full useful set is 48 runs, all unfragmented.

## File formats

### Per-sample wig files

GEO holds 32 `*_merge_genome.wig.gz` files, one per condition, each merged from the
three replicates and attached to the rep1 GSM record only (the other 64 GSMs have no
file). They are also bundled in `GSE146368_RAW.tar` (542 MB). Inspection of
`GSM4379033_wtS_UF` and `GSM4379048_dN_UF` showed the following properties:

- each file is a single `variableStep` track normalised to reads per million genome-mapped reads;
- the signal is read coverage, visible as runs of identical values along each read, and not a count of 5′ ends;
- there is only one track per file, so the strand is lost;
- coordinates are on the BY4742_Toronto_2012 assembly with chromosomes named `chr01`–`chr16`, `chr17` for the mitochondrion and `2-micron`, and chromosome lengths differ from sacCer3 (for example `chr07` ends at 1,086,808 against 1,090,940 for sacCer3 chrVII).

The wigs therefore cannot give strand-aware 5′-end positions in sacCer3
coordinates, and this branch does not use them beyond quick browsing.

### Series-level processed tables

The series supplementary files are tab-separated and indexed by an SGD feature
string on BY4742_Toronto_2012 gene models, where each mRNA is counted over its
5′ region from −120 to +65 around the start codon:

- `GSE146368_rawHits_<strain>_<lib>.txt.gz` hold htseq-count counts with columns `N1 N2 N3 S1 S2 S3` (14 files, for UF of all strains and SF and LF of WT and `npy1Δ`);
- `GSE146368_Stat_<strain>_<lib>_5UTRv.txt.gz` hold the DESeq2 S-versus-N results with columns `baseMean log2FoldChange lfcSE stat pvalue padj`;
- `GSE146368_rawHits_5UTRforAllmRNA_transcriptome_all_strains.txt.gz` and the matching `CDS` file hold transcriptome counts for all 24 libraries;
- `GSE146368_Stat_Transcriptome_<A>VS<B>_5UTRv.txt.gz` hold 12 DESeq2 strain comparisons of the transcriptome.

These tables give a per-gene list of NAD-enriched transcripts, which is useful as a
cross-check, but they say nothing about the exact start position.

### Raw FASTQ

The raw reads are on SRA and ENA as single-end FASTQ and are the only source of
exact NAD 5′ ends. The library was built by 3′-adapter ligation, reverse
transcription, C-tailing of the cDNA by TdT and ligation of an anchored 5′ adapter,
so the read starts at the RNA 5′ end preceded by a G stretch from the C tail. The
3′ adapter is `C NNNNNN AGATCGGAAGAGC…`, in which the six N look like a UMI. Two
read structures exist:

- the 75 nt runs (WT, `npy1Δ`, and every SF and LF library) begin with `GGG` or `GGGG` in about 93 % of reads, since PhiX spike-in was used and the tail is read through;
- the 85–86 nt runs (all strains carrying `dxo1Δ` or `rai1Δ`) were sequenced with the custom primer ending in three Gs, so most reads begin directly at the insert, but position 1 is still G in about 39 % of reads, which suggests leftover tail.

Because the tail is G, the boundary between tail and RNA is ambiguous for any RNA
that starts with G. NAD-RNAs always start with A, so their 5′ ends are unambiguous
once the leading G run is stripped. Soft-clip typing by pyselectal is meaningless
on these reads, because the clipped bases come from the TdT tail and not from the
cap.

## Preliminary observation from existing data

The mitochondrion offers a test that needs no new data. Mitochondrial transcripts
carry no m7G cap, so any CAGE signal at mitochondrial TSSs comes from another cap
or from uncapped background. Pooling the seven yeast libraries in
`results/tss/typed_ctss.tsv.gz` (`yeast_cage_glu_r1` was removed earlier for low
depth) gives the following composition:

| End type | Nuclear (%) | chrM (%) | chrM / nuclear |
|---|---|---|---|
| `1Sg` | 76.15 | 1.33 | 0.02 |
| `1Sa` | 0.45 | 3.28 | 7.3 |
| `1St` | 0.72 | 1.26 | 1.8 |
| `1Sc` | 0.32 | 0.50 | 1.5 |
| `M` | 14.04 | 90.73 | 6.5 |
| `other` | 1.96 | 2.56 | 1.3 |

chrM holds 616,138 typed reads, and 20,212 of them are `1Sa`, which is 20 % of all
`1Sa` reads in the genome. The excess is specific to a clipped A: `1St` and `1Sc`
rise much less. The ten chrM positions with the most `1Sa` reads include five
(73704, 58008, 6470, 46079 and 13357, all on the + strand) at which the read starts
on the final A of the mitochondrial promoter nonanucleotide `ATATAAGTA`, which is
the known mtRNAP initiation site. That is the class of RNA reported to be up to
60 % NAD-capped. Three further top positions start on a G followed by `TGCGC`,
which looks like processing sites and needs separate inspection. The mass of `M`
reads on chrM is most likely uncapped background. This observation fits the
hypothesis, but it does not yet show NAD, because mitochondrial transcripts could
reach the library by another route.

Step 1 of this branch (`scripts/01_chrM_end_types/`) turned the observation into a
figure with a matched control, and it found the following:

- restricted to Y|A starts, the only context in which a clipped A is possible, `1Sa` is 8.42 % of chrM reads against 0.60 % of nuclear reads, a 14-fold difference;
- the nuclear `1Sa` share rises with the length of the A run starting at +1, from 0.21 % at a lone A to 9.6 % at runs of five or more, which points to slippage on A tracts rather than to a cap;
- at a lone A the chrM share is still 7.5 %, which is 36-fold the nuclear value, so the chrM signal is not an A-run artefact;
- nuclear positions where `1Sa` exceeds 20 % of reads carry `CAAAA` or `TAAAA` and never a G at +3, so they lack the `YAAG` motif of NAD-capped nuclear starts, and most nuclear `1Sa` is therefore probably not NAD.

These results change the emphasis of step 5 of the plan. The nuclear test has to
look for an NAD-derived minority within `1Sa` and must stratify by A-run length,
otherwise the slippage signal swamps it.

## Analysis plan

1. Download the 12 priority-1 runs on the vavilov cluster into `nad/data/nad_captureseq/` with `nad/scripts/00_download/download.sh`, which reads the run list, ENA md5 sums and URLs from `nad/scripts/00_download/priority1_runs.tsv`, and add the 36 priority-2 runs only if step 5 gives a signal.
2. Strip the 3′ adapter with its UMI, strip the leading G run, keep reads of at least 18 nt (the authors' threshold for genome mapping), and record the length of the removed G run for each read. This step is `nad/scripts/02_trim/run.sh`: cutadapt removes the adapter (reads that do not reach it are kept with the UMI recorded as `NNNNNN`), and `umi_grun.awk` cuts the `C` and the six-base UMI, strips the G run and writes the read name as `<id>_g<G-run length>_<UMI>`, so that STAR and umi_tools both keep the two values. A check on 20,000 reads per library found that the insert starts with A in 50 % of WT reads and in 72 % of triple-mutant reads.
3. Map to sacCer3 with STAR, forcing end-to-end alignment at the 5′ end (`--alignEndsType Extend5pOfRead1`) so that the first aligned base is the RNA 5′ end, and keep MAPQ 255 as in Stage 1.
4. Count strand-aware 5′ ends per position, per replicate and separately for S and N, call positions or `distclu` clusters enriched in S over N with DESeq2, and require a genomic A at +1. A `YAAG` logo should reappear among the top triple-mutant calls and serves as the positive control for the pipeline.
5. Join the NAD-enriched positions to `results/tss/typed_ctss.tsv.gz` and compare their CAGE end-type composition (`1Sa`, `M` starting on A, `1Sg` starting on A, `2Saa`) against control positions that start with A, are expressed at matched CAGE depth, and are not enriched in NAD captureSeq, which mirrors the authors' background set; the test is a per-position Fisher or a matched permutation test.
6. Reverse the direction by asking whether positions dominated by `1Sa` (Stage 17) are enriched in NAD captureSeq S over N, and whether `1Sa` positions carry the `YAAG` motif more often than matched `1Sg` positions.
7. Measure, in the same NAD captureSeq libraries, the 5′-end signal at the chrM `ATATAAGTA` sites (`chr17` in the authors' assembly, chrM after remapping) to confirm that the mitochondrial `1Sa` sites are NAD 5′ ends.

## Caveats

The caveats below limit what any result can show:

- the CAGE libraries (Börlin et al. 2019) and the NAD captureSeq libraries come from different strain backgrounds and growth conditions, so the closest comparison is the glucose CAGE condition against YPD, and condition-specific TSSs will not overlap;
- NAD-RNAs are short and rapidly decapped in WT cells, so a WT CAGE library may barely sample them, and the triple mutant maps where NAD is incorporated rather than where it survives;
- a nuclear signal is capped at a few percent of reads per TSS, so the test needs pooled positions and depth-matched controls rather than single-site calls;
- whether reverse transcriptase adds T opposite nicotinamide riboside is unknown, and an in vitro check with a synthetic NAD-RNA put through the CAGE reverse-transcription step would settle which signature to expect;
- a `1Sa` excess at NAD sites would show a non-m7G adenosine cap, and FAD caps would look the same.

## Layout

The branch lives in its own folder so that it stays apart from the numbered
pipeline stages:

```text
nad/
  README.md      # this file
  scripts/       # one folder per step, each with a run.sh wrapper
    01_chrM_end_types/
  results/       # tables, one folder per step
  figures/       # figures, prefixed by step number
  reports/       # Quarto reports
  data/          # NAD captureSeq FASTQ, once downloaded
  logs/          # wrapper logs
```

The scripts source `scripts/utils/end_types.R` and `scripts/utils/plot_helpers.R`
from the main pipeline and read paths from `config/params.yaml`, as every stage
does. The project `.gitignore` ignores any folder named `results/`, `reports/`,
`data/` or `logs/`, so those folders under `nad/` are ignored as well.

## Status

Step 1 is done. Step 0 (`nad/scripts/00_download/`) and step 2 (`nad/scripts/02_trim/`) are scripted for the vavilov cluster. It is run with `bash nad/scripts/01_chrM_end_types/run.sh`, which
writes its tables to `nad/results/01_chrM_end_types/` and the figure to
`nad/figures/01_chrM_single_A_clip.{pdf,png}`. The report
`nad/reports/nad_first_finding.qmd` presents that figure and renders to a
self-contained HTML file. No NAD captureSeq data have been downloaded yet, and the
branch has not been given a number in `docs/stages.md`.
