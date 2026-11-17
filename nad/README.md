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
which looks like processing sites; step 9 identifies them as the 5′ ends of the group II introns of COX1. The mass of `M`
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
3. Map to sacCer3 with STAR, forcing end-to-end alignment at the 5′ end (`--alignEndsType Extend5pOfRead1`) so that the first aligned base is the RNA 5′ end, and keep MAPQ 255 as in Stage 1. This step is `nad/scripts/03_map/run.sh`, which builds the index from the same `genome_fasta` and `annotation_gtf` that the CAGE pipeline used, so that NAD 5′ ends and CAGE CTSSs share coordinates, and which caps mismatches at 10 % of the mapped length because inserts can be as short as 18 nt.
4. Count strand-aware 5′ ends per position, per replicate and separately for S and N (`nad/scripts/04_ctss/run.sh` writes per-run tables with raw and UMI-collapsed counts, the genomic base at +1, the number of reads with a mismatched first base, and position-by-run matrices for the next part of the step), call positions or `distclu` clusters enriched in S over N with DESeq2, and require a genomic A at +1. A `YAAG` logo should reappear among the top triple-mutant calls and serves as the positive control for the pipeline.
5. Join the NAD-enriched positions to `results/tss/typed_ctss.tsv.gz` and compare their CAGE end-type composition (`1Sa`, `M` starting on A, `1Sg` starting on A, `2Saa`) against control positions that start with A, are expressed at matched CAGE depth, and are not enriched in NAD captureSeq, which mirrors the authors' background set; the test is a per-position Fisher or a matched permutation test.
6. Reverse the direction by asking whether positions dominated by `1Sa` (Stage 17) are enriched in NAD captureSeq S over N, and whether `1Sa` positions carry the `YAAG` motif more often than matched `1Sg` positions. This item is step 8 (`nad/scripts/08_1Sa_positions/`), and its result is described in its own section below.
7. Measure, in the same NAD captureSeq libraries, the 5′-end signal at the chrM `ATATAAGTA` sites (`chr17` in the authors' assembly, chrM after remapping) to confirm that the mitochondrial `1Sa` sites are NAD 5′ ends. This item is `nad/scripts/07_chrM_sites/`, and its result is described in the next section.

## Mitochondrial promoter starts in NAD captureSeq (plan item 7)

Item 7 refits the DESeq2 model of the S-versus-N step with its size factors
unchanged and adds a window of −3..+3 nt around every chrM promoter start, even
where the step-5 pre-filter dropped the position. The sites are all 14
occurrences of `ATATAAGTA` in chrM, taken with the final A as +1, together with
73704 (COX2) and 6548 (15S rRNA), which step 1 marked as promoter starts by the
`AAGTA` core alone. The cluster part (`run.sh`, `deseq_chrM_sites.R`) writes
`nad/results/07_chrM_sites/{motif_sites,controls,sanity}.tsv`, and
`07_join_cage.R` on the Mac adds the CAGE counts and draws
`nad/figures/07_chrM_sites.pdf`. The refit reproduces step 5 exactly (identical
size factors, r = 1 for log2FC on shared positions), and the chrM tRNA positions
are not enriched (median log2FC −0.94 in WT and −3.44 in the triple mutant).

The table gives every site with CAGE reads, at +1, as mean normalised S reads
and log2FC; N holds 0–3 reads at every site.

| Site | CAGE reads | CAGE `1Sa` (%) | WT S | WT log2FC | triple S | triple log2FC |
|---|---|---|---|---|---|---|
| 58008 (21S) | 120,081 | 2.7 | 728 | 12.0 | 653 | 10.5 |
| 73704 (COX2) | 12,626 | 63.5 | 240 | 10.4 | 57 | 7.0 |
| 6548 (15S) | 6,955 | 12.2 | 0 | – | 0 | – |
| 6470 (15S promoter) | 5,536 | 30.4 | 0.5 | ns | 0.2 | ns |
| 46079 (ATP9) | 1,509 | 37.2 | 1.7 | ns | 0.6 | ns |
| 13357 (COX1) | 1,084 | 29.4 | 35 | 7.6 | 0.9 | ns |
| 82329 | 125 | 64.8 | 179 | 10.0 | 3.2 | ns |
| 54840 (− strand) | 64 | 64.1 | 625 | 11.8 | 69 | 7.3 |
| 32217 | 17 | 41.2 | 37 | 7.7 | 4.1 | 3.2 |

The findings are the following:

- in WT six promoter starts are enriched in S over N (58008, 73704, 13357, 82329, 54840 and 32217), and at each of them the NAD captureSeq reads begin exactly on the final A, which is the position that carries the `1Sa` excess in CAGE;
- the remaining seven `ATATAAGTA` sites carry no NAD captureSeq reads and at most ten CAGE reads, which suggests promoters that are inactive in both growth conditions;
- the CAGE `1Sa` share does not predict the NAD signal everywhere, because 6470 (15S promoter, 30 % `1Sa`) and 46079 (ATP9, 37 % `1Sa`) carry at most three S reads per replicate, and 6548 carries none at +1;
- the triple mutant shows the same sites more weakly, which fits mitochondrial RNA being diluted in its S libraries by the nuclear NAD-RNAs that accumulate without the decapping enzymes;
- the three chrM processing sites that start on a G followed by `TGCGC` (16471, 13987 and 25343) carry no reads in S or N although CAGE sees thousands of reads there, so NAD captureSeq does not pick up every abundant chrM 5′ end, but these sites cannot test enrichment.

The strongest S signal at the 21S rRNA promoter does not sit on the initiating A.
At 58010, two nucleotides downstream, WT S holds 3,617, 1,614 and 1,222 reads
in the three replicates against 1,222, 650 and 500 at 58008, and the three N
libraries together hold three reads. Ensembl annotates the 21S rRNA gene (Q0158)
from 58009, so 58009–58010 is the processed 5′ end of the mature rRNA. A check of
the reads (`first_base_check.sh`) showed that 99.6 % of them start on the genomic
T at 58010 without a mismatch, and that their stripped G run is on average only
0.25 G longer than at 58008 (3.55 against 3.31 in WT), so at most about a fifth of
them can come from the G at 58009 whose base the G-run strip removed. In either
case the end is a G or a T and cannot carry NAD, and premature reverse
transcription stops before the cap are an unlikely cause, because the identical
`AGT` start of COX2 gives a read at +3 for only about 5 % of the reads at +1.
Item 7 therefore supports NAD at 58008 itself and treats the mature 21S end as a
non-NAD RNA that is nonetheless specific to S.

A raw-read check (`raw_read_check.sh`, cluster) asked whether the reads of 6470
and 46079 are lost in trimming or mapping. For every run and site it counts the
5′ ends in the BAM at +1, the trimmed reads whose bases 2–18 match the 18-mer
that starts at +1, and the raw reads that contain the 18-mer anywhere; 82329 is
left out because it shares its 18-mer with 54840. The counts below are pooled
over the replicates.

| Site | WT S: BAM / trimmed / raw | triple S: BAM / trimmed / raw |
|---|---|---|
| 6470 (15S promoter) | 1 / 1 / 1 | 2 / 2 / 4 |
| 46079 (ATP9) | 5 / 5 / 8 | 4 / 4 / 4 |
| 6548 (15S) | 0 / 0 / 381 | 0 / 0 / 205 |
| 58008 (21S) | 2,372 / 2,326 / 2,334 | 5,538 / 5,418 / 5,441 |
| 73704 (COX2) | 824 / 819 / 826 | 489 / 477 / 474 |
| 13357 (COX1) | 109 / 114 / 119 | 9 / 9 / 10 |
| 32217 | 114 / 112 / 114 | 32 / 31 / 31 |

The findings are the following:

- the reads of 6470 and 46079 are absent already from the raw FASTQ, so the pipeline loses nothing there, and these 5′ ends are missing from the NAD captureSeq libraries themselves;
- at the detected sites the three counts agree, which shows that trimming and mapping keep the reads that start on the promoter A;
- at 54840 the trimmed and raw counts (130,625 and 153,958 in WT S) exceed the BAM count (2,160) because the 18-mer occurs at seven loci and the reads that it starts are not uniquely mapped;
- the 18-mer of 6548 occurs in 381 raw WT S reads but starts none of them, and their 5′ ends lie upstream, between 6536 and 6542 for the minority whose prefix the check resolves and more than 12 nt upstream for most of them, so 6548 is covered by NAD captureSeq reads but is not a 5′ end in that assay.

The 18-mers of 6470 (`AATAAATAATAGTTTTAT`) and 46079 (`AATATATATAGTTTATGA`)
are as AT-rich as those of the detected sites 13357 and 32217, and both are
unique in sacCer3, so neither the sequence composition nor mapping ambiguity
explains their absence.

A follow-up asked whether fast processing of the primary transcript hides 6470
and 46079, since in CAGE the mature 15S end at 6546 holds 6.6 times the reads of
its promoter. For each of the eight `ATATAAGTA` promoters with at least 17 CAGE
reads, the check divides the CAGE reads at the promoter by the CAGE reads at all
5′ ends on the same strand within a window downstream of it. The window was 100,
300, 700 and 1,500 nt, and the CAGE and NAD counts come from
`results/tss/typed_ctss.tsv.gz` and the step-4 matrix. The table gives the
300-nt window.

| Promoter | NAD in WT | CAGE at promoter | CAGE downstream | Ratio |
|---|---|---|---|---|
| 73704 (COX2) | detected | 12,626 | 3,805 | 3.3 |
| 13357 (COX1) | detected | 1,084 | 423 | 2.6 |
| 58008 (21S) | detected | 120,081 | 58,133 | 2.1 |
| 54840 (− strand) | detected | 64 | 40 | 1.6 |
| 82329 | detected | 125 | 82 | 1.5 |
| 46079 (ATP9) | absent | 1,509 | 3,286 | 0.46 |
| 32217 | detected | 17 | 101 | 0.17 |
| 6470 (15S promoter) | absent | 5,536 | 65,334 | 0.08 |

The findings are the following:

- 6470 has the lowest ratio at every window, but 46079 ranks in the middle, and 32217 (ratio 0.16–0.63) and, at windows of 700 nt or more, 13357 (0.08–0.24) are detected with lower ratios, so processing does not separate the absent sites from the detected ones;
- the NAD signal does not follow the CAGE depth on chrM either, because 32217 has 17 CAGE reads and 114 WT S reads while 6470 has 5,536 CAGE reads and one WT S read, and the CAGE `1Sa` share of the absent sites (30 and 37 %) lies inside the range of the detected ones (3–65 %);
- the S libraries carry many S-specific starts on T and other non-A bases a few to a few tens of nucleotides downstream of the promoters, such as 73723 (a T, 3,180 WT S reads against 824 at 73704 and 30 in N), 13415–13417 near COX1 and 6535–6539 inside the 15S precursor, and these cannot be NAD 5′ ends.

With eight promoters the check cannot exclude a weaker link, and the most
likely remaining explanation is that 6470 and 46079 carry less NAD under the
growth conditions and in the strain background of the NAD captureSeq libraries,
which differ from those of the CAGE libraries. The data at hand cannot test
this, since a test needs CAGE and NAD captureSeq of the same RNA.

The earlier NAD captureSeq study of Walters et al. 2017 (PNAS 114:480, GEO
GSE90636) points the same way. Its abstract states that NAD-RNA is found only on
those mitochondrial transcripts whose 5′ ends are not processed, and the
literature summary in `nad/summaries/alternative_caps_in_s_cerevisiae_claude_summary.md`
lists ATP9, COX2 and the 15S and 21S rRNAs among them. ATP9 starts at 46079 and
the 15S rRNA precursor starts at 6470, so the two sites that are missing from the
Zhang et al. 2020 libraries would carry NAD in the Walters libraries. The gene
list has been taken from the summary and has not yet been checked against the
full text or the supplementary tables of the PNAS paper. The two studies differ
in size selection, since Walters et al. excluded RNAs shorter than about 200 nt,
so GSE90636 is the natural dataset for a direct test at 6470 and 46079.

The conclusion is that most mitochondrial promoter starts with a `1Sa` excess in
CAGE are NAD 5′ ends in NAD captureSeq, which agrees with the hypothesis that an
NAD cap is read by CAGE as a soft-clipped A. The absence of 6470 and 46079 is a
property of the NAD captureSeq libraries rather than of the pipeline, and
neither repeats, sequence composition nor fast processing explains it.

## Nuclear `1Sa` positions in NAD captureSeq (plan item 6, step 8)

Stage 19's position-level cascade calls only five nuclear positions
`1Sa`-dominated, so step 8 defines the `1Sa` positions by a share threshold
instead. A candidate is a nuclear Y|A position with a pooled CAGE depth of at
least 25 reads, which gives 34,089 candidates. A `1Sa` position is a candidate
whose `1Sa` share is at least 10 %, which gives 554 positions, and the thresholds
5 % and 20 % serve as sensitivity checks. Most of the 554 sit on A tracts (3, 83,
81, 259 and 128 at A runs of 1, 2, 3, 4 and 5 or more), as step 1 found. The NAD
signal of every candidate is read from the step-4 matrix, which was copied to
the Mac, and is scaled by the step-5 size factors, which the script recomputes
exactly. A position counts as NAD-detected when at least two of the three S runs
carry a read and the mean normalised S exceeds twice the mean normalised N.
`1Sa` positions are compared with the other candidates by a stratified label
permutation over the −1 base, the A-run length and the CAGE depth bin, as in
step 6. The script `08_1Sa_positions.R` runs directly with `Rscript` and has no
wrapper. It writes its tables to `nad/results/08_1Sa_positions/` and the figure
to `nad/figures/08_1Sa_positions.pdf`.

The checks pass. The recomputed signal reproduces the `meanS` and `meanN` of step 5
to within 5e-15 on all shared positions. With the NAD values shuffled within
strata, 1 of 60 tests over five seeds gives p < 0.05, and a planted detection of
30 % of the `1Sa` positions is found with p = 0.001.

| Strain | `1Sa` threshold | NAD-detected, `1Sa` positions | NAD-detected, matched controls | p |
|---|---|---|---|---|
| triple mutant | ≥ 5 % | 33.9 % | 35.4 % | 0.005 |
| triple mutant | ≥ 10 % | 32.4 % | 33.1 % | 0.50 |
| triple mutant | ≥ 20 % | 36.8 % | 32.1 % | 0.26 |
| WT | ≥ 5 % | 4.8 % | 4.3 % | 0.60 |
| WT | ≥ 10 % | 5.6 % | 3.6 % | 0.12 |
| WT | ≥ 20 % | 8.3 % | 3.3 % | 0.007 |

The findings are the following:

- the nuclear `1Sa` positions are not NAD 5′ ends more often than matched positions in the triple mutant, which maps NAD incorporation most fully, and at the 5 % threshold they are detected slightly less often;
- the only excess is in WT at the 20 % threshold, where 11 of 133 positions are detected against about 4 expected; these positions lie on A runs of 4 or more, and the excess is the one nominal hit among 12 tests, so it would not survive a correction for multiple testing;
- relative to the rate expected from its stratum, detection stays near 1 across `1Sa` share bins in the triple mutant, and it drops to 0.67–0.80 at A runs of 2 for `1Sa` shares between 2 and 20 %, which is the only run length at which `YAAG` can occur;
- among A-run-2 candidates, the G at +3 that completes `YAAG` is present at 1 of 83 `1Sa` positions (1.2 %) against 40 % of matched `1Sg`-majority positions (p = 0.001), so the nuclear `1Sa` positions avoid the motif of NAD-capped nuclear starts;
- the chrM `1Sa` positions serve as the positive control: 40 % of them are NAD-detected in WT and 67 % in the triple mutant, against 6 % and 32 % of the nuclear `1Sa` positions.

The conclusion is that nuclear `1Sa` reflects slippage on A tracts rather than
NAD caps, which agrees with step 1 and with the absence of a `1Sa` excess at the
nuclear NAD starts in step 6. The link between `1Sa` and NAD holds only on chrM.
The detection rule is lenient, since it counts 30–47 % of all candidates as
detected in the triple mutant, so it separates NAD starts from other positions
less sharply than a strict call would.

## Is the extra A on chrM specific to NAD? (step 9)

Step 7 left open whether the extra A of chrM CAGE reads marks NAD or any
non-m7G 5′ end, because three chrM 5′ ends that NAD captureSeq never sees
(16471, 13987 and 25343) carry `1Sa` in 28–46 % of their CAGE reads. Step 9
compares the extra base across classes of chrM 5′ ends within the same base at
−1, since a clipped base can never equal that base. `fine_clip_types.py` counts
the 5′ clip types of the seven subsampled CAGE BAMs without collapsing rare types
into `other` and writes `nad/results/09_chrM_clip_spectrum/fine_types.tsv.gz`, and
`09_chrM_clip_spectrum.R` classifies every chrM 5′ end with at least 50 CAGE reads
and writes `positions.tsv`, `spectrum_by_class.tsv`, `tests.tsv`,
`pol3_controls.tsv` and `nad/figures/09_chrM_clip_spectrum.pdf`. Both scripts run
on the Mac. The checks pass: the per-position totals and `1Sa` counts equal those
of `results/tss/typed_ctss.tsv.gz`, and labels shuffled within the −1 base give
p < 0.05 in 0 of 20 tests.

The classes are the step-7 promoter starts called in NAD captureSeq, the step-7
promoter starts that are not called (6470, 46079 and 6548), the first bases of
the chrM introns that begin with the group II consensus `GTGYG`, A starts that
are NAD-detected but are not promoters, A starts without S reads, and all other
5′ ends. Every promoter start has a T at −1, so the comparison below uses that
stratum, in which `1St` is impossible.

| Class (−1 = T) | Positions | `M` | `1Sa` | `1Sg` | A among single-base clips |
|---|---|---|---|---|---|
| promoter, NAD | 5 | 54 % | 45 % | 0.3 % | 92 % |
| promoter, no NAD | 3 | 44 % | 27 % | 0.7 % | 97 % |
| group II intron 5′ end (16471) | 1 | 30 % | 46 % | 23 % | 66 % |
| A, NAD-detected, not a promoter | 3 | 94 % | 2.7 % | 1.3 % | 54 % |
| A, no S reads | 15 | 65 % | 0.4 % | 21 % | 16 % |
| other 5′ ends | 93 | 98 % | 0.2 % | 0.8 % | 8 % |
| nuclear 5′-ppp control SCR1 | 1 | 20 % | 0.15 % | 75 % | 0.2 % |

The findings are the following:

- the three chrM 5′ ends with a `1Sa` excess that NAD captureSeq does not see are exactly the 5′ ends of the three group II introns of COX1 (aI1, aI2 and aI5γ), whereas the 5′ ends of the group I introns of COX1 and COB carry 0–6 CAGE reads each, so these reads most likely come from excised lariats, in which reverse transcription can cross the 2′–5′ branch and copy the branch A;
- the lariat ends receive A together with T and G (A is 50–66 % of their single-base clips), while the promoter starts receive almost only A (92–97 %), so within the T stratum the two sources are told apart by the share of G, which is about 1 % at promoters and 23 % at 16471;
- A starts that are NAD-detected but lie inside a transcript, such as 58011–58060 in the 21S rRNA, carry almost no extra A (2.7 % against 45 % at the promoters, p = 0.048), which fits read starts inside NAD-RNAs, for example premature reverse-transcription stops in the S library, and shows that the extra A belongs to the primary 5′ end itself;
- SCR1, a nuclear Pol III transcript with a 5′ triphosphate that starts in the same `T|A` context, receives a G in 75 % of its reads and an A in 0.15 %, so the extra A at the chrM promoters is not a generic property of a triphosphorylated A start (NME1, which starts in `C|A`, gives a mixture of G, T and A, so the nuclear control rests mainly on SCR1);
- the promoters 6470 and 46079, which NAD captureSeq does not detect, carry the same almost pure A signature as the detected promoters, which agrees with Walters et al. 2017 having found NAD on ATP9 and the 15S rRNA;
- the A starts without S reads that receive a G in 80–100 % of their reads (13628, 22246, 37737) carry almost no `M`, which looks like a sequence difference between the CEN.PK mitochondrial genome and sacCer3 rather than a 5′ structure.

The conclusion is that the extra A on chrM marks the primary transcription start
and not merely a non-m7G 5′ end. It separates the promoter starts from internal
NAD-detected starts, from ordinary processed ends, and from a nuclear 5′-ppp start
in the same context, and it differs in spectrum from the lariat ends, which are
its only other source on chrM. NAD remains the only known structure at those
starts that explains it, but the test cannot exclude a mitochondria-specific
behaviour of 5′-ppp RNA, which only a spike-in of synthetic 5′-ppp and NAD-RNA
would settle.

`09b_fine_clip_figure.R` draws `nad/figures/09_nad_fine_clip_types.pdf` from the
same fine types and from the run-A positions of step 6. It shows that the NAD
signature is `1Sa` together with its T-extended relatives `2Sta` and `3S·ta`,
which form the same extra-T ladder as `1Sg`, `2Stg` and `3S·tg` at m7G caps,
whereas `2Saa`, `2Sca` and `2Sga` are depleted on chrM, and that no A-ending clip
type exceeds the nuclear background at the nuclear NAD starts.

## How many NAD TSSs produce a `1Sa` read in CAGE? (step 10)

Step 10 (`nad/scripts/10_nad_tss_1Sa/10_nad_tss_1Sa.R`, Mac) counts the NAD TSSs
that carry at least one `1Sa` read in pooled CAGE and compares the count with an
expectation. The NAD TSSs are the step-5 positions enriched in S over N with an A
at +1, unmasked in the nucleus and of any biotype on chrM, and each set is also
cut to its top quartile by mean normalised S. Only `Y|A` starts are counted,
because a clipped A is impossible after a genomic A. The expectation sums, over
the NAD TSSs, the share of non-enriched nuclear `Y|A` positions with at least one
`1Sa` read in the same stratum of −1 base, A-run length and log2 CAGE depth. One
chrM position, 58011 (`T|AAAAA` at 13,253 reads), lies in a stratum without
background positions and takes the background of its −1 base and A-run at depth
bins of at least 10. The script writes `summary.tsv` and `positions.tsv.gz` to
`nad/results/10_nad_tss_1Sa/` and the figure `nad/figures/10_nad_tss_1Sa.pdf`.

| Strain | Compartment | Set | NAD TSSs | In CAGE | Median CAGE depth | `Y\|A` in CAGE | Observed ≥ 1 `1Sa` | Expected | p |
|---|---|---|---|---|---|---|---|---|---|
| WT | nucleus | all | 212 | 184 | 756 | 167 | 118 | 112.4 ± 4.2 | 0.31 |
| WT | nucleus | top quartile | 53 | 48 | 4,831 | 46 | 38 | 37.4 ± 1.7 | 0.48 |
| WT | chrM | all | 32 | 27 | 16 | 13 | 9 | 4.9 ± 1.0 | 0.06 |
| WT | chrM | top quartile | 8 | 6 | 95 | 6 | 5 | 3.6 ± 0.8 | 0.30 |
| triple | nucleus | all | 3,367 | 3,099 | 242 | 2,215 | 1,153 | 1,169.0 ± 17.8 | 0.68 |
| triple | nucleus | top quartile | 842 | 806 | 867 | 600 | 422 | 423.8 ± 8.3 | 0.54 |
| triple | chrM | all | 51 | 44 | 12 | 13 | 6 | 3.5 ± 0.8 | 0.15 |
| triple | chrM | top quartile | 13 | 12 | 72 | 5 | 4 | 3.3 ± 0.6 | 0.41 |

The findings are the following:

- the nuclear NAD TSSs are not missing from CAGE, since 87–96 % of them are present and the top quartile has a median depth of 870–4,800 reads, which fits mRNA 5′ ends that are mostly m7G-capped at the same position;
- the number of nuclear NAD TSSs with at least one `1Sa` read matches the expectation almost exactly in both strains, so the NAD-RNA share of a nuclear TSS leaves no trace in `1Sa`, in agreement with steps 6 and 8;
- chrM holds only 32 (WT) and 51 (triple mutant) NAD TSSs, many of which are read starts just downstream of 58008 and 73704 with an A at −1, so only 13 of them are `Y|A` starts in CAGE, and five WT NAD TSSs (73721, 73722, 73724, 13416 and 26573) have no CAGE reads at all;
- the chrM excess (9 against 4.9 and 6 against 3.5) is not significant, and chrM background positions themselves carry at least one `1Sa` read in 20–30 % of cases at moderate depth, so a count of positions with at least one `1Sa` read cannot separate NAD starts on chrM, where the per-read share of steps 7 and 9 does.

## Browser tracks of NAD 5′ ends (step 11)

Step 11 (`nad/scripts/11_tracks/11_nad_bigwigs.R`, Mac) turns the step-4
raw-read matrix into stranded bigWigs in `nad/results/11_tracks/`, so that the
NAD captureSeq signal can be browsed beside the CAGE end-type tracks of Stage 21.
No BAM is opened, so the step runs where the matrix is. Each strain and strand
gets three tracks:

- `S` is the mean over the three S replicates of reads divided by the step-5 size factor of each run;
- `N` is the same quantity for the three N replicates;
- `diff` is S minus N at every position, and it is negative wherever N exceeds S.

The size factors are those of step 5, recomputed as in step 8, so the `S` and `N`
tracks equal the `meanS` and `meanN` columns of `<strain>_S_vs_N.tsv.gz` at every
tested position, which the script checks, and the difference track shows what
DESeq2 compared. Reads per million are deliberately not used, because the N
libraries hold mostly rRNA and background, and depth scaling would make S and N
incomparable at mRNA starts. A log2 ratio track is not written either, because
single-nucleotide counts are mostly between 0 and 3 and their ratio is noise. A
bigBed per strain (`nad_<strain>_calls.bb`) marks the step-5 calls with padj below
0.05 and log2FC above 1, red where the base at +1 is A and grey otherwise.

`scripts/13_tracks/02_make_hub.R` adds these files to the Stage 21 hub as one
composite per strain and strand. Within a composite the three tracks share one
y scale (`autoScale group`), so that S and N heights are directly comparable,
and only `S` is switched on by default. The files must be uploaded into the same
Hub Space folder as the other Stage 21 files. In the triple mutant the N
libraries carry more normalised signal than S over the rDNA and other masked
ncRNA (the total of the difference track is negative), so a window holding rDNA
should be read on N as well as on S.

## Caveats

The caveats below limit what any result can show:

- the CAGE libraries (Börlin et al. 2019) and the NAD captureSeq libraries come from different strain backgrounds and growth conditions, so the closest comparison is the glucose CAGE condition against YPD, and condition-specific TSSs will not overlap;
- NAD-RNAs are short and rapidly decapped in WT cells, so a WT CAGE library may barely sample them, and the triple mutant maps where NAD is incorporated rather than where it survives;
- a nuclear signal is capped at a few percent of reads per TSS, so the test needs pooled positions and depth-matched controls rather than single-site calls;
- whether reverse transcriptase adds T opposite nicotinamide riboside is unknown, and an in vitro check with a synthetic NAD-RNA put through the CAGE reverse-transcription step would settle which signature to expect;
- a `1Sa` excess at NAD sites would show a non-m7G adenosine cap, and FAD caps would look the same;
- the literature summaries in `nad/summaries/` disagree on whether cap trapping captures NAD caps: `alternative_caps_cage_capture_claude_summary.md` correctly states that periodate reads any free 2′,3′-*cis*-diol, which the nicotinamide riboside of an NAD cap carries, whereas the caveats of `alternative_caps_in_human_claude_summary.md` claim that CAGE selects on the m7G diol and cannot see metabolite caps, which is wrong, and the `1Sa` excess at the chrM NAD starts agrees with NAD-RNAs reaching the CAGE libraries;
- enrichment in S over N does not prove an NAD cap, because the N libraries are nearly empty at mRNA starts and any 5′ end that reaches S looks enriched, as the mature 21S rRNA end at chrM:58010 shows;
- an NAD cap can sit only on an A, yet 22 % of the enriched unmasked nuclear positions of step 5 in WT and 23 % in the triple mutant start on C or T, and G starts are invisible after the G-run strip, so at least about a fifth of the enriched positions are not NAD 5′ ends and the A-starting set used in step 6 is diluted by an unknown non-NAD share as well.

## Layout

The branch lives in its own folder so that it stays apart from the numbered
pipeline stages:

```text
nad/
  README.md      # this file
  scripts/       # one folder per step, each with a run.sh wrapper
    00_download/ 01_chrM_end_types/ 02_trim/ 03_map/ 04_ctss/
    05_deseq/ 06_cage_end_types/ 07_chrM_sites/ 08_1Sa_positions/
    09_chrM_clip_spectrum/ 10_nad_tss_1Sa/ 11_tracks/
  results/       # tables, one folder per step
  figures/       # figures, prefixed by step number
  reports/       # Quarto reports
  summaries/     # literature summaries on alternative and absent caps (from a colleague)
  data/          # NAD captureSeq FASTQ, BAMs and the Ensembl GTF (on vavilov)
  logs/          # wrapper logs
```

The scripts source `scripts/utils/end_types.R` and `scripts/utils/plot_helpers.R`
from the main pipeline and read paths from `config/params.yaml`, as every stage
does. The project `.gitignore` ignores any folder named `results/`, `reports/`,
`data/` or `logs/`, so those folders under `nad/` are ignored as well.

## Status

Steps 0–11 have been run. Step 1 (`nad/scripts/01_chrM_end_types/run.sh`, Mac)
writes `nad/figures/01_chrM_single_A_clip.{pdf,png}`, which the report
`nad/reports/nad_first_finding.qmd` presents. Steps 0 and 2–5 run on the vavilov
cluster through `sbatch`, and the S-versus-N test of step 5 runs in the conda
env `r_nad` (`envs/r_nad.yml`). Step 5 found 300 enriched unmasked positions in WT
and 4,402 in the triple mutant, and the `YAAG` positive control passed. Step 6
(`nad/scripts/06_cage_end_types/`, Mac) found no `1Sa` excess in CAGE at the
nuclear NAD starts of the triple mutant (0.30 % against 0.36 % at matched
controls, ratio 0.82, p = 0.017). Step 7 (`nad/scripts/07_chrM_sites/`, cluster
and Mac) and step 8 (`nad/scripts/08_1Sa_positions/`, Mac) are described in
their own sections above. The raw-read check of step 7 showed that 6470 and
46079 are missing from the raw NAD captureSeq reads, and a comparison of the
eight CAGE-supported chrM promoters found that fast processing of the primary
transcript does not explain their absence either. Step 9
(`nad/scripts/09_chrM_clip_spectrum/`, Mac) showed that the extra A on chrM
marks promoter starts, including 6470 and 46079, and that its only other chrM
source is the 5′ end of the group II intron lariats. Step 10
(`nad/scripts/10_nad_tss_1Sa/`, Mac) found that NAD TSSs with at least one `1Sa`
read in CAGE are as frequent as expected from matched background in the nucleus
and not significantly more frequent on chrM. Step 11
(`nad/scripts/11_tracks/`, Mac) wrote the NAD captureSeq bigWigs and added them
to the Stage 21 hub. The branch has not been given a
number in `docs/stages.md`.
