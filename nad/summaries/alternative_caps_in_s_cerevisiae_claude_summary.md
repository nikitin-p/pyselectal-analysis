## Baseline: what's *not* in yeast

Worth stating first, since it shrinks the "alternative cap" space:

- *S. cerevisiae* mRNAs are **cap0 only** — no CMTR1/CMTR2 orthologue, so no cap1/cap2 2′-O-methylation.
- No **PCIF1** orthologue, so no m⁶Am.
- No **MePCE/Bin3** orthologue, so no γ-monomethyl phosphate cap: the MePCE ortholog present in fission yeast is absent from the *Saccharomyces cerevisiae* genome, implying U6 snRNA is stabilised by some other route (Páez-Moscoso et al., *Nat Commun* 2022, 13:1030). Contrast *S. pombe*, where U6 transcription by Pol III is accompanied by a 5′ γ-monomethyl phosphate cap added by the Bin3/MePCE-family enzyme Bmc1 (Porat et al., *NAR* 2023, 51:8805).
- **TMG (m2,2,7G)** caps are present, via Tgs1, on snRNAs/snoRNAs and telomerase RNA (TLC1) — the one well-established non-m⁷G nucleotide cap in budding yeast that isn't metabolite-derived.

So in practice, "alternative caps in *S. cerevisiae*" means: NAD⁺/NADH, plus a much thinner literature on FAD/dpCoA, plus the unmethylated-GpppN quality-control substrate class.

## 1. NAD⁺ / NADH — the well-characterised case

**Discovery and mechanism.** Walters et al. found 5′ NAD-RNA on subsets of nuclear- and mitochondrial-encoded mRNAs in *S. cerevisiae*, apparently produced co-transcriptionally: NAD-RNA is present on pre-mRNAs, and among mitochondrial transcripts only on those that are not 5′-end processed. Their estimate was on the order of 1–5% of a given mRNA population, with EGD2, GPP1 and RPL16b at ~10%, 7.5% and 5.5% respectively, and enrichment for mRNAs encoding components of the translation machinery, plus mitochondrially encoded ATP9, COX2 and the 15S/21S rRNAs. Mechanism is NCIN incorporation — Bird et al., *Nature* 2016, 535:444 showed yeast RNAP II can initiate with NAD⁺, NADH and dpCoA.

**The comprehensive landscape.** Zhang et al. redid this with the original NAD captureSeq protocol including the small-RNA fraction, and the picture changed substantially: NAD-RNAs are abundant (>1400 species), short (<170 nt), and mostly correspond to mRNA 5′-ends; modification is <5%; NAD incorporation occurs mainly at initiation by RNAP II using distinct promoters with a YAAG core motif; most NAD-RNAs are 3′-truncated; Rai1, Dxo1 and Npy1 guard against NAD-RNA at different cellular locations; and NAD-mRNAs are not translatable in vitro. Two specifics worth having:
- For nearly half of the species with sharp 5′ ends, the NAD-RNA 5′ end differs from the published canonical TSS, with transcript leaders both longer and shorter — in some cases by more than 100 nt.
- Promoter mutagenesis at the TDH3 locus showed the G at +3 modulates NAD incorporation specifically, without affecting overall expression.

Their conclusion is deliberately deflationary: NAD incorporation into nuclear yeast RNA appears largely accidental, reflecting NAD/ATP competition at initiation, with the cell disfavouring full-length NAD-RNA synthesis, decapping it via a multi-tiered compartmentalised machinery, and rejecting it from ribosomes.

**Mitochondria are the exception.** Yeast mtRNAP caps RNA with NAD⁺ and NADH more efficiently than nuclear RNAPs, with up to ~60% NAD⁺/NADH capping of yeast mitochondrial transcripts, set by promoter sequence at and upstream of the TSS and by intracellular NAD⁺/NADH levels (Bird et al., *eLife* 2018, 7:e42179; motif RA̲). Zhang et al. explicitly contrast this with the nuclear situation, noting the redox state of mitochondrial NAD caps varies with metabolic growth conditions, and that the three nuclear/cytosolic decapping enzymes lack mitochondrial targeting sequences and are absent from yeast mitochondrial proteomes — i.e. a dedicated mitochondrial deNADding activity is still unidentified.

**Decapping enzymes (four, non-redundant):**

| Enzyme | Class | Chemistry | Location |
|---|---|---|---|
| Rai1 | DXO/Rai1 | removes entire NAD moiety | nuclear, RNAP II-associated; acts first |
| Dxo1 | DXO/Rai1 | removes entire NAD moiety | nuclear + cytosolic |
| Npy1 | Nudix (NudC homolog) | cleaves NAD pyrophosphate → NMN + 5′-p-RNA | cytosolic (previously called peroxisomal) |
| Xrn1 | 5′→3′ exonuclease | releases intact NAD, then degrades | cytosolic |

Refs: Jiao et al., *Cell* 2017, 168:1015 (fungal Rai1/Dxo1 deNADding in vivo; *rai1Δ dxo1Δ* accumulates NAD-RNA); Zhang et al. 2020 for the in vivo hierarchy; Sharma et al., *Nat Commun* 2022, 13:889 — yeast Xrn1 is a deNADding enzyme that releases intact NAD and subsequently degrades the RNA, with a deNADding-deficient mutant that retains 5′-monophosphate exonuclease activity, linking mitochondrial NAD-capped RNA to overall NAD regulation.

**Translation — the live disagreement.** Zhang et al. found NAD-mRNAs untranslatable in a yeast extract across seven constructs with both long and short leaders. Sharma et al. 2023 reached a partly different conclusion in vivo: a dichotomy in which NAD-RNAs are detected with mitochondrial ribosomes but minimally on cytoplasmic ribosomes, indicating a propensity to be translated in mitochondria. These aren't strictly contradictory (in vitro cytosolic extract vs. in vivo mitoribosome association), but they're the two poles of the current debate.

## 2. FAD and dpCoA — enzymology only, no yeast profiling

Both are chemically plausible in yeast (yeast RNAP II accepts dpCoA as an NCIN; Bird 2016), and the relevant hydrolases exist, but there is no published transcriptome-wide FAD- or dpCoA-cap map for *S. cerevisiae*.

- DXO/Rai1-family enzymes remove FAD and dephospho-CoA caps ("deFADding", "deCoAping"), demonstrated structurally with *S. pombe* Rai1 bound to 3′-FADP (Doamekpor et al., *NAR* 2020, 48:6136; the study also used *Kluyveromyces lactis* Dxo1).
- Sharma et al., *NAR* 2022, 50:8807 — deFADding activity in human, yeast and bacterial 5′→3′ exoribonucleases (Xrn1/Rat1).
- Npy1 is a plausible but uncharacterised yeast deFADding enzyme: an earlier report that Npy1 dinucleotide hydrolysis extends beyond NAD to FAD warrants a thorough biochemical characterisation to define the set of metabolite caps it targets (AbdelRaheim et al., *Arch Biochem Biophys* 2001, 388:18).
- CapQuant (Wang et al., *NAR* 2019, 47:e130) included yeast RNA in its LC-MS/MS survey of 26 cap species, which is the closest thing to a quantitative yeast metabolite-cap inventory.

## 3. Unmethylated GpppN — yeast's own founding contribution

Not a metabolite cap, but a genuine non-canonical 5′ end and the system where 5′-cap quality control was defined. Rai1 (with Rat1) removes incomplete/unmethylated caps and triggers co-transcriptional degradation — Jiao et al., *Nature* 2010, 467:608; Xiang et al., *Nature* 2009, 458:784 (Rat1–Rai1 structure); Chang et al., *NSMB* 2012, 19:1011 (Dxo1 as a decapping/5′→3′ exonuclease hybrid). This is mechanistically the same surveillance machinery later shown to act on NAD caps.

## 4. Reported in other systems, *not* in yeast

Ap4A/NpnN caps, ADP-ribose caps, UDP-Glc/UDP-GlcNAc caps, and m⁷Gpppm⁶A have all been reported in bacteria, plants, or mammals but have no published *S. cerevisiae* evidence. Worth noting that Arabidopsis NudiX enzymes cleave Ap4A- and NADH-capped RNA in vitro, and the authors argue such caps should be searched for in plants — the same argument applies to yeast, where Npy1 is the obvious candidate to test.

## S. cerevisiae alternative-cap profiling datasets

| Dataset | Study | Method / content |
|---|---|---|
| **GSE90636** (GSM2409726–GSM2409737) | Walters et al., *PNAS* 2017, 114:480 | NAD captureSeq (+/− ADPRC) on BY4742, YEPD and synthetic media, triplicate, plus matched RNA-seq. Small-RNA fraction excluded. |
| **GSE146368** | Zhang et al., *Nat Commun* 2020, 11:5508 | The most comprehensive resource. Unfragmented + small-fragmented + large-fragmented NAD captureSeq, plus rRNA-depleted transcriptomes, for WT and all seven combinations of *npy1Δ*, *dxo1Δ*, *rai1Δ*. Ligation-based 5′ adapter → single-nucleotide 5′-end resolution. Source Data (64.8 MB zip) is attached to the paper. Companion proteomics: **PXD017893** (PRIDE). |
| **GSE217259** | Sharma et al., *Commun Biol* 2023, 6:406 | NADcapPro (copper-free SPAAC) and circNC in WT and *xrn1-H41A*, with side-by-side SPAAC-NAD-seq and standard RNA-seq for comparison. SPAAC-NAD-seq identified 268 transcripts, NADcapPro 769 in WT and 830 in xrn1-H41A. All source data are contained in the GSE217259 repository. |

Method papers you'd want alongside these: Winz et al., *Nat Protoc* 2017, 12:122 (NAD captureSeq); Vvedenskaya et al., *Mol Cell* 2018, 70:553 (CapZyme-seq — *E. coli*, but the enzymology generalises); Grudzien-Nogalska et al., *RNA* 2018, 24:1418 (NAD-capQ, bulk quantification).

## Practical caveats if you plan to reanalyse

**The three datasets are not directly comparable, and the discrepancies are methodological, not biological.** Walters excluded RNAs below ~200 nt; Zhang deliberately captured them and found the population is overwhelmingly short — only 5 of the 12 genes explicitly reported by Walters overlapped Zhang's unfragmented library. Sharma's NADcapPro then reported full-length polyadenylated NAD-RNAs that Zhang's protocol concluded were essentially absent. All three used the same organism and broadly the same medium. Size selection and adapter strategy dominate the result.

**The ADPRC background problem is real.** ADPRC has residual activity on m⁷G-capped RNA in eukaryotes; anti-m⁷G antibody depletion mitigates but does not eliminate it. NADcapPro addresses this with Dcs1 pretreatment — ADPRC showed noticeable activity on m⁷GpppA-capped RNA, and Dcs1 treatment specifically abolished that reactivity without affecting NAD-RNA detection. Any reanalysis of the 2017 or 2020 datasets should treat the −ADPRC control, not RNA-seq, as the baseline, exactly as the original papers did.

**Directly relevant to your TSS work.** The Zhang dataset is the single best resource for the promoter question, because the ligation-based 5′ adapter gives nucleotide-resolution 5′ ends rather than the heterogeneous ends from random priming. Two findings that bear on promoter-shape quantification: NAD-capped and m⁷G-capped species of the *same gene* frequently initiate at different positions, and NAD initiation is associated with a specific core-promoter variant (YA̲AG, +3 G critical) that sits inside the ordinary yeast Pol II initiator space. In practical terms this means a fraction of the initiation-site heterogeneity you'd measure at a yeast promoter reflects cap-identity-dependent start-site selection — a component that cap-trapping-based TSS data cannot resolve by construction, since it selects on the m⁷G diol. Whether that materially perturbs a shape score depends on the ~1–5% modification rate, which is probably small for most promoters, but the *dxo1Δ rai1Δ npy1Δ* libraries in GSE146368 give a direct way to test the magnitude.

**Reviews for orientation:** Sharma & Kiledjian, *JBC* 2022, 298:102171; Potužník & Cahová, *RNA Biol* 2024, 21:739. Neither is yeast-specific, but both cover the yeast enzymology accurately.
