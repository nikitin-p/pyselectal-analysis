## 1. Variants of the canonical m⁷G cap

Not alternative caps in the strict sense, but part of the same 5'-end landscape:

- **cap0 / cap1 / cap2** — 2'-O-methylation of the +1 and +2 nucleotides by CMTR1/CMTR2; principal role is self/non-self discrimination by RIG-I and IFIT proteins.
- **m⁷Gpppm⁶Am** — N6-methylation of a +1 adenosine by **PCIF1/CAPAM**. PCIF1 is the mammalian methyltransferase generating m⁶Am at the transcription start nucleotide; Pcif1-mutant mice are viable and fertile but show reduced body weight, and transcriptome analysis supports a transcript-stabilising role. Function remains contested — its effect on mRNA stability and translation has been reported inconsistently, and the FTO-as-eraser model (Mauer et al., *Nature* 2017) has not held up cleanly.
- **m⁷Gpppm⁶A** (i.e. m⁶A at +1 *without* 2'-O-methylation) — a genuinely new species, found by mass spec: CapQuant applied to human and mouse RNA revealed FAD, UDP-Glc, UDP-GlcNAc and m⁷Gpppm⁶A as new cap structures, along with surprisingly high proportions of caps lacking 2'-O-methylation.
- **TMG (m2,2,7G)** — trimethylguanosine hypermethylation by TGS1, on Sm-class snRNAs, snoRNAs and telomerase RNA; drives snRNP nuclear import.
- **Unmethylated GpppN** — an incomplete cap, cleared by DXO/Rai1 in "no-cap decay" 5'-end quality control.

## 2. Metabolite / cofactor caps (NCIN-derived)

These arise when RNA polymerase uses a metabolite as a **non-canonical initiating nucleotide** rather than ATP/GTP, so cap identity is coupled to metabolite pool size and to promoter sequence at/upstream of the TSS.

**NAD⁺ / NADH** — the best-characterised.
- Human discovery: Jiao et al., *Cell* 2017, 168:1015 — 5' NAD cap in human cells promotes decay via DXO-mediated deNADding.
- Targets: NAD⁺-capped RNAs are widespread in human cells and consist mainly of mRNAs and snoRNAs.
- Function is the opposite of bacteria: in bacteria the NAD cap protects the 5' end, whereas in mammals and plants it promotes decay; the NAD cap does not appear to support translation in mammalian cells.
- Decapping enzymes: **DXO** (removes the whole NAD moiety), **NUDT12** (cleaves the pyrophosphate, releasing NMN), **NUDT16**. NUDT12 is structurally and mechanistically distinct from DXO and targets different RNAs; its endogenous NAD-capped mRNA targets are enriched for proteins involved in cellular energetics, and NAD-cap levels shift under nutrient or environmental stress in a DXO- vs NUDT12-selective manner. NUDT16 was identified as a third mammalian deNADding enzyme.
- Mitochondria are a hotspot: human mitochondrial RNA polymerase caps RNA with NAD⁺ and NADH more efficiently than nuclear RNAPs, with up to ~15% NAD⁺ capping of human mitochondrial transcripts, set by promoter sequence and intracellular NAD⁺/NADH levels. This is the strongest case for caps acting as a metabolic sensor.
- Key refs: Grudzien-Nogalska et al., *Nat Chem Biol* 2019, 15:575; Sharma et al., *NAR* 2020, 48:6788; Bird et al., *eLife* 2018, 7:e42179; Grudzien-Nogalska et al., *Cell Rep* 2019 (NUDT12–BLMH cytoplasmic surveillance).

**FAD**
- FAD-capQ detected FAD-capped RNAs in human cells, restricted to short RNAs, and these are upregulated when DXO is absent — indicating DXO-dependent regulatory control. NUDT2 and NUDT16 also hydrolyse FAD-capped RNA in vitro, with NUDT16 regulating cellular FAD-cap levels.
- Ref: Doamekpor et al., *NAR* 2020, 48:6136.

**dpCoA (3'-dephospho-coenzyme A)**
- Present in mammalian RNA (Wang et al. 2019 CapQuant; Sharma et al. 2020). Seven mammalian Nudix proteins — NUDT2, 7, 8, 12, 15, 16, 19 — have deCoAping activity in vitro, but no dedicated cellular deCoAping enzyme has been pinned down in human cells.
- The most recent advance is largely plant-based: Hu et al. identified Arabidopsis NUDT11 as a dpCoA-specific decapping enzyme and used it to build quantification and transcriptomic profiling methods, showing dpCoA-RNAs exist across species with tissue- and condition-specific variation; in Arabidopsis they reach up to 15% of m⁷G-capped RNA abundance, have distinct TSSs, and associate with translating ribosomes (*Nat Biotechnol*, 2026; note there is an author correction to Extended Data Fig. 2i). Whether human dpCoA-RNA is comparably abundant or translatable is untested.

**UDP-glucose and UDP-GlcNAc** — detected in human and mouse RNA by CapQuant (Wang et al. 2019). No function, writer, or eraser established. Treat as "chemically demonstrated, biologically uncharacterised."

**ADP-ribose** — mechanistically distinct: added post-transcriptionally by transferases, not by NCIN initiation.
- RNA ADP-ribosylation occurs in human cells, mediated by TRPT1, PARP10, PARP11, PARP12 and PARP15, and reversed by TARG1, PARG and ARH3; ADPr-capped mRNA is protected against XRN1-mediated degradation but is not translated, and cellular stressors modulate ADPr-RNA levels. A family-wide macrodomain screen found PARP9macro1 and PARP14macro1 act as ADPr-RNA decapping domains in vitro, but PARG is the major ADPr-RNA decapping enzyme in HEK293 cells.
- Refs: Weixler et al., *NAR* 2022, 50:9426; Munnur & Ahel, *NAR* 2019, 47:5658.

## 3. Alarmone (dinucleoside polyphosphate) caps

**Ap4A** — diadenosine tetraphosphate was reported as a 5' RNA cap in human and rat cell lines; Ap4A-RNA levels are independent of the intracellular Ap4A concentration, and NUDT2 and DXO were identified as decapping enzymes. Ap4A-RNA co-fractionates with mRNA but is not translated, and is recognised as self by the innate immune system — it elicits no interferon response. Transfected Ap4A-capped RNAs are stabilised by the cap and do not interact with eIFs.
- Ref: Potužník et al., *Angew Chem Int Ed* 2024, 63:e202314951. Given Ap4A's alarmone status, a stress-responsive role is the obvious hypothesis but is not yet demonstrated.

## 4. Non-nucleotide, cap-like 5' structures

- **γ-monomethyl phosphate** — the oldest known non-canonical cap. It is a simple cap-like structure with no flanking nucleotide moiety, first identified on the Pol III transcript U6 snRNA; installed by **MePCE/BCDIN3**, also the 7SK capping enzyme (Jeronimo et al., *Mol Cell* 2007, 27:262). Protects Pol III transcripts from 5' degradation; on 7SK it is part of the 7SK snRNP that sequesters P-TEFb.
- **BCDIN3D** — the paralogue, 2'-O-monomethylphosphate capping of pre-miR-145, suppressing Dicer processing (Xhemalce et al., *Cell* 2012, 151:278); also described as a cytoplasmic tRNA^His 5'-monomethylphosphate capping enzyme (*NAR* 2020, 48:1572).

## 5. Viral caps encountered in human cells

- **HCV / FAD** — the HCV RNA-dependent RNA polymerase uses FAD as a non-canonical initiating nucleotide, giving a 5'-FAD cap at ~75% frequency, the highest metabolite-cap frequency observed in any organism; it is present in vivo on RNA from patient samples and provides protection against the innate immune response. (Sherwood et al., *Nature* 2023, 619:811; structural follow-up in *Nat Commun* 2025, 16, s41467-025-62609-w.)
- **HIV-1 / host NAD caps** — NAD captureSeq in HIV-1-infected cells showed four snRNAs (including U1) and four snoRNAs losing their NAD cap on infection, and the NAD cap decreases the stability of the U1/HIV-1 pre-mRNA duplex (Benoni et al., *ACS Chem Biol* 2024, 19:1243).

## Caveats

**Abundance numbers are soft.** Commonly quoted figures like 1–6% of mRNA for NAD caps come from methods with known artefacts. A significant limitation of NAD captureSeq in eukaryotes is residual promiscuity of the ADPRC enzyme toward m⁷G-capped RNA, partially mitigated by anti-m⁷G antibody depletion. NADcapPro (copper-free click) and circNC (intramolecular circularisation) were developed to resolve these limitations, and in yeast overturned prior conclusions: NAD-RNAs can be full-length and polyadenylated, NAD-capped and m⁷G-capped transcripts can have different TSSs, and NAD caps can be added after transcription initiation. The equivalent re-analysis has not been fully redone in human cells, so I'd treat pre-2023 human NAD-RNA target lists as provisional.

**Post-transcriptional capping exists.** The discovery of NAD-capped processed RNAs in human cells implies non-NCIN capping pathways — relevant because NCIN logic (metabolite pool × promoter sequence) can't explain caps on intronic snoRNAs.

**Relevance to 5'-end sequencing.** Cap-trapping CAGE selects on the m⁷G diol, so essentially none of section 2 or 3 is captured — non-canonically capped 5' ends are systematically invisible in standard CAGE/NET-CAGE. Since the circNC work found NAD and m⁷G TSSs can be genuinely distinct, at least in yeast, there may be a class of initiation sites that cap-trapping data structurally cannot see. If that matters for your promoter-shape work, the orthogonal methods to pair with CAGE are NADcapPro/circNC, CapZyme-seq, and CapQuant LC-MS/MS.

**Best entry points to the literature:** Potužník & Cahová, *RNA Biol* 2024, 21:739 (broadest recent survey); Sharma & Kiledjian, *JBC* 2022, 298:102171; Mattay, *WIREs RNA* 2022, 13:e1730; Wiedermannová, Julius & Yuzenkova, *R Soc Open Sci* 2021, 8:201979 (capping mechanisms beyond NCIN).
