# Eukaryotic RNAs Produced Uncapped

The m⁷G cap (m⁷GpppN) is the hallmark 5′ modification added co-transcriptionally to RNA polymerase II transcripts, discovered by Shatkin in 1976 [1] and subsequently shown to be recruited to the elongating Pol II complex via phosphorylation of the C-terminal domain (CTD) of the largest Pol II subunit [2, 3]. Because the capping machinery is physically coupled to Pol II, the polymerase identity of a transcript largely determines whether it receives an m⁷G cap. RNAs produced by **RNA polymerase I** and **RNA polymerase III** are the principal classes produced without an m⁷G cap.

## 1. RNA Polymerase I transcripts — 45S/47S pre-rRNA

The 45S (47S in mouse) pre-rRNA, which is processed into the 28S, 18S, and 5.8S mature rRNAs, is transcribed by Pol I in the nucleolus and **lacks an m⁷G cap**, bearing a 5′ triphosphate instead. Direct biochemical evidence in *Trypanosoma brucei* showed that Pol I–transcribed rRNAs (and Pol I–driven protein-coding genes in that organism) are uncapped, with cap addition only occurring subsequently via trans-splicing [4]. Genome-wide mapping of transcription start sites by ReCappable-seq, which specifically captures uncapped 5′-triphosphorylated ends that cap-based methods miss, confirmed that Pol I rRNA precursors carry triphosphate 5′ termini [5]. Reviews of Pol I transcription make no mention of cap addition as part of the initiation/elongation cycle, consistent with the absence of CTD-mediated capping enzyme recruitment that operates at Pol II [6, 7].

## 2. RNA Polymerase III transcripts (the majority)

Most Pol III products are synthesized with a **5′ triphosphate and no m⁷G cap** [8, 9, 10]. This class includes:

- **Pre-tRNAs** — transcribed with a 5′ triphosphate; the 5′ leader is subsequently removed by RNase P, generating a mature 5′ monophosphate, but no m⁷G cap is ever added [8, 9].
- **5S rRNA** — the only rRNA made by Pol III; produced with a 5′ triphosphate and no cap [10].
- **SINE/Alu transcripts and other short repetitive element RNAs** — transcribed by Pol III with 5′ triphosphates; their uncapped triphosphate ends are relevant to innate immune sensing by RIG-I-like receptors [8].
- **Vault RNAs, Y RNAs, RMRP/RPPH1 (RNase P RNA), and other small non-coding Pol III products** — fall under the same general rule of triphosphate, uncapped 5′ ends [11, 10].

### Important exception within Pol III: U6 snRNA and 7SK RNA

U6 snRNA and 7SK RNA are Pol III transcripts but are **not strictly "uncapped"** — they receive a distinct, non-m⁷G **γ-methylphosphate cap** (CH₃-ppp-N) added by the methylphosphate capping enzyme (MePCE/BICD-C3). This was established biochemically by Gupta & Reddy (1990), who showed U6 capping can be uncoupled from transcription and is catalyzed by a phosphate-specific methyltransferase [12], and by Shimba & Reddy (1994), who purified the ~130 kDa enzyme responsible for γ-monomethyl cap formation on U6 and related small RNAs [13]. Structural work by Yang et al. (2018) confirmed MePCE methylates the 5′ γ-phosphate of both U6 and 7SK and retains them through a slow product-release mechanism [14]. So while these RNAs lack the canonical m⁷G cap, they carry a different cap structure and should be distinguished from the genuinely uncapped Pol III transcripts above.

## What is NOT uncapped (common misconceptions)

- **Enhancer RNAs (eRNAs)** are Pol II transcripts and are **capped** (but typically non-polyadenylated and unstable). The original eRNA discovery paper by Kim et al. (2010) explicitly described "bidirectional, capped enhancer RNA transcription" at neuronal activity-regulated enhancers [15], and this is reiterated in subsequent reviews [16, 17]. eRNAs are unstable because they lack poly(A) tails and are rapidly degraded, not because they lack caps.
- **PROMPTs / uaRNAs** (promoter-upstream transcripts) are also Pol II products that are **capped and polyadenylated** but extremely short-lived due to rapid nuclear exosome-mediated degradation, as shown by Preker et al. (2008) [18] and characterized further in Preker et al. (2011) [19]. Their instability is a surveillance phenomenon, not a capping deficit.
- **Replication-dependent histone mRNAs** are Pol II transcripts that are capped but lack poly(A) tails (ending in a stem-loop instead).

## Summary

| RNA class | Polymerase | 5′ end status | Key evidence |
|---|---|---|---|
| 45S/47S pre-rRNA (→ 28S/18S/5.8S) | Pol I | Uncapped, 5′-ppp | [4, 5] |
| Pre-tRNAs | Pol III | Uncapped, 5′-ppp (→ 5′-p after RNase P) | [8, 9] |
| 5S rRNA | Pol III | Uncapped, 5′-ppp | [10] |
| SINE/Alu transcripts, vault RNAs, Y RNAs, RMRP | Pol III | Uncapped, 5′-ppp | [11, 8, 10] |
| U6 snRNA, 7SK RNA | Pol III | **γ-methylphosphate cap** (non-m⁷G) — a partial exception | [14, 12, 13] |

The organizing principle: **the m⁷G capping machinery is coupled to Pol II via CTD phosphorylation** [2, 3], so any transcript made by Pol I or Pol III is produced without an m⁷G cap, with U6/7SK being the notable case of a Pol III product that acquires a different, γ-methylphosphate cap post-transcriptionally.
References
1. Capping of eucaryotic mRNAs. Aaron J. Shatkin. Cell. 1976. doi:10.1016/0092-8674(76)90128-8.
2. mRNA capping enzyme is recruited to the transcription complex by phosphorylation of the RNA polymerase II carboxy-terminal domain. Cho EJ, Takagi T, Moore CR, et al.. Genes & development. 1997.
3. CTD-dependent and -independent mechanisms govern co-transcriptional capping of Pol II transcripts | Nature Communications. Melvin Noe Gonzalez, Shigeo Sato, Chieri Tomomori-Sato, et al.. Nature Communications. 2018.
4. The 5' end structure of transcripts derived from the rRNA gene and the RNA polymerase I transcribed protein coding genes in Trypanosoma brucei.. T. Bruderer, L. Tu, M. G. Lee. Molecular and biochemical parasitology. 2003. doi:10.1016/s0166-6851(03)00095-1.
5. Comprehensive determination of transcription start sites derived from all RNA polymerases using ReCappable-seq. Yan B, Tzertzinis G, Schildkraut I, et al.. Genome research. 2022.
6. Basic Mechanisms in RNA Polymerase I Transcription of the Ribosomal RNA Genes. Goodfellow SJ, Zomerdijk JC. Sub-cellular biochemistry. 2013.
7. RNA-polymerase-I-directed rDNA transcription, life and works. Russell J, Zomerdijk JC. Trends in biochemical sciences. 2005.
8. The molecular language of RNA 5′ ends: guardians of RNA identity and immunity. R. G. Avila-Bonilla, S. Macias. RNA. 2024. doi:10.1261/rna.079942.124.
9. RNA methyltransferases involved in 5′ cap biosynthesis. Byszewska M, Śmietański M, Purta E, et al.. RNA biology. 2014.
10. The Pol III transcriptome: Basic features, recurrent patterns, and emerging roles in cancer. Zhou S, Van Bortle K. Wiley interdisciplinary reviews. RNA. 2023.
11. RNA Polymerase III-Transcribed RNAs in Health and Disease: Mechanisms, Dysfunction, and Future Directions. Longjie Sun, Mingyue Chen, Xin Wang. International Journal of Molecular Sciences. 2025. doi:10.3390/ijms26125852.
12. Capping of U6 small nuclear RNA in vitro can be uncoupled from transcription.. Shashi Gupta, Ravinder Singh, R. Reddy. The Journal of biological chemistry. 1990. doi:10.1016/s0021-9258(19)38876-3.
13. Purification of human U6 small nuclear RNA capping enzyme. Evidence for a common capping enzyme for gamma-monomethyl-capped small RNAs.. Shigeki Shimba, R Reddy. Journal of Biological Chemistry. 1994. doi:10.1016/s0021-9258(18)99890-x.
14. Structural basis of 7SK RNA 5′ γ-phosphate methylation and retention by MePCE. Yuan Yang, C. D. Eichhorn, Yaqiang Wang, et al.. Nature chemical biology. 2018. doi:10.1038/s41589-018-0188-z.
15. Widespread transcription at neuronal activity-regulated enhancers. Tae-Kyung Kim, Martin Hemberg, Jesse Gray, et al.. Nature. 2010. doi:10.1038/nature09033.
16. Enhancer-derived RNA: A Primer. Liu F. Genomics, proteomics & bioinformatics. 2017.
17. Enhancer RNAs: a class of long noncoding RNAs synthesized at enhancers.. Tae-Kyung Kim, M. Hemberg, Jesse M. Gray. Cold Spring Harbor perspectives in biology. 2015. doi:10.1101/cshperspect.a018622.
18. RNA Exosome Depletion Reveals Transcription Upstream of Active Human Promoters. Pascal Preker, Jesper S. Nielsen, Susanne Kammler, et al.. Science. 2008. doi:10.1126/science.1164096.
19. PROMoter uPstream Transcripts share characteristics with mRNAs and are produced upstream of all three major types of mammalian promoters. Pascal J. Preker, Marianne S. Christensen, Eivind Valen, et al.. Nucleic Acids Research. 2011. doi:10.1093/nar/gkr370.