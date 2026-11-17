## What cap-trapping actually selects for

The chemistry is three steps, and none of them touches the guanine, the methyl, or the triphosphate:

1. **Reverse transcription** first (random-primed in standard CAGE, so cDNA synthesis runs toward the 5′ end).
2. **NaIO₄ oxidation** at pH ~4.5, 4 °C. Periodate cleaves the C2′–C3′ bond of any **vicinal 1,2-*cis*-diol**, converting it to a dialdehyde. In CAGE the diol groups at both ends of the mRNA are oxidised in preparation for the biotin hydrazide reaction, and the 5′ end can only be biotinylated if it carries an intact cap structure.
3. **Biotin hydrazide** forms a hydrazone with those aldehydes; RNase ONE then digests any RNA not protected in an RNA:cDNA duplex (removing the 3′ biotin and, if the cDNA fell short, the 5′ biotin too); streptavidin capture keeps only full-length hybrids.

So the selection criterion is precisely: *a periodate-oxidisable vicinal diol at the 5′ terminus, on a molecule whose cDNA reached that terminus.* Nothing else.

Why the m⁷G cap qualifies is exactly what the diagram shows. Because the cap guanosine is joined 5′→5′, its own 5′-OH is consumed by the triphosphate bridge, leaving C2′-OH and C3′-OH both free and *cis* on the ribofuranose ring. Every other nucleotide in the transcript has its 3′-oxygen locked in a phosphodiester, so no internal position is reactive. This is old and well-established chemistry: the 5′-5′ pyrophosphate linkage and the free 2′-3′ *cis*-hydroxyl group in the m⁷G ribose are what permit isolation of cap-adjacent sequences — capped 5′-terminal oligonucleotides were originally purified on borate resin that binds specifically to 2′-3′ *cis*-OH groups in the m⁷G ribose.

## m⁷G-derived caps: all captured, all indistinguishable

This part of your question has a clean answer. **Every known modification of the m⁷G cap leaves the diol untouched**, because they all fall in one of two places:

| Cap species | Where the modification sits | m⁷G ribose 2′,3′-diol | Cap-trapped? |
|---|---|---|---|
| m⁷GpppN (cap0) | — | intact | ✔ |
| m⁷GpppNm (cap1) | 2′-O of the **+1** ribose | intact | ✔ |
| m⁷GpppNmNm (cap2) | 2′-O of +1 and +2 | intact | ✔ |
| m⁷Gpppm⁶Am | N6 of the +1 adenine + its 2′-O | intact | ✔ |
| m⁷Gpppm⁶A | N6 of the +1 adenine only | intact | ✔ |
| m2,2,7GpppN (TMG) | N2, N2 of the **cap guanine base** | intact | ✔ |
| **GpppN (unmethylated)** | nothing — no N7 methyl | intact | **✔** |

The last row is the one worth dwelling on. Cap-trapping reads the *ribose*, not the *methyl*. An unmethylated GpppN cap — the DXO/Rai1 "no-cap decay" substrate — is chemically identical from periodate's point of view and is captured with equal efficiency. CAGE therefore does not distinguish canonical mRNA caps from incompletely capped surveillance intermediates, and it captures TMG-capped snRNAs/snoRNAs too (which is why they show up in CAGE libraries as strong non-mRNA peaks). If you need m⁷G specificity, you need a different reagent entirely — anti-m⁷G antibody, eIF4E affinity, or a differential DcpS/Cap-Clip treatment.

## Where your premise needs correcting

The claim that non-m⁷G caps lack diols is true for only a minority of them. Most metabolite caps carry a vicinal diol on the capping moiety, for exactly the same structural reason as m⁷G: the moiety is attached through its 5′-phosphate, so its own 2′/3′ (or equivalent) hydroxyls stay free. This is not speculative — it is the basis of boronate affinity electrophoresis: boronyl groups form complexes with 1,2-*cis*-diols occurring naturally at the 3′-end of RNA and also in the nicotinamide riboside of NAD-modified RNA at the 5′-end, retarding NAD- or FAD-modified RNAs relative to unmodified ones. APB electrophoresis distinguishes capped RNA containing a vicinal-diol moiety — m⁷G, NAD⁺, NADH, FAD, or NpₙN — from uncapped RNAs. Note that m⁷G sits in that list alongside NAD: they are the *same* chemical class as far as diol reagents are concerned.

Here is the same contrast at the level of the two capping moieties themselves.## Full audit of the non-canonical caps

In each case, the capping moiety attaches through a phosphate/pyrophosphate, and the question is only whether *that moiety* has two adjacent free hydroxyls.

| Cap | Capping moiety | Vicinal diol? | Chemically cap-trappable? |
|---|---|---|---|
| NAD⁺ | nicotinamide riboside | 2′,3′-*cis* on the nicotinamide ribose | ✔ |
| NADH | dihydronicotinamide riboside | same ribose, unaffected by reduction | ✔ |
| FAD | FMN / riboflavin | 2′,3′,4′ of the acyclic ribityl chain | ✔ (weaker) |
| Ap4A, other NpₙN | adenosine 5′-polyphosphate | 2′,3′-*cis* on the capping adenosine | ✔ |
| ADP-ribose | terminal ribose | 2′,3′ free | ✔ |
| UDP-Glc / UDP-GlcNAc | glucose / GlcNAc | 2,3,4-OH on the pyranose | ✔ |
| **dpCoA** | pantetheine-4′-phosphate | none — one isolated C2-OH, a thiol, a *gem*-dimethyl quaternary C3 | **✘** |
| **γ-monomethyl phosphate** | a methyl on the γ-phosphate | no sugar at all; the +1 ribose has its 3′-O esterified | **✘** |
| **2′-O-methylphosphate (BCDIN3D)** | methyl on a 5′-monophosphate | same — no terminal free ribose | **✘** |
| 5′-ppp / 5′-p (uncapped) | — | 3′-O in phosphodiester | ✘ |

The FAD row deserves a footnote: FAD-RNA is retarded on APB gels relative to triphosphorylated RNA, but less strongly than NAD-RNA — the ribityl polyol binds boronate more weakly than a rigid *cis*-diol on a furanose ring, and periodate would cleave it at multiple points rather than opening a single ring. So FAD capture is plausible but likely inefficient and I'd treat it as untested rather than established.

Note also the internal consistency check: dpCoA is the one metabolite cap for which the field's detection chemistry is **thiol-based** (dpCoA tagSeq uses a maleimide–thiol reaction, and Hu et al. 2026 needed a dedicated decapping enzyme, AtNUDT11, to build a profiling method). Nobody uses diol chemistry on dpCoA — because there is no diol.

## So why *doesn't* CAGE report these caps?

Four reasons, in descending order of importance. Only the last is chemical.

**1. CAGE gives no read-level cap identity.** This is the real answer. If an NAD-capped transcript is oxidised, biotinylated, captured and sequenced, the resulting CTSS is a normal-looking 5′ end at a normal-looking coordinate. There is no barcode distinguishing it. The failure mode is not *loss* — it is *silent misattribution*. Any NAD, Ap4A, ADPr or UDP-sugar-capped 5′ ends that survive the protocol are simply counted as canonical TSSs. That matters for your promoter-shape work in a specific way: it does not create missing signal, it creates a small contaminating mixture with a *different* start-site distribution, since NAD initiation in yeast has its own core-promoter preference (YA̲AG) and frequently uses a distinct TSS from the m⁷G species of the same gene.

**2. Abundance.** At 1–5% modification for most transcripts, the non-canonical population contributes a tail, not a peak. Standard depth won't separate it from ordinary TSS heterogeneity.

**3. Size and integrity selection.** CAGE library prep involves RNase ONE digestion, gel/bead size selection, and (in many workflows) an RT that must traverse the full 5′ UTR. Yeast NAD-RNAs are predominantly short and 3′-truncated; short, fragmented species are disproportionately lost. NET-CAGE compounds this — it profiles chromatin-associated nascent RNA, so anything already decapped by Rai1 co-transcriptionally never enters the library.

**4. Chemical exclusion.** Genuinely applies only to dpCoA caps, γ-monomethylphosphate caps (U6, 7SK), 2′-O-methylphosphate ends, and uncapped 5′-p/ppp ends. These are invisible to cap-trapping by construction, and no depth of sequencing recovers them.

## If you want to test this on your own data

The clean experiment is a differential decapping pretreatment before oxidation, exploiting the fact that each enzyme destroys a different diol:

- **DcpS or Cap-Clip** removes m⁷G, collapsing the canonical population. What still cap-traps afterward is diol-bearing but non-m⁷G. This is exactly the logic ONE-seq uses in reverse: after yDcpS treatment and 3′-adaptor blocking with dideoxycytidine, NAD-RNA but not m⁷Gppp-RNA oligos were retained by boronic acid beads.
- **NudC or DXO** removes NAD; comparing CAGE ± NudC gives you the NAD-attributable fraction of your CTSS signal directly.
- **3′-end blocking with a ddC adaptor** eliminates the competing 3′ diol so the 5′ signal is unambiguous.

Applied to a NET-CAGE-style protocol, that would give you a per-CTSS estimate of non-canonical cap contamination, which is the number you'd actually want before deciding whether it perturbs a promoter shape score.
