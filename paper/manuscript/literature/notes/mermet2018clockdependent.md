# mermet2018clockdependent
**Ref:** Mermet, Yeung, Hurni, Mauvoisin, Gustafson, Jouffe, Nicolas, Emmenegger, Gobet, Franken, Gachon, Naef, Genes Dev 32:347 (2018). DOI 10.1101/gad.312397.118. Type: biology. Lab: Naef (EPFL).
**Source:** full text (main text with Methods)

## Question and design
Are promoter-enhancer contacts static or rhythmic along the 24 h cycle, and do they matter for circadian transcription and behaviour? Mouse liver and kidney, 4C-seq at two
genes (Cry1, peak ZT20; Gys2, peak ZT08), WT, Bmal1 KO and a new Cry1 intronic enhancer deletion mouse (Cry1Δe, 300 bp), smRNA-FISH for burst fraction, locomotor activity in constant darkness.
**Sampling:** timed. Main contrast at two times: "livers of wild-type mice collected at ZT08 and ZT20 (n = 4 per time point)." Time course: liver "every 4 h in wild" type and Cry1Δe, "type and Cry1Δe (n = 3 per time point)."

## Method (if any)
4C-seq with locally weighted multilinear regression, signed -log10 P differential contacts between conditions, harmonic regression for oscillation of contacts; smRNA-FISH on Cry1 pre-mRNA for burst fraction and intensity. No gene-expression model selection.

## Key results
- Contacts oscillate: Cry1 promoter contacts a +26 kb intronic and a -7 kb upstream element more at ZT20 (P < 5.5e-17 at the peak); Gys2 promoter contacts exon 8 more at ZT08 (P < 8.7e-5), absent in kidney where Gys2 is flat. Reciprocal 4C confirms.
- Clock dependent: in Bmal1 KO the loops are static, Cry1 locked closed, Gys2 locked open: "the promoter at comparable frequencies at ZT20 and" ZT08. Cry1 effect is "Bmal1 knockout most likely reflects an indirect effect" via REV-ERB and ROR activities (REV-ERBα and RORγ bind the Cry1 enhancer; Rev-Erbα/β low and Rorγ high in Bmal1 KO).
- Enhancer deletion (Cry1Δe): contact frequencies lower at all times and oscillation compromised; peak Cry1 mRNA reduced; burst fraction 2.2-fold higher at ZT20 than ZT08 in WT, reduced by 28 % in Cry1Δe at ZT20, burst intensity unchanged (burst frequency, not size).
- Behaviour: Cry1Δe period is 15 min shorter than WT littermates (P < 1.1e-5, t-test), compared with the 1.2 h shortening of Cry1 loss of function; the deletion is a hypomorph.

## Not addressed / limits
- Two loci only (Cry1, Gys2); two time points in the main comparison, 4 h resolution time course with n = 3; liver and kidney; mouse.
- Bmal1 KO result for Cry1 is indirect (ROR/REV-ERB), so the data do not show that BMAL1 itself opens the loop at Cry1; Rev-erb is not manipulated here.
- Does not address tissue-specific phase differences, cross-organ synchrony or disease.

## Link to us
- **Intro claim(s):** refs_circadian_desynchrony.md claim: clock-dependent chromatin topology (Bmal1 KO, Rev-erb) modulates rhythmic transcription. PARTLY OVERSTATED: Bmal1 KO and a Cry1 enhancer deletion are the manipulations; Rev-erb is only a binding factor, not perturbed. Clock-dependent oscillatory promoter-enhancer loops are supported at two loci.
- **Relation to clocksync:** background (mechanism for clock-gene and clock-controlled gene rhythms); no inference method.

## Ideas
- Replace the Rev-erb attribution in the reference table by a statement of what was done (4C-seq, Bmal1 KO, Cry1Δe, burst frequency).
- Context reference for the same lab: Yeung 2018 (tissue-specific loops, already a card); Xu 2016 (Rev-erbα dynamically modulates chromatin looping) is cited here and is the Rev-erb manipulation paper we may mean.
