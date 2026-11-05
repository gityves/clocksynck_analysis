# palumaa2026local
**Ref:** Palumaa, Cherry, Palta, Burns, bioRxiv (2026). DOI 10.64898/2026.03.24.713830. Type: data catalogue.
Lab: not Naef (authors' affiliations not checked).
**Source:** full text

## Question and design
Rhythmic transcriptome of 45 GTEx v11 tissues (14,886 samples) using per-tissue phases. **Sampling:** post-mortem,
whole day (GTEx).

## Method
CHIRAL (Talamanca 2023) with the same 12 clock genes: per-tissue local phase (LP) = CHIRAL tissue internal phase
(TIP); donor global phase (GP) = DIP. Rhythmicity by harmonic regression against LP (Bonferroni); KEGG enrichment.
No new inference method.

## Key results
- "Local phase inference identified a median of 5,747 rhythmic transcripts across tissues, including 5,531 across brain tissues"; claim that "prior donor-level approaches have substantially underestimated the extent of rhythmic transcription in the human brain".
- LP vs GP: "the estimated LPs were systematically phase-shifted relative to GP, and the strength of alignment was low to moderate (r range = 0.20-0.70, median r = 0.50, Fig. 1A)".
- >200 genes of Alzheimer's, Parkinson's, Huntington's and prion pathways rhythmic in brain (SNCA, PRKN, PSEN1, APP, HTT...).

## Not addressed / limits
- More rhythmic genes against a phase fitted within the same tissue is partly circular (LP is estimated from that tissue's data; regression on it inflates rhythmicity: cf. Anafi 2017 "too liberal" caveat). No null by permutation of LPs reported in the main text.
- LP-GP concentration r mixes biological tissue scatter with LP estimation noise (no technical-noise model), so it cannot be read as desynchrony.
- No disease contrast in donors (disease genes are annotated, not compared between patients and controls).

## Link to us
- **Intro claim:** tissue-level phases in GTEx reveal tissue-specific rhythmic programmes, incl. brain; cite briefly in internal desynchrony / circadian cycle, not as evidence of desynchrony. Our refs table overstated it as closest competitor: WRONG, it is the CHIRAL TIP step plus catalogue.
- **Relation:** background; the LP-GP r is the quantity clocksync estimates properly (coupling with Poisson noise separated).

## Ideas
- Contrast in the paper: their LP-GP r (0.20-0.70) vs our coupling c per tissue on GTEx; we separate measurement noise, they do not: brain low r may be noise (low amplitude) rather than desynchrony.
