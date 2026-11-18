# gibbs2014epithelial
**Ref:** Gibbs et al., Nature Medicine 2014 (20:919-926). DOI 10.1038/nm.3599. Type: biology (cell-type-specific knockout, immunity). Lab: Loudon / Ray (Manchester).
**Source:** full text

## Question and design
Which lung cell type carries the circadian gating of pulmonary neutrophil recruitment? Bmal1 flox crossed with a club-cell Cre (codon-improved Cre, iCre, under the CCSP/Scgb1a1 promoter) gives "Ccsp-Bmal1-/-"; littermate Bmal1 flox/flox controls without Cre. A myeloid knockout (LysM-Cre) is used as a comparator, adrenalectomy and CXCL5-null mice for mechanism. Readouts: Per2-luciferase ectopic lung slices, clock genes in whole lung and liver, in situ hybridisation of the bronchiolar epithelium at CT0 and CT12, neutrophil counts, cytokines, GR ChIP, bacterial infection. **Sampling:** timed (CT0 / CT12 for most readouts; 24 h profile for Cxcl5; ZT for infection experiments); n 3-12 per group depending on the panel; clock-gene time series of whole lung are in Supplementary Fig 3a (not in the text read).

## Method (if any)
None. Conventional qPCR, bioluminescence imaging, ChIP.

## Key results
- **Driver and targeting.** "expression of iCre was limited to the lung and trachea (Fig. 2b)" (RT-PCR over cerebrum, cerebellum, gut, heart, uterus, liver, bladder, kidney). Recombination efficiency in club cells and recombination in other lung cell types are not quantified in the text read; specificity rests on the promoter and on the in situ readouts.
- **Targeted cells' clock lost (ex vivo, bronchioles).** PER2 bioluminescence in bronchioles of control lung slices oscillates, "(25.8 ± 0.17 h, mean ± s.e.m., Fig. 2e)", and is disrupted in the knockout (32.5 ± 0.8 h, Fig 2f). In situ: Bmal1 and Nr1d1 in bronchial epithelia at CT0 / CT12 and loss of Nr1d1 induction (n = 4 and 3 per genotype).
- **Whole tissue, non-targeted cells measured.** "Loss of rhythmicity was confined to bronchioles, as clock gene (Bmal1, Nr1d1, Per2 and Cry1) oscillations were similar in both genotypes in the whole lung (Supplementary Fig. 3a)." So the lung clock-gene rhythm in bulk is retained when the bronchiolar clock is cut. Numbers (amplitude, phase, n, time points) are in the supplement and not in the text read; "similar" is the authors' word.
- Earlier work cited by the authors made the opposite argument: "Our earlier studies showed that these club cells express clock proteins and are essential for the maintenance of circadian rhythmicity in lung" (their ref 25, ex vivo slices).
- **Systemic drive.** After adrenalectomy, pulmonary Per2, Bmal1 and Nr1d1 show "similar high-amplitude, time-of-day" changes in adrenalectomized and intact animals (Supplementary Fig 10b): the lung clock does not need the adrenal glucocorticoid rhythm.
- **Other tissues and behaviour.** "Gene targeting had no impact on oscillations of clock genes in the liver (Supplementary Fig. 3b)" and none on behavioural rhythms (Supplementary Fig 3c). Corticosterone rhythm unchanged in the knockout ("persistent glucocorticoid rhythm", Fig 6).
- **Other cell type.** Myeloid Bmal1 knockout (LysM) "showed suppressed expression of clock genes Nr1d1 and Nr1d2" in bronchoalveolar macrophages, yet neutrophil and cytokine responses were like controls: "suggesting that cells of monocytic origin are not involved in the circadian gating we observed". Whole-lung clock genes of the LysM knockout are not reported in the text read.
- Output: club-cell Bmal1 loss makes Cxcl5 non-rhythmic and fourfold elevated (24 h profile); three- to sixfold more neutrophilia after LPS.

## Not addressed / limits
- The tissue-level clock-gene result is a supplementary panel, described in one sentence; no amplitude or phase estimate, no test of whether whole-lung phase shifted.
- Whole lung mixes alveolar epithelium, endothelium, immune and stromal cells; which of them sustain the retained rhythm is not tested (no sorted or single-cell clock data). The retained rhythm could be cell-autonomous in those cells, coupled, or driven by systemic cues (SCN-driven behaviour and temperature, adrenal ruled out).
- Bronchiolar loss is shown on ectopic lung slices (ex vivo) and in situ, not on isolated club cells; CCSP-Cre recombination efficiency not reported.
- Whether the whole-lung series of Supplementary Fig 3a is free-running or in LD is not stated in the text read.

## Link to us
- **Intro claim(s):** within a tissue, the clock does not depend on one cell population being rhythmic: **evidence against pacing by club cells** (bulk lung rhythm survives club-cell loss), compatible with autonomy of other lung cells or with systemic drive; it cannot separate autonomy from coupling because bronchiolar cells are targeted, and neighbours are not read separately. It does show that a cell-type clock can be removed without removing the organ-level clock-gene rhythm, as in our brain-dead donors (tissue clock amplitude kept) and in our "no conserved leader" result.
- **Relation to clocksync:** background and supports (indirect). Human lung single-cell atlases could test whether club or epithelial cells share more or less host state than endothelium.

## Ideas
- Request or locate Supplementary Fig 3a source data (whole-lung clock gene series of both genotypes) to quantify residual amplitude and phase; no GEO accession is given in the text read.
- The 32.5 h versus 25.8 h bronchiolar PER2 periods are a within-tissue contrast that a coupling model (club cells entrained by neighbours) would not predict for loss of rhythm; worth comparing with Bmal1-stopFL lung data if available (cf. Koronowski for liver).
