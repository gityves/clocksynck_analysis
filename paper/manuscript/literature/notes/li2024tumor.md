# li2024tumor
**Ref:** Li, Hammarlund, Wu, ..., Hogenesch, Anafi, Meng, PNAS 121:e2311854121 (2024). DOI 10.1073/pnas.2311854121.
Type: method (CYCLOPS 2.0) + data + biology. Labs: Ron Anafi (Penn) / Qing-Jun Meng (Manchester).
**Source:** full text (SI not in export). The method paper for CYCLOPS 2.0 (supersedes hammarlund2025cyclops thesis).

## Question and design
Circadian rhythms of non-cancerous breast and breast tumours, and their clinical relevance.
- 43 mastectomy pairs (non-cancerous ≥4 cm away + tumour, same patient), Manchester; resection times recorded;
  26 time-stamped pairs with >20 M reads used. Luminal A 29, luminal B 3, HER2 2, TNBC 9.
- Public: GTEx breast (167), TCGA non-cancerous (106) and luminal A (193). Organoids (BMAL1-Luc), MCF-7 invasion.
- **Sampling:** surgery hours (Manchester recorded; TCGA inferred working day); GTEx autopsy.

## Method (CYCLOPS 2.0)
"We modified the CYCLOPS ( 25 ) neural network to accommodate explicit confounding variables, simultaneously learning confounder adjustments and a common circular structure" + optional supervision from time-stamped samples (weight set by 10-fold CV). Tumours projected onto the eigengene space of non-cancerous samples, tumour status as covariate. CYCLOPS magnitude (CMag) = radius, "a weighted sum of the amplitudes of the individually cycling seed genes".

## Key results
- Benchmark: "We benchmarked CYCLOPS 2.0 on actual, semisynthetic, and fully synthetic data with different temporal biases" (details in SI only).
- Time-stamped validation: "The CYCLOPS-predicted sample phases show a significant correlation with the known sample collection times of the 26 subjects (Corrcirc = 0.7, P < 0.005)".
- **Hospital hours:** "the sample phases assigned to the TCGA biopsies were concentrated during the presumed working day"; GTEx "more uniformly distributed, but still showed a significant bias".
- **Tumour vs paired normal phase (within patient):** "The sample phases assigned to the tumors and their noncancerous matches are poorly correlated (Fig. 3E)", but for high-CMag tumours "the circadian molecular phases assigned to high CMag luminal A samples generally match the phases assigned to their noncancerous pair".
- Clock organisation by subtype (Mantel Zstat): luminal A > luminal B (P < 0.01); organoids: luminal A "persistent (P < 0.01) but weakened rhythms", TNBC no consistent rhythms.
- Prognosis: high-CMag luminal A tumours, higher 5-y mortality "(P = 0.047, ANOVA)", survives age, metastasis, MKI67.
- BMAL1 knockdown / KL001 reduce MCF-7 invasion in 3D.
- Batch: "If different centers have different biases in collection time, that approach may remove true circadian signal."
- On CHIRAL: "Their approach assumes that tissues obtained from the same individual at the same time are at the same molecular phase."

## Not addressed / limits
- Tumour-normal phase discordance is shown on few TCGA pairs (Fig. 3E), with tumour phases from a projection onto normal eigengenes; no uncertainty, no specificity control (could be noise in low-amplitude tumours: the high-CMag concordance is consistent with this).
- Survival result marginal (P = 0.047) and flagged by the authors as needing an independent cohort.
- Still one phase per sample: needs the cycle covered (shift workers / time-stamped seeds).

## Link to us
- **Intro claims:** surgical cohorts sit in working hours (TCGA); tumour clocks persist but dampened (luminal A), lost in TNBC; strong-clock tumours = worse prognosis (direction!); tumour vs adjacent-normal desynchrony within patient is an open question with weak prior evidence (Fig. 3E).
- **Relation:** closest competitor on the tumour-normal question; CYCLOPS 2.0 = orderer with covariates + timed seeds. Their quote on CHIRAL motivates our coupling c (we do not assume same phase; we estimate it).

## Ideas
- Direct comparison: clocksync tumour-vs-adjacent coupling and offset in TCGA BRCA pairs, stratified by our amplitude (|z|) as analogue of CMag: do strong tumours couple to their normal? Prediction from their Fig. 4E: yes.
- Our weak-tumour low-coupling must be separated from estimation noise (they could not): the Poisson V_ik + c makes this explicit; a selling point.
- Check the CMag-survival link with our |z| in TCGA luminal A (pre-register).
