# cadenas2014loss
**Ref:** Cadenas et al., Cell Cycle 13:3282-3291 (2014). DOI 10.4161/15384101.2014.954454. Type: clinical (bulk tumours). Lab: Hengstler (IfADo Dortmund), Oster.
**Source:** full text (main; supplementary tables not in export)

## Question and design
Are 17 canonical clock genes (plus TIMELESS) prognostic in untreated node-negative breast cancer, and is co-expression of
clock genes (a proxy for a functional clock) lost in aggressive tumours? 766 patients in three Affymetrix cohorts:
Mainz (n = 200), Rotterdam (280), Transbig (286). Endpoint metastasis-free survival (MFS).
**Sampling:** one surgical tumour sample per patient, time of day unknown and not addressed; no normal comparison.

## Method (if any)
Univariate and multivariate Cox (dichotomised at median), Mann-Whitney vs clinical variables, Spearman correlation with
metagenes, and pairwise clock-gene correlations stratified by ER, HER2, grade and metastasis.

## Key results
- Prognosis: higher CLOCK, PER1, PER2, PER3, CRY2, NPAS2 and RORC expression go with longer MFS (HR < 1, FDR < 0.05); TIMELESS shorter (HR = 1.405).
- Multivariate: "only PER3 (HR = 0.66; P = 0.016) and RORC (HR = 0.42; P = 0.003)" independent of clinical parameters (abstract).
- Proliferation link: "Most clock genes (PER1, PER2, PER3, CRY2, RORC, BHLHE40, BHLHE41) showed a negative correlation (R < -0.3; P < 0.001) with the proliferation metagene", strongest CRY2 (R = -0.612), TIMELESS +0.674; clock genes track ER-positive status, low grade, low pT.
- Co-expression: PER2-PER3 R = 0.537 in ER+ but 0.054 in ER- tumours; "correlations between the expression of canonical circadian clock genes decrease in more aggressive tumor types".
- ARNTL vs PER anti-phase correlation (R = -0.2 for ARNTL-PER1) only in low-grade tumours.
- Residual distance from the PER2-PER3 regression as a loss-of-synchrony score: only a trend, HR = 1.72, P = 0.155.

## Not addressed / limits
- Clock-gene levels are not separated from proliferation, ER status, grade or immune content; most genes lose independence in multivariate models, as the authors say.
- Cross-sectional bulk: correlation across patients mixes clock coupling with composition and with unknown sampling time; no matched non-clock pairs, no gate.
- The loss-of-synchrony score is not significant; no tumour-vs-normal contrast.

## Link to us
- **Intro claim(s):** higher clock gene expression predicts longer MFS; loss of clock co-expression with progression. Our refs row is SUPPORTED but incomplete: the paper itself stresses ER/grade/proliferation and the synchrony-residual test is non-significant.
- **Relation to clocksync:** data we could use and already used.

## Ideas
- Done: A192 (fresh Mainz GSE11121 / Rotterdam GSE2034 / Transbig GSE7390, 684 tumours): clock score HR 0.809 per SD, not clock-specific against proliferation-matched pseudo-scores; falls to 0.883 (p 0.115) after G2M, ESR1, ERBB2; RORC stays 0.805 after G2M. A166/A167 on METABRIC registered.
- Not done: tumour-level clock coupling across ER+ vs ER- using the host-referenced readout, with matched non-clock gene pairs, to test the "correlation drop" at equal expression spread.
