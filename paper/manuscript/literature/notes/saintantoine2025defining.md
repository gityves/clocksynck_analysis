# saintantoine2025defining
**Ref:** Saint-Antoine M, Anafi RC. "Defining and detecting global transcriptional amplitude in circadian gene expression." bioRxiv 2025 (posted 16 Nov 2025), DOI 10.1101/2025.11.14.688559. Preprint, not peer reviewed. Type: method + re-analysis of public timed human data. Lab: Anafi (U. Pennsylvania). Code: github.com/ranafi/Subject-Global-Amplitude.
**Source:** full text (bioRxiv v1, fetched 2026-11-24); Table 1 and the limitations section not read verbatim.

## Question and design
Do people differ in a "global" amplitude of their circadian transcriptome? Re-analysis of timed human series: skin EPIDERMIS_1 / DERMIS_1 (del Olmo 2022, GSE205155: 11 subjects, 7 samples), EPIDERMIS_2 / DERMIS_2 (Wu 2020, GSE139301 / GSE139300: 19 men, 4 samples at 6 h), ADIPOSE GSE87761, MUSCLE GSE108539, BLOOD GSE113883, GSE56931 and Möller-Levet 2013. **Sampling:** timed, per person.

## Method
Collective cosinor (F-test, BH q <= 0.05) selects rhythmic genes; per-subject cosinor; per-gene relative amplitudes percentile-ranked; a subject's "global amplitude" is the mean over rhythmic genes; row-shuffle permutation test (1,000 shuffles) for a subject-level trend.

## Key results
- "distributions of gene amplitudes tended to cluster by subject"; most subjects significant in three of four skin datasets (DERMIS_2: 8 of 19).
- Internal consistency: global amplitude from before-noon and after-noon peaking genes, R^2 0.843-0.855 (0.426 in DERMIS_2).
- Shared across layers of the same person: "R^2 = 0.649, p = 0.002 in EPIDERMIS_1/DERMIS_1; R^2 = 0.345, p = 0.008 in EPIDERMIS_2/DERMIS_2"; "subjects with higher epidermal rhythm strength also had higher dermal rhythm strength."
- Core-clock amplitude was "not strongly correlated" with global amplitude (abstract); sleep restriction did not track it in blood (R^2 0.019).

## Relevance for us
- The closest published counterpart of our person volume, in living timed tissue: a person-level amplitude shared by two skin layers. Ours spans 27 organs (GTEx deaths), uses the clock plane, has expression-matched and responsive companions and a glucocorticoid correlate (A360, A361).
- Caveats against equating them: their global amplitude is over all rhythmic genes (mostly outputs), and they report it is not strongly tied to core-clock amplitude; two layers of one biopsy share handling; no matched-gene companion.
- Priority: "persons differ in a shared rhythm amplitude" is reported here for skin; cite it, and claim the multi-organ, clock-specific version.

## Limits
Small cohorts; young, healthy; low sampling density in some sets; mechanism not identified [summary of the discussion, not verbatim].

## Our re-analysis (2026-11-24, post hoc; `results_preliminary.md`, "Skin shared amplitude against random genes")
Their statistic re-implemented from the repository: del Olmo layer R^2 0.608 (published 0.649). Random non-rhythmic,
expression-matched gene sets of the same size give r 0.699 (q95 0.720) against 0.780: most of the layer-sharing is a
non-rhythmic scale of the shared biopsy. Between two separate skin sites (GSE326530, not in their paper) the global
amplitude is not shared (r 0.13). Cite with this caveat.
