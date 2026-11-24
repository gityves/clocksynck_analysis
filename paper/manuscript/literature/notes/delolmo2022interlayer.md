# delolmo2022interlayer
**Ref:** del Olmo M, Spörl F, Korge S, Jürchott K, Felten M, Grudziecki A, de Zeeuw J, Nowozin C, Reuter H, Blatt T, Herzel H, Kunz D, Kramer A, Ananthasubramaniam B. "Inter-layer and inter-subject variability of diurnal gene expression in human skin." NAR Genom Bioinform 2022;4(4):lqac097. DOI 10.1093/nargab/lqac097. Type: clinical / data / method (mixed models). Lab: Kramer, Herzel, Ananthasubramaniam (Charité / HU Berlin). Data: GEO GSE205155 (public since 2022-12-12, 154 samples, Agilent 4x44k raw TXT). Code: github.com/bharathananth/skin-data-analysis.
**Source:** full text (PMC9803873), first ~40 % verbatim, the rest from a model summary (marked [summary]); GEO record via the GEO API (2026-11-24).

## Question and design
How much do diurnal gene expression rhythms vary between skin layers and between subjects? "Eleven healthy volunteers (six males, five females, aged 20–30 years in 2011)", intermediate chronotypes; "one skin punch (3 mm diameter) was taken from the lower back" "every 4 h over a 24 h period (8:00, 12:00, 16:00, 20:00, 0:00, 4:00, 8:00)" (seven samples); epidermis and dermis separated from each biopsy. **Sampling:** timed, full day, per person; time converted to internal time by mid-sleep (MSFsc); DLMO in six subjects only.

## Method
Per-gene linear mixed cosinor models with subject and layer as random effects; amplitude and phase variances by error propagation.

## Key results
- "10% of the expressed genes showed such diurnal rhythms at the population level, of which only a third differed between the two layers."
- "Amplitude and phases of diurnal gene expression varied more across subjects than layers, with amplitude being more variable than phases." (abstract)
- "Expression amplitudes in the epidermis were larger and more subject-variable, while they were smaller and more consistent in the dermis." (abstract)
- Layers in phase: mean phase difference -0.08 h for common diurnal genes [summary].
- No per-subject correlation of amplitude between layers is reported in the part read (Saint-Antoine & Anafi 2025 compute it on these data: R^2 0.649).

## Relevance for us
- The person-level clock amplitude ("person volume", A312, A336 part 2) in living, timed tissue: subjects differ in rhythm amplitude more than layers do. Two layers of one biopsy, so a shared person / handling scale is not separated (no companion genes).
- Data public and timed (7 samples per person and layer): a candidate living test of the person volume with true amplitude (A364 on GSE326530 was not decidable: 4 samples, power 0.54).

## Limits
11 subjects, young, healthy; two layers of one punch; 4 h sampling; sex and subject not separable [summary]; variance estimates from two layers "likely less accurate" [summary].
