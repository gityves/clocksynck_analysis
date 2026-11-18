# manella2021liverclock
**Ref:** Manella, Sabath, Aviram, Dandavate, Ezagouri, Golik, Adamovich, Asher, Nat Metab 2021 (3:829-842). DOI 10.1038/s42255-021-00395-7. Type: biology. Lab: Asher (Weizmann).
**Source:** full text

## Question and design
Male C57BL/6 mice (3-4 months), 12 h light-dark, fed ad libitum (AL) or only in the light phase (daytime feeding, DF) for 30 d; controls Alb-Cre (AlbCre) and hepatocyte-specific Bmal1 knockout (BLKO, liver clock off). Do clocks of different peripheral tissues respond alike to DF, and does the liver clock shape the other tissues? Tissues for clock genes: liver, white adipose tissue (WAT), lung, quadriceps, kidney, heart (qPCR of Arntl, Dbp, Nr1d1, Per1, Per2, Cry1). Transcriptome (MARS-seq, 3' end) of liver, WAT, lung in AlbCre and BLKO, AL and DF; data GSE159135.
**Sampling:** timed, around the clock. qPCR: every 2 h for 24 h, n = 2 mice per time point and condition (tissues collected at death). RNA-seq: per-ZT averages (time-point spacing and n not stated in the methods text; check GSE159135). Second cohort for metabolic phenotyping (glucose, glycogen at ZT11 and ZT23).

## Method (if any)
JTK_CYCLE for rhythmicity and phase; cosine fit for amplitude (peak-to-trough, linear scale); two-stage JTK plus harmonic regression for differential rhythmicity of transcriptomes (rhythmic in AL only, DF only, both). Polar summary per tissue: mean phase shift of rhythmic clock genes (angle) and relative amplitude change (length).

## Key results
- Tissue-specific response of the clock to DF. "the expression of liver-clock genes was 12-h phase-shifted on DF" and WAT similar; lung "barely affected"; kidney and heart "approximately 4–8-h phase shift"; quadriceps lost rhythmicity. Summary: "While the clocks of some tissues (liver, WAT) were strongly phase-shifted, other tissues were only partially shifted (kidney), showed no response (lung) or even exhibited loss of rhythmicity (quadriceps)". Per-tissue values are given graphically (Fig. 1h polar plot, Ext. Data Fig. 1); no table of shifts or amplitudes in the text.
- Amplitude falls under DF in every tissue, including the phase-inverted ones. "We also noticed a marked reduction in the amplitude of clock genes on DF in all tissues." (Results). The amplitude drop is not quantified per tissue in the text (polar vector length, 1 = no change).
- Liver transcriptome: comparable number of rhythmic genes in AL and DF, "similar amplitude", mean 12 h shift. WAT clock inverted but its rhythmic transcriptome mostly advanced (median -5.5 h, spread over 12 h). Lung: DF shrinks the rhythmic transcriptome with a strong amplitude reduction, phases unchanged.
- Cue interaction is named explicitly. "The DF paradigm is, by its nature, a conflicting signal paradigm" (Discussion); signals interact "beyond the additive effects of each signal" and "it can explain differences in phase shift (direction and extent), clock-gene-specific effects and loss of rhythmicity".
- Loss of rhythmicity is left open between amplitude reduction and desynchrony: "One option is that statistical loss of rhythmicity is a result of amplitude reduction." Another is cell or mouse desynchrony.
- Liver-clock knockout leaves other tissues' clock phases and DF-induced amplitude reduction essentially unchanged (Fig. 3a-c), but alters their transcriptome phases on DF, with impaired glucose homeostasis.
- Interpretation of tissue differences: tissues "can differ in their sensitivity or their level of exposure to a given signal", i.e. a per-tissue weight on the feeding cue.

## Not addressed / limits
- Only AL vs DF (inversion); no feeding-free or arrhythmic arm, so amplitude retention without a feeding rhythm is not tested here.
- n = 2 per time point; per-tissue shifts and amplitude ratios only graphical; JTK phase of a handful of genes, not a population-level amplitude.
- No model; the additive-versus-interacting statement is verbal. Lung and quadriceps phase and amplitude are not related to a weight parameter.
- Amplitude drop under DF not separated from tissue-level desynchrony (pooled mice per time point).

## Link to us
- **Intro claim(s):** supports peripheral clocks differ in compliance with feeding and that conflicting light/SCN and food cues lower clock-gene amplitude. Shift per tissue (liver 12 h, WAT inversion, kidney/heart 4-8 h, lung about 0, quadriceps arrhythmic) is the best single-experiment ordering of a tissue weight b for the behavioural drive.
- **Relation to clocksync:** data we could use (GSE159135 liver, WAT, lung RNA-seq, AL vs DF, AlbCre and BLKO). Feeds A297 (liver AlbCre AL vs DF) and the conflicting-cue model: monotone pattern of phase shift with tissue (lung < kidney/heart < liver/WAT) with amplitude lowest where conflict is strongest fits weight near 1 cancels, but quadriceps loss of rhythm also fits tissue-level desynchrony. It does not test removal of the cue.

## Ideas
- Extract per-tissue phase shift and amplitude ratio from Supplementary Data 1 (qPCR) and GSE159135 (liver, WAT, lung) and fit b per tissue; test whether amplitude ratio is minimal at intermediate shift (kidney, heart) as the cancellation model predicts, versus lowest at full inversion.
- Cross-check with Wefers 2018 human muscle (shift 2.7 h, amplitude 0.23): quadriceps here also loses rhythm, same tissue.
- Compare with Greenwell 2019 (arrhythmic feeding keeps the liver clock) in the same pipeline.
- Their amplitude drop in liver under DF (clock genes) is consistent with our A297 liver amplitude ratio for AlbCre AD-lib vs DF; check it is not only a JTK artefact with n = 2 per time point.
