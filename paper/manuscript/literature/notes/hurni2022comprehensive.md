# hurni2022comprehensive
**Ref:** Hurni C, Weger BD, Gobet C, Naef F, PLoS Genet 18(8):e1009903, 2022. DOI 10.1371/journal.pgen.1009903. Type: biology (mouse, mechanism). Lab: Naef.
**Source:** full text

## Question and design
Whether nuclear export and cytoplasmic degradation shape rhythmic mRNA profiles in mouse liver. Male C57BL/6J mice (8-12 weeks), 12:12 light-dark, night-restricted feeding for 3 days before sampling. Nuclear and cytoplasmic fractions (plus unfractionated control), total (rRNA-depleted) and poly(A) RNA-seq, bulk liver. **Sampling:** timed and full 24 h: ZT0, 4, 8, 12, 16, 20, n = 2 mice per time point (12 mice); the Results text speaks of two consecutive days, so 6 phases are covered. "Two mice were sacrificed by decapitation every 4 hours (ZT0, 4, 8, 12, 16 and 20)." (Methods). smFISH for Per1, Arntl, Cry1, Mlxipl, Actb. Data: GEO GSE182856.

## Method (if any)
ODE model of pre-mRNA, nuclear polyA mRNA and cytoplasmic mRNA with constant or rhythmic (cosine: phase, relative amplitude) production, export and degradation; four models M1-M4 chosen by likelihood plus BIC. Uses known period (24 h) and phase / amplitude relations between successive species; relative rates only (underdetermined). Outputs: export times, half-lives, model class per gene.

## Key results
- Subcellular fraction dominates the variance, time is minor: "readily split the types of samples" (PC1 nuclear vs cytoplasmic = 64 %); time effects sit on "smaller (PC4 and PC5, < 1.5% each)." (Fig 1C).
- 1,424 genes with rhythmic nuclear RNA, 79 also rhythmically degraded; 386 rhythmic only in cytoplasm (M3): "confidently classified 1,424 genes as having a rhythmic accumulation of nuclear RNA" (Fig 4D).
- About 15 % of expressed genes cycle in at least one compartment: "Thus, approximately 15% of expressed genes have a rhythmic temporal profile in at" (Results).
- Rhythmic export: 1,632 genes with rhythmic nuclear mRNA, 13 % exported rhythmically: "We confidently classified 1632 genes in one of the three models with rhythmic NAE" (Fig 5A); M3 peaks cluster near ZT8: "33% of the M3 genes peak in the time window between ZT6 and ZT10" (Fig 5C).
- Headline: "nuclear export contributes to the modulation and generation of rhythmic profiles of 10% of the cycling" nuclear mRNAs (Abstract). Median half-life 2.45 h; most export < 30 min.
- Core clock genes differ in nuclear fraction (Per1 more nuclear than Arntl, Cry1; p < 1e-4, smFISH; Fig 6).

## Not addressed / limits
- Mouse liver only, one tissue, no human data, no disease: "Although our observations were made only on mouse liver".
- Sex: males only; light-dark plus night-restricted feeding, so rhythms mix clock and feeding.
- n = 2 per time point; gene-level (isoforms summed), no alternative splicing; relative rates only.
- Nothing about tissue-tissue coordination or inter-individual phase variation.

## Link to us
- **Intro claim(s):** our current claim (refs_*.md): mouse liver nuclear vs cytoplasmic time series, post-transcriptional contributions (Naef lab). Supported as stated. Fits only as lab background (rhythms in mRNA are shaped by post-transcriptional steps, so pseudobulk clock-gene amplitudes and phases are not pure transcription readouts); weak support for any cancer / desynchrony statement.
- **Relation to clocksync:** background.

## Ideas
- Caveat for the method: clock-gene phase / amplitude read from total or cytoplasmic mRNA (bulk or sc) includes export and half-life effects (Per1 is more nuclear than Arntl and Cry1); snRNA vs scRNA (nuclear vs whole-cell) pseudobulks may therefore have systematically different offsets arg z; a check on paired sn / sc datasets is possible.
- Compare M2/M3/M4 gene lists with the clock gene panel used by clocksync to see whether any panel gene is post-transcriptionally rhythmic.
