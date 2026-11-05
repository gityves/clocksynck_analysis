# hollis2025reconstructed
**Ref:** Hollis, Sharma, Sheehan, Maggi, Weber, Hammarlund, Bennett, Menon, Musiek, Anafi, Neuron 113:2822-2838 (2025). DOI 10.1016/j.neuron.2025.07.010. Type: method application + biology (human data, APP/PS1 mice, model). Labs: Anafi (Penn), Musiek (WashU).
**Source:** full text.

## Question and design
Which molecular rhythms of the human cortex change in Alzheimer's disease, by cell type. ROSMAP DLPFC snRNA-seq, 409 subjects (159 AD dementia, 250 aged controls), pseudobulk per cell type and subject; mouse APP/PS1 microglia validation (polysome profiling, western) and a clock-translation model.
**Sampling:** no clock time. "Studying inaccessible tissues from human subjects, researchers often use post-mortem samples and time of death (TOD) to reconstruct rhythms", but ROSMAP lacks TOD, so one phase per subject is inferred (cohort-level ordering of 409 subjects around the circle).

## Method (if any)
- CYCLOPS 2.0 (autoencoder, covariates sequencing batch and disease state) applied to pseudobulk of excitatory neuron subclusters 3-5, chosen with clock correlation distance (CCD) as the compartment with the clearest clock structure; CHIRAL gives a similar order (circular rank correlation 0.84).
- Unit of inference: subject (one phase per donor), then gene-wise cosinor and nested-model differential rhythm tests per cell type, adjusting PMI, sex, batch.
- Needs full-cycle coverage across the cohort (subjects spread over the day), fixed gene acrophases from a curated gene list, and one shared phase per donor. Other cell types are placed on the neuron-derived phase.

## Key results
- Phases: CYCLOPS 2.0 "predicted TOD with consistent accuracy across age groups" in bulk human cortex data; acrophases agree with earlier bulk BA11 data.
- Cycling transcripts in control cells: 5,001 excitatory, 2,209 inhibitory, 885 astrocyte, 504 microglia (BH.q<0.1, relative amplitude >= 0.2).
- Core clock genes preserved in AD: "Comparing core clock rhythms in CTL and AD across cell types, we found almost no differences" (only TEF amplitude up in microglia).
- Output rhythms disrupted: ribosome and OXPHOS genes lose amplitude and mean level in AD across cell types; stronger link to CERAD plaque load than to Braak.
- Cross-cell-type phase: "Rhythms remain synchronized across excitatory neurons (green), inhibitory (orange), astrocytes (blue), and microglia (red) in AD", a qualitative reading from acrophase plots, with no per-donor offset statistics.
- Model: lower translation rate and more translational noise raise cycle-to-cycle variability of oscillator period; mouse APP/PS1 circadian jet-lag reduces ribosomal protein.

## Not addressed / limits
- Cell types are not compared within donor: all use the neuron-derived phase, so a donor-level neuron-vs-glia offset cannot appear by construction; the authors argue a fixed offset "would have detected this pattern", but shown only as pooled acrophases.
- No disease-specificity control (no other disease, no non-AD cohort); condition is a covariate in the ordering, so AD-specific phase structure can be absorbed.
- Informatic ordering is "a central weakness" per the authors; poorly ordered diseased tissue could lead to false conclusions. Subject-level phase scatter is not reported as an outcome.
- Phases are relative and cannot be tied to external time; diurnal (behaviour) and circadian rhythms are not separated.
- Glia underpowered (fewer nuclei) so null glial effects are weak evidence.

## Link to us
- **Intro claim(s):** refs_cancer_phase_inference.md describes it as cohort-level CYCLOPS reconstruction of cell-type rhythms (bulk + deconvolution). The bulk plus deconvolution part is wrong: the data are snRNA-seq pseudobulk per cell type, not deconvolved bulk. Otherwise supported; also supports that AD preserves the core clock but disrupts outputs and the need for phase-ordering.
- **Relation to clocksync:** phase-ordering competitor on the same data (ROSMAP/RADC); it needs a cohort that spans the day, here via ordering, and assumes cell types share one clock phase. We estimate the neuron/glia offset and coupling per donor.

## Ideas
- Done: clocksync RADC (AD_radc, 26 controls) is in A135 and A141; AD contrasts in RADC are noise-driven, not estimable (results_preliminary.md, A141).
- Test the paper's claim that core clock phases are shared across cell types in AD: compare arg z between excitatory neurons and glia in RADC AD vs control donors at common spread d* (needs power check first; A141 says low).
- Use their CYCLOPS phases (if released) as an external, partly timed truth for the donor-level host state in RADC.
