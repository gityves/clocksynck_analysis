# preitner2002orphan
**Ref:** Preitner et al., Cell 2002 (110:251-260). Type: biology. Lab: Schibler (Geneva).
**Source:** full text

## Question and design
Mouse liver and SCN; EMSA of liver nuclear extracts, RNase protection/Northern/TaqMan, Western blots, in situ in SCN, Rev-erb-alpha
knockout mouse (null allele, F2 progeny n = 383: 23.5% wt, 48.6% het, 27.9% hom), wheel-running behaviour. **Sampling:** timed. Liver
and SCN sampled every 4 h around the clock (ZT in LD; CT in DD for protein and SCN); four animals per time for qPCR; behaviour in DD
and LL. Full 24 h, single tissue per assay (liver, SCN).

## Method (if any)
None.

## Key results
- The Bmal1 promoter has two ROREs; EMSA shows a strongly circadian RORE binding complex with a maximum at ZT 8, matched to the REV-ERB-alpha profile (Fig. 1).
- REV-ERB-alpha and Bmal1 pre-mRNA are antiphasic (Bmal1 pre-mRNA maximum at ZT 20). Abstract: "Bmal1 and Per transcription cycles display nearly opposite phases and are thus governed by different mechanisms".
- Knockout abolishes the Bmal1 rhythm: "the amplitudes of Bmal1 mRNA and premRNA oscillations are at least 20-fold in wild-type animals, but less than 2-fold" in the mutant (Fig. 3).
- Clock mRNA amplitude is about 2.5-fold in wild type and nearly flat in the knockout; Cry1 about 7-fold vs about 2.5-fold; Per2 and Cry2 barely change (Fig. 4).
- Behaviour persists: "does not result in arrhythmic behavior when mice are placed in a constant environment"; period shorter in DD (23.76 vs 23.38 h) and LL (24.59 vs 24.01 h), p < 0.002.
- "light pulses delivered at CT 22 provoked dramatic phase advances" (5.4 +/- 1.8 h vs 1.1 +/- 0.5 h in wild type; Fig. 5).
- Rev-erb-alpha is controlled by the negative limb: abstract "it is not required for circadian rhythm generation".

## Not addressed / limits
- Wording matters: REV-ERB-alpha is shown to be a major but not essential regulator; the paper itself titles it the positive limb, and "stabilising loop" is our phrase, not theirs.
- Mouse, liver and SCN only; no other peripheral organ is examined, and no human data.
- Single genotype comparison at a mixed 129/Sv x C57BL/6J background; no double REV-ERB-alpha/beta knockout (beta is expressed in liver).

## Link to us
- **Intro claim(s):** REV-ERB represses Bmal1 and sets the antiphase between Bmal1 and Per/Cry (Circadian cycle); also justifies using Bmal1 and Per/Cry/Dbp as opposite-phase readouts.
- **Our current claim (refs_circadian_desynchrony.md):** 'REV-ERBalpha represses Bmal1 transcription, the stabilising loop of the mouse clock'. SUPPORTED for repression of Bmal1, but "stabilising loop" is an interpretation; the paper says the loop is dispensable for rhythm generation and modulates period and phase shifting.
- **Relation to clocksync:** background. The Bmal1-Per antiphase and the reduced Clock/Cry1 rhythm in the knockout are the biological reason a fixed gene-wise phase offset is plausible; mutants with a damped Bmal1 limb are a test case for the |z| (relative amplitude) readout.

## Ideas
- Check in GTEx/single-cell data whether NR1D1/NR1D2 pseudobulk covariance with BMAL1 (ARNTL) is antiphase and whether the gate flags compartments with a weakened Bmal1 limb.
- Missing: Sato et al. 2004 (ROR/REV-ERB loop, stabilising role, Neuron) and Cho et al. 2012 (REV-ERB alpha/beta double knockout, Nature).
