# koike2012transcriptional
**Ref:** Koike et al., Science 2012 (338:349-354). Type: data / biology. Lab: Takahashi (UT Southwestern).
**Source:** full text (main text; supplement not included).

## Question and design
Mouse liver, constant darkness. About 200 ChIP-seq samples (more than 9 billion reads) for BMAL1, CLOCK, NPAS2, PER1, PER2,
CRY1, CRY2, RNAPII (8WG16, Ser5P), p300, CBP and histone marks; whole-transcriptome RNA-seq (exon and intron). **Sampling:**
timed, free-running: ChIP-seq at six time points (CT0 to CT20, every 4 h); RNA-seq every 4 h over 48 h in DD. Full cycle, one tissue
(liver), mouse only; knockout mice as antibody controls.

## Method (if any)
Peak calling (MACS, PeakSplitter) to a master peak list per factor; cyclic binding and cycling transcripts by ARSER (P < 0.05);
intron signal read as nascent transcription, exon signal as mRNA. Cohort-level phases, not per-individual.

## Key results
- Sites: "We find ~5900 BMAL1 sites, ~4600 CLOCK sites, and ~2300 NPAS2 sites for the activators" and ~4600 PER1, ~7300 PER2, ~16,000 CRY1, ~10,000 CRY2 sites; 1444 bound by all six factors (Fig. 1D).
- Binding phases: activators "peak at CT6 to 8"; "PER1, PER2, CRY1, and CRY2 show binding rhythms that peak at CT15.9, CT17.3, CT0.4, and CT15.4, respectively" (Fig. 1C).
- Rhythmic transcripts: "we detected 1371 intron and 2037 exon RNA cycling transcripts" (Fig. 2); intron cycling peaks cluster at CT15.1 whereas exon cycling genes spread over all phases.
- Only 458 genes (22% of exon cycling) overlap: abstract "Only 22% of messenger RNA (mRNA) cycling genes are driven by de novo transcription", so post-transcriptional control shapes most mRNA rhythms.
- Three phases of the cycle: poised (CRY1 with CLOCK:BMAL1, CT0), de novo transcription (peak CT15), repressed (CT16-20); RNAPII 8WG16 peaks at CT14.5, Ser5P at CT0.6.
- Targets of all seven factors are enriched in metabolic pathways and "pathways in cancer" (tables S3/S4).

## Not addressed / limits
- Liver only; the transcriptional cycle of other organs or cell types, and of humans, is not measured (the discussion states the architecture for murine liver).
- Mouse housed in DD with six time points per factor: no inter-individual variation in phase, so no statement on coupling between organs or donors.
- Mean binding and mRNA phases are population averages from pooled samples; no per-gene uncertainty on acrophase.

## Link to us
- **Intro claim(s):** genome-wide architecture of the core loop and phases of its factors (Circadian cycle); also that mRNA rhythms in a tissue are only a fraction of transcriptional rhythms (Technical challenges: why a small fixed set of clock genes, not all cycling genes, is the readout).
- **Our current claim (refs_circadian_desynchrony.md):** 'Genome-wide ChIP-seq of core clock factors in mouse liver; defines the transcriptional architecture and phase-specific binding'. SUPPORTED. Add that it also covers RNAPII/chromatin and the 22% result.
- **Relation to clocksync:** background / data we could use (liver acrophases of Dbp, Per, Nr1d1 as priors).

## Ideas
- Use the Fig. 1C phase table (CT of binding) to define an expected activator-repressor offset (~9-10 h CLOCK:BMAL1 vs PER/CRY2) and compare with the offset the model recovers between Bmal1-type and Per/Cry-type genes.
- Missing: Zhang et al. 2014 PNAS (mouse circadian atlas of 12 organs) for tissue-resolved phases if not already cited.
