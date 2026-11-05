# masri2016lung
**Ref:** Masri S, Papagiannakopoulos T, Kinouchi K, Liu Y, Cervantes M, Baldi P, Jacks T, Sassone-Corsi P, Cell 2016;165:896-909. DOI 10.1016/j.cell.2016.04.039. Type: biology (mouse). Lab: Sassone-Corsi with Jacks. Data: GEO GSE73222.
**Source:** full text

## Question and design
KrasLSL-G12D;p53fl/fl (Cre, tumour-bearing TB) versus FlpO-treated p53fl/fl controls (WT), C57BL/6J, adenoviral intratracheal delivery. Host liver studied for the effect of a distant tumour; neither host clock nor tumour clock was perturbed experimentally and the tumour clock was not measured; hepatic clock versus tumour not compared. **Sampling:** timed and complete 24 h: livers every 4 h at ZT0, 4, 8, 12, 16, 20 under LD12:12 (microarray, metabolomics 600 metabolites, qPCR, western); locomotor free-running in DD for 2 weeks (n = 12 per group); indirect calorimetry n = 7 WT, 8 TB; per-ZT n and time after tumour initiation not recorded in the text read.

## Method (if any)
JTK_cycle rhythmicity (p < 0.01 genes, p < 0.05 metabolites), set logic (WT-only, TB-only, BOTH), DAVID GO, two-way ANOVA, cytokine array at ZT12, ChIP.

## Key results
- Rewiring: "peak in phase of expression was around ZT 8 in the WT category" while TB-only oscillating genes are bi-phasic at ZT0 and ZT12; the shared set "505 genes that retain oscillation in both WT and TB mice".
- "Strikingly, 46% of circadian genes exhibited a phase change, with 68% of these genes being phase advanced and 32% were phase delayed by at least 1 hr" (Fig 1E).
- Metabolome: "though the oscillation of 159 metabolites persisted, 90 were rhythmic exclusively in WT and 84 exclusively in TB mice"; 53 % of shared metabolites shift phase; lipids reduced; NAD+, ATP, acetyl-CoA lower.
- Core clock unchanged: "were unchanged in the livers of TB animals" (Clock, Bmal1, Per1-3, Cry1/2, Nr1d1, Dbp; qPCR Fig 3C), "no change in the free-running period was observed between WT and TB mice", food intake and RER still rhythmic.
- Authors conclude "the tumor does not appear to function as a classical zeitgeber" and call it an endogenous circadian organizer (ECO).
- Mediator: "Of the 31 cytokines assayed, a number were increased, decreased, or unchanged in TB mouse serum versus WT"; IL-6 significantly up; STAT3-Socs3 up (Socs3 peaks ZT8), insulin signalling (pAKT, IRS1) and AMPK/SREBP altered.
- Organ specificity: "this transcriptional rewiring observed in the liver is representative of a localized response, as the WAT and muscle gene expression profiles differ".

## Not addressed / limits
- Tumour (lung) clock is never measured; no phase or amplitude comparison between tumour and liver or between organs beyond liver vs WAT/muscle snapshots.
- Clock readout is per-gene qPCR of a few genes; an unchanged core clock is shown by overlay of means, with no formal phase or amplitude test (offset below about 1 h not excluded).
- Set-based comparison (WT-only, TB-only) is threshold dependent (JTK p < 0.01) and not a direct test of altered rhythms; the 46 % phase-change figure depends on a 1 h cutoff and 4 h sampling.
- One model (KP, Cre-induced), no metastasis; no human samples.

## Link to us
- **Intro claim(s):** the refs_*.md row (tumour-to-host clock effects) is overstated if read as host clock desynchronisation: the paper reports that the hepatic core clock and behaviour are unaffected, while clock-controlled metabolic transcripts and metabolites are rewired. It supports the reading that a tumour reorganises downstream circadian output of another organ without moving its core clock.
- **Relation to clocksync:** supports the premise that a host reference built on core clock genes is not shifted by a distant tumour; also data we could use (GSE73222, timed 6 ZT liver). Not yet used in paper/docs (searched METHOD.md, registered_experiments.md, results_preliminary.md).

## Ideas
- Register a positive-control test on GSE73222: core clock-gene phase and amplitude in TB vs WT liver with known ZT (expected offset about 0), plus the same readout on non-core rhythmic genes (expected shifts), to show the clock-gene panel is the robust host reference.
- Missing reference: Hojo 2017 (remote reprogramming of the hepatic circadian transcriptome by breast cancer), Hadadi 2020 (chronic jet lag in a mammary tumour model); Altman 2015 (MYC and clock).
