# guan2020hepatocyte
**Ref:** Guan et al., Science 2020 (369:1388-1394). DOI 10.1126/science.aba8984. Type: biology (cell-type-specific clock knockout, liver). Lab: Lazar (UPenn).
**Source:** full text

## Question and design
Adult mouse liver. REV-ERBa/b deleted in hepatocytes only (HepDKO; AAV8-TBG-Cre into floxed mice, adult, so no developmental effect). Question: is the hepatocyte clock needed for rhythms of other liver cell types, and how does it rank against feeding? Readouts: bulk liver RNA-seq every 3 h (8 time points, n = 3 per time point); single-nucleus RNA-seq at ZT8 (18,239 nuclei; hepatocytes, endothelial cells ECs, Kupffer cells KCs, stellate, immune); sorted ECs and KCs, RNA-seq every 6 h (n = 4 per time point) and metabolomics (n = 3-4 per time point); reverse-phase feeding (RPF, food only in the light phase) for 3 weeks. **Sampling:** timed, LD, full 24 h; knockout analysed 2 weeks after AAV.

## Method (if any)
JTK_CYCLE rhythmicity (adj. P < 0.01, 21-24 h, peak-to-trough > 2 for liver; adj. P <= 0.05, > 1.5 for sorted cells); TF-motif and enhancer-RNA analysis (IMAGE); NicheNet ligand-receptor inference hepatocyte to EC/KC. Rhythmic genes classed as HepDKO-dominant, RPF-dominant, both, or retained.

## Key results
- Perturbed: hepatocyte clock only. Hepatocyte core clock broken: "both of these genes were constitutively derepressed in the REV-ERB HepDKO" (Bmal1, Npas2); other core clock genes lose rhythmicity. Behaviour, feeding and plasma insulin rhythms not much affected (fig. S3), so the effect is not behavioural.
- Other cells keep their own core clock but their transcriptome is remodelled. "Rev-erba/b gene expression was virtually unchanged in the ECs and KCs from the HepDKO livers, although the amplitude of Rev-erba/b rhythms was muted in KCs" (Fig. 2C). "the diurnal rhythmic transcriptomes were extensively remodeled in ECs" (Fig. 2D): rhythms lost and gained, widespread enhancer reprogramming. KCs likewise. Largest number of changed transcripts in ECs, then KCs (Fig. 2B); cell composition unchanged.
- Cell-type specific: "little overlap between different cell types" of remodelled transcripts (fig. S6). Candidate mediators: ligands (e.g. Csf1 loses rhythm in HepDKO hepatocytes, receptor Csf1r in KCs, target Cxcl10 loses rhythm) and metabolites; correlative, not tested.
- Hepatocyte clock does not simply set the phase of other cells (no phase-offset numbers given); the effect is on which genes are rhythmic.
- Feeding hierarchy: RPF shifts core clock genes of control livers by 12 h, "nearly all rhythmic transcripts exhibited a 12-hour phase shift in the livers of control mice under RPF". Of rhythmic transcripts 11.5% HepDKO-dominant, 30.7% RPF-dominant, 45.0% both, 12.8% retained in HepDKO plus RPF (rhythm independent of the hepatocyte clock and feeding; phase follows food). In ECs about 74% and in KCs about 65% of cell-specific rhythmic genes are regulated by both HepDKO and RPF. "Thus, nonautonomous signals resulting from feeding and communication from hepatocytes play vital roles" in non-hepatocyte rhythms.
- Authors' summary: "We demonstrate that some but not all hepatocyte diurnal rhythms are controlled by the core clock in a cell-autonomous manner in vivo."

## Not addressed / limits
- Never measures whether EC or KC clock phase moves relative to hepatocytes; reports rhythmic gene sets, not clock-gene phase. Core clock of ECs/KCs "virtually unchanged" means these cells keep a clock and presumably receive the same systemic cues.
- Mouse, LD, AAV knockout (2 weeks); sorted-cell RNA-seq at 6 h resolution, n = 4; sNuc-seq only at ZT8. No time course beyond 2 weeks (3 weeks for RPF), no SCN or extrahepatic manipulation.
- Mediators inferred computationally (NicheNet); no single-cell dispersion or cell-to-cell coupling measurement.

## Link to us
- **Intro claim(s):** hepatocytes influence rhythms of neighbouring cell types (non-parenchymal liver cells); peripheral rhythm is both cell-intrinsic and imposed (feeding, hepatocyte signals). Supports intercellular dependence of output rhythms, and a hepatocyte lead on output programs, not on the clock phase of other cells.
- **Relation to clocksync:** background (supports intercellular influence with a dominant population, hepatocytes, on non-hepatocyte rhythmic outputs; core clock of ECs/KCs stays). For the question: autonomy of the core clock in each cell type (kept), but output rhythms of ECs/KCs depend on hepatocytes and on feeding. Does not bear directly on the brain-death result; its time scale (2 to 3 weeks) is far longer than ventilation; the RPF result shows feeding alone moves all liver cell rhythms (12 h), consistent with a shared tissue phase.

## Ideas
- Prediction for our sc atlases: if hepatocytes pace non-hepatocyte output but not clock phase, class phases stay within about 1 h while non-hepatocyte rhythmic output genes change; test only in liver sc data with timed samples.
- Cited but missing in our references: Koronowski 2019, Saini 2013 (see sinturel2021circadian card).
