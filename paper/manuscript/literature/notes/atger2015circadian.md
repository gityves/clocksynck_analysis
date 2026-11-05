# atger2015circadian
**Ref:** Atger, Gobet, Marquis, Martin, Wang, Weger, Lefebvre, Descombes, Naef, Gachon, PNAS 112:E6579 (2015). DOI 10.1073/pnas.1515308112. Type: biology + method. Lab: Naef (EPFL) with Gachon (Nestle).
**Source:** full text (main text with Methods; SI not in the export)

## Question and design
How do the circadian clock and feeding rhythms shape rhythmic transcription (pre-mRNA from intronic reads), mRNA accumulation (exonic reads) and
translation (ribosome footprints, RFP) in mouse liver? Male C57BL/6 mice 10 to 14 weeks old, individual livers, 84 samples for RNA-seq and Ribo-seq.
**Sampling:** timed, around the diurnal cycle under 12 h light : 12 h dark. "livers of individual mice every 2 h under" ad libitum feeding (ALF, WT); the same every 4 h in WT and Bmal1 KO under night-restricted feeding (RF).
All samples timed (mouse, zeitgeber time known); nothing untimed.

## Method (if any)
Harmonic regression y(t) = mu + alpha cos(2 pi t/T) + beta sin(2 pi t/T) on log2 signals (T = 24 h, and 12 h in ALF), with model selection across conditions: alpha and beta
set to zero (no rhythm) or shared by any subset of N conditions, number of models given by the Bell number B(N+1). Each model fitted by linear regression, BIC and Schwarz weights;
threshold "An arbitrary threshold of 0.4 is set on the BIC weight." (0.1 for the 877-model Bmal1 comparison). Intron/exon/RFP signals classify genes by level of rhythmicity.
This is the model-selection framework later packaged as dryR (Weger 2021).

## Key results
- Translation follows mRNA: for rhythmically translated mRNAs, "those (72%) exhibit rhythmic transcription (Fig. 2C)."; RFP shows no delay vs mRNA, about 6 h average delay to protein.
- Bmal1 KO affects transcription and post-transcription: "indicating that Bmal1 deletion affected the" general rhythmic pattern; phases shift "of both pre-mRNA and mRNA was shifted in the KO compared" with WT in both directions; "amplitudes of such rhythms at the mRNA level were, on average," decreased in KO (P = 1e-32) while pre-mRNA amplitude slightly increased (P = 1e-4).
- Genes moving from rhythmic mRNA only (WT) to rhythmic pre-mRNA and mRNA (KO) had 12 h pre-mRNA rhythms in WT that became 24 h: "described that 12-h-period genes changed to 24-h period in circadian clock-deficient mice".
- Translation efficiency rhythms: 5'-TOP genes (translation machinery, peak near ZT17) and TISU genes (mitochondria, near ZT10), mostly feeding-driven; sharper phase and larger amplitude under RF.
- Summary: "circadian clock and feeding rhythms have broad impact on rhythmic" gene expression.

## Not addressed / limits
- Liver only, male mice, single tissue; no cross-organ phase relations, no human data.
- Two sampling resolutions (2 h for ALF WT, 4 h for RF WT and KO), so the feeding comparison (ALF vs RF) is made at different resolution and replicate structure; the clock comparison (WT vs KO) is within RF only.
- Ribosome profiling gives translation, not protein activity; most rhythmic RFP encode non-rhythmic proteins.
- The method needs known sampling times and a full 24 h cycle; it classifies rhythms and cannot infer unknown time.

## Link to us
- **Intro claim(s):** refs_cancer_phase_inference.md: circadian vs feeding rhythms differentially control rhythmic transcription and translation, dryR precursor (Naef lab). SUPPORTED on both parts: the model-selection scheme here ("model selection approach. In the case of N conditions, we allowed α and β to" ...) is the dryR scheme. The statement that feeding contributes rhythms independent of the core clock is supported by the Bmal1 KO plus restricted feeding design.
- **Relation to clocksync:** background. Establishes that part of tissue rhythms is systemic (feeding), clock-independent, which bears on what shared host state clocksync measures when feeding times shift (shift work, hospital).

## Ideas
- A dryR-style comparison per compartment needs time; test whether clocksync's coupling c is lower where WT-vs-KO analyses show a larger systemic (feeding) fraction.
- Clock-independent systemic rhythms (Per1,2 persist in Bmal1 KO) could bias a clock-gene-based host state; check on Atger/Weger KO livers (GSE73554) whether clock-gene pseudobulks of KO mice are flagged as gate-failures.
