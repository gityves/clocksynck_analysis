# yeung2018transcription
**Ref:** Yeung, Mermet, Jouffe, Marquis, Charpagne, Gachon, Naef, Genome Res 28:182 (2018). DOI 10.1101/gr.222430.117. Type: biology + method. Lab: Naef (EPFL).
**Source:** full text (main text; Supplemental Methods not in the export)

## Question and design
How do tissue-specific rhythms arise on top of a shared clock? Mouse liver and kidney (RNA-seq in Bmal1 KO and WT littermates), 11 WT tissues from
Zhang 2014, MARA motif activity, DNase-seq, ChIP-exo, 4C-seq (liver vs kidney at the time of peak mRNA of Mreg, Pik3ap1, Slc44a1, Slc45a3).
**Sampling:** timed. Zhang set: "was obtained under dark–dark (DD), ad libitum feeding, sampled every 2 h." (24 time points over 2 days, 11 tissues). New liver (GSE73554, from Atger 2015) and kidney
Bmal1 KO/WT: "obtained under light–dark (LD), night-restricted feeding conditions, sampled every 4 h." Replicates n = 2 for TF activities and 4C.

## Method (if any)
Model selection (MS) extending harmonic regression (Fisher 1929) to assign each gene to a module of rhythmic profiles shared across subsets of
tissues (tissue-wide, tissue-restricted, tissue-specific), BIC-based. Modules summarised by complex-valued SVD: genes by log2 peak-to-trough amplitude and phase, and tissues by a
vector with relative amplitude and phase offset: "tissue vectors are scaled such that the highest amplitude tissue has length of 1 and a phase offset of 0 h". Clock- vs system-driven rhythms by WT vs Bmal1 KO.

## Key results
- Variance: "Most of the variance concerned differences in expression between tissues"; 24 h is the largest temporal component, 12 h close to background.
- Clock required for high amplitude: "mRNAs was reduced by threefold compared to WT littermates" in Bmal1 KO liver and kidney, more so at large amplitudes; "rhythms across diverse tissues, while systemic signals regulate lower amplitude".
- Shared clock module: tissue-wide transcripts (clock genes, Dbp, Nfil3) "oscillated in synchrony across all tissues" at fixed times of day, with amplitudes differing (brain smallest).
- Tissue-restricted modules with fixed offsets: "offsets, such as the nearly 12-h shifted rhythms in brown fat and" skeletal muscle (n = 20 genes, first SVD component 81 % of variance).
- Per1/Per2: "clock genes Per1,2 continued to oscillate in Bmal1 KO in multiple" tissues (clock-independent, systemic input).
- Liver: system-driven module n = 1395 (84 % variance, feeding-linked), clock-driven liver-specific n = 991 (83 %), kidney-specific clock-driven n = 156 (80 %). TF activities (MAFB, EGR1, TFCP2, ELF, RORE peaking near ZT21, E-box near ZT7) explain tissue-specific rhythms; liver-specific chromatin loops recruit clock-bound enhancers, promoter-specific at below 10 kb.

## Not addressed / limits
- Cohort-level time courses per tissue; tissue vectors are properties of gene modules (not of an individual animal), so no inter-individual variation in tissue offsets.
- Zhang (DD, ad libitum) and the Bmal1 KO sets (LD, night-restricted feeding) are analysed separately because designs differ.
- Only mouse; liver and kidney for loops; 4C-seq at a few loci; no disease.
- Replicate numbers are small (n = 2 for TF and loops).

## Link to us
- **Intro claim(s):** refs_circadian_desynchrony.md: TF activity rhythms and tissue-specific chromatin loops explain tissue-specific circadian expression (liver/kidney). SUPPORTED. It also supports (not claimed so far) a shared clock core with tissue-specific amplitudes and offsets.
- **Relation to clocksync:** background and closest conceptual ancestor in the lab: rank-one complex SVD of genes x tissues with per-tissue amplitude and phase relative to a reference tissue is formally close to clocksync's per-compartment z (|z|, arg z) relative to a shared host state. Not the same: here time is known and genes are modules, there time is unknown and the shared state is learned from the donor's other compartments.

## Ideas
- Cite as the timed, genes-x-tissues counterpart of clocksync; compare tissue vectors of the tissue-wide module (11 tissues, Zhang GSE54652) with clocksync's |z|, arg z estimated without time: agreement would be a ground-truth validation.
- Brown fat vs muscle antiphase (about 12 h) is a testable positive control for arg z in untimed data.
- Reference that cites it: Sobel 2017 (liver H3K27ac) and Korenčič 2014 (tissue-specific rhythms); neither is in our batch.
