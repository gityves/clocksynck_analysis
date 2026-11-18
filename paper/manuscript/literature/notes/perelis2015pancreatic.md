# perelis2015pancreatic
**Ref:** Perelis, Marcheva et al., Science 2015 (350:aac4250). DOI 10.1126/science.aac4250. Type: biology (islet genomics, enhancers; cell-type-specific knockouts). Lab: Bass (Northwestern).
**Source:** full text

## Question and design
Which genes does the cell-autonomous clock drive in the pancreatic beta cell, and how? Isolated mouse islets (wild type, forskolin-synchronised ex vivo) with RNA-seq over two circadian cycles, CLOCK/BMAL1/PDX1 and histone-mark ChIP-seq, comparison with liver ChIP data, human islets (RNA-seq every 4 h for 48 h, groups of 200 islets), the pancreas-specific knockout PdxCre;Bmal1flx/flx (islet RNA-seq at one time, ZT2) and an inducible knockout PdxCreER;Bmal1flx/flx (tamoxifen in vivo in adults, or ex vivo in islets). **Sampling:** WT islets ex vivo every 4 h over 72 h after forskolin shock (insulin secretion) and RNA-seq over two cycles; knockout islets at a single ZT (ZT2); in vivo metabolic tests at ZT2 and ZT14.

## Method (if any)
eJTK_CYCLE for cycling transcripts; differential expression (FDR < 0.05) between genotypes; ChIP-seq peak classification (promoter proximal versus distal, within 2 kb of the TSS).

## Key results
- **Rhythms in a mixed islet, targeted cells inferred.** "We detected a total of 3905 cycling transcripts (Bonferroni-corrected P < 0.05), which accounted for ~27% of all expressed transcripts within the islet". The text calls this the beta-cell clock; islets are read as a unit (alpha and other endocrine cells, endothelium are not separated in the text read).
- **Driver and targeting.** PdxCreER with tamoxifen "abrogates BMAL1 expression exclusively within the b cell (fig. S6)" in vivo; ex vivo tamoxifen on islets gives ">60% reduction in Bmal1 expression". Specificity and efficiency are shown in supplementary figures not in the text read.
- **Targeted cells' output.** Ex vivo, vehicle-treated PdxCreER;Bmal1flx/flx islets have higher GSIS at the zenith (48 h after shock) than at the nadir (36 h) (P = 0.038), tamoxifen-treated islets do not (P = 0.974): a rhythm of secretion that needs islet Bmal1 and persists without the body.
- **Knockout transcriptome, single time.** In PdxCre;Bmal1flx/flx islets at ZT2, 1757 genes change (1074 down, 683 up); "a total of 720 oscillating genes exhibited altered expression in animals with disrupted pancreatic clock function". Clock genes Per2, Rev-erb-alpha (Nr1d1), Tef, Nfil3 are among the strongest changes. These are single-time comparisons, not rhythm tests.
- **Other tissues and behaviour.** The inducible knockout "displayed normal wheel-running rhythms, period length, food intake, and body weight (fig. S7)"; hyperglycaemia and hypoinsulinaemia within 10-14 days after tamoxifen at both ZT2 and ZT14, "These results establish that circadian disruption in fully differentiated cells is sufficient to induce metabolic disease". BMAL1 binding sites at cycling genes differ between beta cells and liver; PDX1 co-binds active enhancers.
- Human islets: about 1800 cycling RNAs (Bonferroni P < 0.05), 481 orthologous to cycling mouse islet genes.

## Not addressed / limits
- Whole-islet RNA-seq: non-beta cells' clocks, and the tissue-level rhythm after the knockout, are not measured (knockout islets at one ZT, no time series).
- Ex vivo islets are removed from systemic cues; whether in vivo islet rhythms are autonomous, coupled or driven is not tested.
- No other tissue's clock reported in the knockout (activity and food intake normal only).

## Link to us
- **Intro claim(s):** clock-driven output of one cell type can be cell-autonomous (isolated islets rhythmic, rhythm needs Bmal1); mouse and human islets share cycling genes. For within-tissue organisation: **not informative about coupling or pacing** (islets mixed, knockout without time series).
- **Relation to clocksync:** data we could use. GEO GSE69889 (RNA-seq including "lists of cycling genes in mouse and human islets") and GSE70960 (ChIP-seq); the paper does not say whether the knockout ZT2 islets are in GSE69889. The human islet RNA-seq (every 4 h, 48 h) is a timed non-cancer time series of one pancreatic compartment.

## Ideas
- Pull GSE69889 and check whether the PdxCre knockout versus control islet RNA-seq is in it; with its single time point, only mean shifts in clock genes (e.g. whether Per2 or Nr1d1 drop to the residual level expected for non-beta islet cells) can be read.
- Use the human islet cycling list as an islet-specific reference in human atlases containing pancreas.
