# mure2018diurnal
**Ref:** Mure, Le, Benegiamo, Chang, Rios, Jillani, Ngotho, Kariuki, Dkhissi-Benyahya, Cooper, Panda, Science 359:eaao0318 (2018). DOI 10.1126/science.aao0318. Type: data. Lab: Panda (Salk) with Cooper (Inserm). NOT a Naef-lab paper.
**Source:** full text (main text and the one-page summary; Methods and supplementary tables not in the export)

## Question and design
Diurnal transcriptome of the baboon Papio anubis, a diurnal primate. 64 tissues including 22 brain regions, poly(A) RNA-seq (average 15.4 million
SE50 mapped reads per sample); 768 samples from young male baboons housed individually in seminatural conditions, 12 h light : 12 h dark.
**Sampling:** timed, full 24 h, "every 2 hours across the 24-hour day" (ZT0 to ZT22, 12 time points). 768 = 12 x 64, which implies one animal per time point
with all tissues taken from that animal (12 animals, average activity and temperature periods 23.98 and 23.97 h); to be confirmed in the supplementary tables.

## Method (if any)
MetaCycle (integration of several periodicity algorithms) for 24 h rhythms per tissue; phase set enrichment analysis for pathways; comparison with
published mouse tissue time courses (11 tissues; mice in constant darkness for 1 to 2 days).

## Key results
- Rhythmicity is pervasive: "with up to 81.7% of protein-coding genes showing daily rhythms in expression" (15,269 of 18,684 expressed PCGs); 16,442 genes (65.5 % of expressed) rhythmic in at least one tissue.
- Tissue specific: even in tissues with many cyclers, "less than 1% (<10) rhythmic genes were shared" between tissues (more than 1000 rhythmic genes); number cycling per tissue from below 200 (pineal, mesenteric lymph nodes, supraoptic nucleus, lateral hypothalamus, bone marrow) to above 3000 (thyroid, stomach fundus, gastrocnemius, PVN, prefrontal cortex).
- 97 % of 10,989 ubiquitously expressed genes cycle somewhere: rhythms add a tissue-specific layer to housekeeping functions.
- Global timing: peaks at ZT8 (11,146 genes) and ZT22 (5,127 genes), a quiescent zone at ZT14-19 (696 genes at ZT16).
- Clock genes: "the peak phase of expression of each clock component was mostly restricted" to about a 6 h window across tissues; e.g. Bmal1 ZT13.3 +/- 0.2 and Npas2 ZT14.3 +/- 0.5. Relative phases "components are remarkably conserved across" tissues; no clock gene cycles in all 64 tissues (Ciart in 52).
- Primate vs mouse: clock-gene phases are shifted "by ~12 hours between mouse and baboon tissues (with the notable exception of the SCN)"; the number of cyclers per tissue is uncorrelated (R2 = 0.086, P = 0.38).

## Not addressed / limits
- One animal per time point per tissue (as implied), males only, young, captive; no inter-individual variation in phase can be separated from time.
- Baboon, not human; no disease. Mouse comparison uses constant-darkness mouse data against LD baboon: "mice that were released into constant darkness" for 1 to 2 days.
- Dissection heterogeneity (e.g. brain nuclei, adjacent kidney cortex vs medulla with ZT2 vs ZT7 main phase clusters) makes tissue offsets partly compositional.
- Clock-gene phase spread across tissues is reported as a window (about 6 h), not decomposed into shared and tissue-specific offsets.

## Link to us
- **Intro claim(s):** refs_circadian_desynchrony.md claim: baboon atlas, 64 tissues, 2 h sampling, tissue-specific rhythmic programs. SUPPORTED. The claim "~82 %" is the 81.7 % of expressed PCGs. WRONG attribution in refs_circadian_desynchrony.md: "Naef lab, with Mure" is not right: Mure is at the Panda lab, Salk; no Naef author.
- **Relation to clocksync:** data we could use; background for tissue-specific offsets in a diurnal primate (human-like phase relative to light).

## Ideas
- Best available multi-tissue, known-time primate data: since each animal contributes all tissues (to confirm), it is a donor x compartment table with ground truth time; test clocksync offsets (arg z per tissue) against Mure's per-tissue clock-gene phases (about 6 h window) and the PVN, kidney cortex vs medulla offsets.
- Check whether clocksync on day-only subsets (ZT2-ZT10) still recovers tissue offsets and the SCN exception.
