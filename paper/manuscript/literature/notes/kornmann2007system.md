# kornmann2007system
**Ref:** Kornmann, Schaad, Bujard, Takahashi, Schibler, PLoS Biol 2007 (5:e34). DOI 10.1371/journal.pbio.0050034. Type: biology (conditional hepatocyte clock). Lab: Schibler (Geneva).
**Source:** full text

## Question and design
Which rhythmic liver genes need the local hepatocyte oscillator and which are driven by systemic cues? Mouse: LAP-tTA x TRE-Rev-erba (double homozygous), hepatocyte-specific, Dox-dependent HA-REV-ERBa that represses Bmal1; without Dox the hepatocyte clock is arrested, with Dox it runs. Readouts: whole-liver Affymetrix arrays on 12 RNA pools (6 time points ZT0-ZT20, each of Dox and no-Dox mice, pools of 3 mice), "representing two entire days", LD12:12; Taqman and Western for Per2; mPer2::luc liver and lung explants in culture. **Sampling:** timed, full 24 h over two days; explants continuous.

## Method (if any)
Rhythmic transcripts called in Dox-fed mice (about 350), then compared in mice without Dox by amplitude and phase (permutation, 50,000 time-point permutations).

## Key results
- Perturbed: hepatocyte clock only; non-parenchymal cells (endothelial, bile duct, Kupffer) do not express the transgene and so keep their clock; the authors estimate "around 10% of the liver RNA is contributed by cells that do not express the LAP-tTA transgene" (rat-based volume estimate), which could contribute residual Bmal1.
- In vivo, hepatocyte clock off: the rhythmic transcription of most genes depended on functional hepatocyte clocks; fewer than 10% of the 351 circadian transcripts stayed rhythmic with high amplitude. "we discovered 31 genes, including the core clock gene mPer2" (41 feature sets) with unchanged amplitude and phase. Several heat-shock genes (Hspca, Hspa8, Hsp105) share the Per2 phase; Nocturnin, Cirbp too (Northern validated).
- In vitro, same livers: "circadian cycles of mPer2::luciferase bioluminescence could only be observed when hepatocyte oscillators were operational". Without Dox no rhythm (5 homozygous and 3 heterozygous mice); with Dox rhythms returned with a phase delay of approximately 6 h relative to mPer2::luc controls (time to decay of HA-REV-ERBa and restore Bmal1). Lung explants (transgene not expressed) rhythmic regardless of Dox.
- Reading: a liver without a working hepatocyte clock still shows rhythmic Per2 in the animal, driven by systemic cues (SCN-controlled: temperature, feeding, hormones), but no rhythm in isolation. Heat shock (150 min at about 40 C) induces mPer2 about 2-fold in liver explants without Dox.
- Authors' view of cells: peripheral oscillators are "self-sustained and cell autonomous", needing the SCN for phase coherence (abstract).

## Not addressed / limits
- Does not measure other liver cell types separately (whole liver pools of 3 mice; 10% non-transgene RNA); no sorted cells, no single cells, so cell-type pacing and cell-cell coupling are not tested.
- Cannot separate cell autonomy from coupling: explants without Dox lack rhythm, but this tests the hepatocyte oscillator, not whether hepatocytes couple. Duration of Dox withdrawal before sampling not stated in the passages read.
- One tissue (liver); in-animal data are under LD, so systemic cues are present throughout.

## Link to us
- **Intro claim(s):** liver rhythm in vivo mixes a local oscillator and systemic cues; ~10% of rhythmic genes (incl. Per2) follow the host, the rest need hepatocyte clocks. Systemic signals set the phase (Per2 gating) of peripheral clocks: the basis of a host reference.
- **Relation to clocksync:** background / supports. Direct support for "host-referenced" reading: some clock-gene rhythm in a tissue is imposed from outside and is shared across tissues. For the question: hepatocyte oscillators are autonomous (rhythm in explants) and produce most rhythmic output; coupling and dominance are not addressed. Bears on A315 in the sense that systemic rhythms alone give clock-gene (Per2) rhythm in a tissue with arrested clock, whereas brain death keeps tissue amplitude but loses the arrangement: the paper gives no time scale for rhythm persistence without systemic input in vivo.

## Ideas
- Clock-gene amplitude in tissues of brain-dead donors (A315) may partly reflect Per2/heat-shock/temperature-driven genes plus local oscillators; check which clock genes carry amplitude.
- Pair with sinturel2021circadian (hepatocyte-only clock keeps weak 24 h rhythm without SCN).
