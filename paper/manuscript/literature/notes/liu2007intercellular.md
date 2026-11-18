# liu2007intercellular
**Ref:** Liu, Welsh, Ko, Tran, Zhang, Priest, Buhr, Singer, Meeker, Verma, Doyle, Takahashi, Kay, Cell 2007 (129:605-616). DOI 10.1016/j.cell.2007.02.047. Type: biology (clock mutants, tissue vs single-cell). Lab: Kay (Scripps) / Takahashi (Northwestern).
**Source:** full text

## Question and design
Which clock genes are needed for sustained rhythms in cells, in peripheral tissues and in the SCN, and does intercellular coupling in the SCN compensate for cell-level defects? Per1, Per2, Per3, Cry1, Cry2 knockout mice crossed with the mPer2Luc reporter (germline knockout, whole body). Readouts: bioluminescence of SCN, lung, liver and cornea explants (recorded 2+ weeks, medium changes at day 8 and 15), dissociated primary and immortalised fibroblast populations, single fibroblasts (20 cells per genotype), single dissociated SCN neurons (at least 100 per genotype, 6-7 day imaging), plus a coupled-oscillator simulation. **Sampling:** ex vivo free-running.

## Method (if any)
FFT relative power and sine-fit goodness-of-fit to call a cell rhythmic; Rayleigh test of end-of-recording phases; spatial correlation of phase and distance; Leloup-Goldbeter model with noise, with and without coupling.

## Key results
- Peripheral tissue explants (lung, liver, cornea): wild type persistent; Cry1-/- and Per1-/- explants have much less persistent or absent rhythms (lung most), Cry1-/-Cry2-/- arrhythmic, Cry2-/- longer period. Gradual damping in wild type is "consistent with progressive cellular desynchrony due to period variation among cells".
- Cell autonomy in peripheral cells: dissociated fibroblasts, "Rhythmicity was detected in 80% of WT fibroblasts (n = 16/20)", 100% of Cry2-/- (20/20); only 5% of Per1-/- (1/20) and 30% of Cry1-/- (6/20) showed even weak or transient rhythms. So Per1 and Cry1 are needed cell-autonomously and the tissue phenotype equals the single-cell phenotype: the peripheral tissue is no more robust than its cells.
- SCN contrast: Per1-/- and Cry1-/- SCN slices stay rhythmic, but dissociated neurons are largely arrhythmic (rhythmic: 10% of Per1-/- n = 12/125; 6% of Cry1-/- n = 7/118; versus 64% of WT, n = 106/165, and 84% of Cry2-/-). Coupled slices synchronise and stabilise weak cells; "unlike in peripheral tissues, the SCN circadian clock is more than just the sum of its cells".
- Dissociated SCN neurons "drifted out of phase from one another by the end of the recording" (Rayleigh P = 0.97, n = 121, dissociated; P < 0.0001, n = 40, slice); no local coupling in dispersal: "closer cells did not tend to have more similar phases" (r = 0.04).
- Simulation: a simple Per-mRNA-dependent coupling restores persistent rhythms to noisy arrhythmic cells.

## Not addressed / limits
- Perturbation is germline and in all cells; it does not remove one cell type inside a tissue, so says nothing on pacing by a dominant cell population in peripheral tissue. Coupling is shown for the SCN only; for peripheral tissue the conclusion that coupling does not rescue Per1/Cry1 loss is an inference from the explant versus single-cell comparison (no uncoupling experiment inside a peripheral tissue).
- Explants lose cell composition and in vivo input; dissociated fibroblasts are not in a tissue; reporter is Per2Luc only.
- Time scale: days to weeks in vitro; explant rhythms damp gradually over this period.

## Link to us
- **Intro claim(s):** peripheral cells are cell-autonomous oscillators whose tissue rhythm is the sum of cells that desynchronise slowly in vitro, without strong coupling (the classic view that Finger 2021 revisits); SCN differs by intercellular coupling.
- **Relation to clocksync:** background. Supports cell autonomy for peripheral cells; no coupling or dominant population shown. Relevance to A315: explant tissue rhythms persist for days to weeks without any host signal and damp by cell dephasing, giving a time scale (days) much longer than ventilation, so brain-dead donors retaining clock amplitude is expected from autonomous oscillators.

## Ideas
- Contrast with noguchi2013fibroblast and finger2021intercellular (density dependence, paracrine factors) which show that fibroblast cultures are not purely independent.
