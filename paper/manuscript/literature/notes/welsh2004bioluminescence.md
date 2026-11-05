# welsh2004bioluminescence
**Ref:** Welsh, Yoo, Liu, Takahashi, Kay, Curr Biol 2004. DOI 10.1016/j.cub.2004.11.057. Type: biology. Lab: Kay / Takahashi.
**Source:** full text

## Question and design
Is the damping of peripheral circadian rhythms due to single-cell damping or to loss of synchrony? Bioluminescence
imaging of single cells: Rat-1 fibroblasts acutely transfected with Bmal1::luc (188 cells in two cultures, 5 and 24 days after
transfection; 144 rhythmic throughout) and primary tail fibroblasts from mPER2::LUC-SV40 knock-in mice (178 cells,
three cultures, 7-11 days). Luminometer records of whole dishes for comparison. **Sampling:** timed, continuous in vitro
recording of 1-2 weeks, medium change at the start as the only synchronising cue.

## Method (if any)
Sine fits to the first and last 3 days for start and end phases (circular statistics), chi-square periodogram (Q(P)),
amplitude factor per day; pairwise cell distance versus phase and period differences (6078 pairs).

## Key results
- Rat-1: mean period 23.14 +/- 1.98 h (range 19.1-30.3, n = 144). Start phases clustered (p < 0.0001, 125 deg), end phases random (p = 0.9). Amplitude factor 0.98/day (no damping) (Fig. 1).
- Primary fibroblasts: period 25.65 +/- 1.40 h (range 22.4-29.7, n = 178); start phases clustered (p < 0.0001), end phases random (p = 0.2); amplitude factor 1.03/day; Q(P) 296 versus 156 for Rat-1 (Figs 2, 3).
- Conclusion: "single fibroblasts can oscillate robustly and independently with undiminished amplitude and diverse circadian periods." (Summary)
- No functional coupling: no correlation between distance and start or end phase difference (6078 pairs; 444 within 500 um); distance versus period difference r = 0.09 (r2 = 0.01, p < 0.0001).
- Ensemble damping appears even while individual cells grow in amplitude, so "the rapid ensemble damping is much better explained by loss of synchrony among cells than by damping of individual cells."

## Not addressed / limits
- Fibroblasts only (Rat-1 and tail-derived primary cells); the authors write "we cannot rule out the possibility that other cell types may be damped oscillators or arrhythmic."
- Dissociated cultures: no coupling is found here, whereas SCN or tissue cells may couple; Izumo et al. challenge is
  addressed only qualitatively.
- Reporter plasmid dynamics (Rat-1) and low photon counts (primary cells, 2.2 photons/min/cell).
- Says nothing about in vivo phase relations between organs.

## Link to us
- **Intro claim(s):** cells are independent oscillators and population rhythm amplitude depends on synchrony.
  refs summary (independently phased, population damped by desynchrony) is supported, with the caveat above.
- **Relation to clocksync:** background. Explains why a bulk clock-gene signal in a compartment is a sum over cells
  of varying coupling: the same reading in our |z| (relative amplitude) can be low from dispersion, not from a weak clock.

## Ideas
- Single-cell clock-gene datasets (scRitmo) could test cell-to-cell phase dispersion directly: dispersion of cell phases within
  a donor-compartment versus its pseudobulk |z|.
- Same-lab companion for tissue: yoo2004periodluciferase (tissue explants persist; SCN lesion desynchronises organs).
