# finger2021intercellular
**Ref:** Finger, Jaeschke, del Olmo, Hurwitz, Granada, Herzel, Kramer, Sci Adv 2021 (7:eabg5174). DOI 10.1126/sciadv.abg5174. Type: biology (coupling, in vitro). Lab: Kramer (Charite, Berlin).
**Source:** full text

## Question and design
Do peripheral circadian oscillators couple intercellularly, and by what mediator? Mainly human U-2 OS cells (Bmal1:Luc, Per2:Luc, NR1D1::VNP reporters), 2D and 3D spheroids, plus mouse lung explants, PER2::LUC organotypic slices, liver organoids, primary hepatocytes (as conditioned-medium source). Readouts: population bioluminescence, single-cell fluorescence imaging over 5 days, RNA-seq, mass spectrometry of active fractions. **Sampling:** in vitro, continuous recording of free-running rhythms after synchronisation; no in vivo data.

## Method (if any)
Coculture of reporter cells with phase- or period-different non-reporter cells (phase and period pulling); culture-density series; inhibition of secretion; conditioned medium (CM) stimulation; two-step chromatography plus mass spectrometry; TGF-beta antibody depletion, recombinant TGF-beta, receptor blocker LY2109761, RNAi screen; stochastic amplitude-phase oscillator model (200 oscillators). Single-cell phase coherence by pyBOAT.

## Key results
- Phase pulling in 2D: reporter cells cocultured with 6 h advanced or delayed non-reporter cells moved toward them, but "the extent of phase adjustment was at most 2 hours (despite a 6-hour phase difference)": weak coupling. Arrhythmic (Bmal1 knockdown) partners attenuate the pulling.
- Period pulling only in 3D spheroids: "we did not observe period changes of peripheral clock cells under conventional 2D culture"; in spheroids with a fivefold excess of WT cells, the long period of CRY2-/- reporter cells shortened and the short period of TNPO1-/- lengthened.
- Density: "Reducing culture density led to amplitude reduction, increased damping, and period lengthening"; single-cell phase coherence is lower at low density although cell division is also lower (so not division noise); rescue by coculture with non-reporter cells, also when physically separated (diffusible factors).
- Secretory pathway: inhibitors and RNAi (GBF1, SEC13, TMED2/10) in U-2 OS and mouse lung explants cause amplitude reduction, increased damping and period lengthening.
- Mediator: CM from human and mouse cells, including primary hepatocytes, gives "robust phase delays of 4 to 8 hours upon CM treatment" at the PER2 trough; mouse explants and liver organoids also delayed by U-2 OS CM (conserved across tissues and species). Mechanism: CRE-driven immediate-early induction of PER2 (2-4 h), then suppression of DBP and NR1D1. Mass spectrometry: 431 proteins in active fractions, 37 secreted; TGF-beta2 in all active fractions.
- TGF-beta as coupling factor: immunodepletion abolishes CM activity; recombinant TGF-beta induces CRE and PER2; LY2109761 abolishes CM phase shifts. Inhibition raises the speed of single-cell desynchronisation (NR1D1::VNP; solvent 215 cells, inhibitor 180 cells; linear fit over days 1-3, R2 0.95 and 0.96), reduces amplitude and increases damping of U-2 OS and of PER2::LUC lung slices. Robustness: a 20 C cold pulse of 8 h gave no significant shift with solvent but shifted inhibitor-treated cells by about 8 h.

## Not addressed / limits
- Cell line (U-2 OS, osteosarcoma) and immortalised or primary cultures; coupling shown as weak (2 h of 6 h); no cell-type-specific coupling, no in vivo test, whether TGF-beta signalling is itself rhythmic is unknown ("currently unknown").
- Between cell types: CM from hepatocytes, fibroblast lines etc. acts on U-2 OS, so the factor is not cell-type specific; no experiment asks which cell type dominates.
- Time scale: days of free-running recording.

## Link to us
- **Intro claim(s):** peripheral cells are coupled by paracrine factors (TGF-beta), so tissue coherence is not purely cell-autonomous; coupling is weak and raises amplitude, lowers damping and sets sensitivity to zeitgebers.
- **Relation to clocksync:** background. Supports intercellular coupling and its observable consequences (amplitude, coherence). Mediator conserved across cell types and species, consistent with a shared cell-class phase (within about 1 h) and no strict leader. TGF-beta is secreted by fibroblasts and mural cells in vivo, a possible link to why stromal-vascular cells share most clock state in our atlases (hypothesis only).

## Ideas
- Test whether TGF-beta pathway scores in stromal cells correlate with their shared-clock-state fraction across tissues (post hoc, registered first).
- Companion: finger2021peripheral (commentary) and sinturel2021circadian.
