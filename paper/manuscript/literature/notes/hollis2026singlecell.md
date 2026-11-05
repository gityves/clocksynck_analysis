# hollis2026singlecell
**Ref:** Hollis, Veltri, Korac, Menon, Bennett, Ronnekleiv-Kelly, Kim, Anafi, bioRxiv 2026.03.23.713759 (posted 25 Mar 2026; not peer reviewed). Type: method (ORPHEUS) + biology. Lab: Ron Anafi (Penn).
**Source:** full text (preprint).

## Question and design
Bulk rhythm amplitude can fall because cellular amplitude falls or because cells drift out of phase; ORPHEUS separates the two. Data: mouse SCN scRNA-seq, mouse liver scRNA-seq (public time courses), human DLPFC snRNA-seq from ROSMAP, 409 subjects (250 aged controls, 159 AD dementia).
**Sampling:** mouse sets are timed, uniform full-day time courses (needed for the method); ROSMAP has no recorded time of death, so subject phases come from CYCLOPS 2.0 ordering (earlier work), not from clock time.

## Method (if any)
- Unit of inference: the single cell within a cell type and phase bin, pooled per gene; the output is a gene-wise phase dispersion (SD of cellular phase offsets, h), summarised by the median over rhythmic genes per cell type, condition or pathway bin.
- Principle: desynchrony imprints a 12 h rhythm on intercellular expression variance; a negative-binomial model solves for cellular amplitude, phase dispersion and overdispersion from pseudobulk and variance rhythms.
- Needs: time-course samples covering 24 h with enough points per cycle (bulk phases known or CYCLOPS-ordered). It does not need clock-gene coupling between cell types. "The algorithm currently relies on single-cell time-course data with relatively fine sampling to detect and accurately describe 12hr rhythms in intercellular variance."
- Accurate only for phase SD of 1-5 h (in silico).
- Cell types of the same donor are not compared: each cell type is analysed separately; within-donor cell-type offsets are not estimated.

## Key results
- Validation in mouse SCN: de novo dispersion 2.82 h (VIP+/Nms+, 10 genes) and 2.87 h (Cck+/C1ql3+, 7 genes), "in good accord with the past experimental observations" (PER2::LUC imaging).
- Mouse liver hepatocytes: median dispersion 3.09 h (67 genes); no difference between four hepatocyte subtypes; higher MTORC1 activity gives tighter synchrony (p=0.002, BH.q=0.088).
- ROSMAP excitatory neurons (subtypes 3 and 5), AD vs CTL: "We observed a significant reduction (p<2.2x10-16, Wilcoxen rank-sum test) in circadian synchrony within the AD excitatory neuron population" (3.76 h, 526 genes, vs 3.23 h, 464 genes); paired genes p=5.52e-14 (278 genes). The p-values are across genes, not across donors.
- Dose-response: dispersion rises with CERAD and with Braak stage (no numbers in text); within AD, MTORC1 (p=0.032, BH.q=0.26) and protein secretion (p=0.025, BH.q=0.26) activity bins differ; in controls no pathway does.
- OXPHOS genes: higher dispersion in AD (3.17 h, 10 genes) vs CTL (2.76 h, 5 genes), p=0.005.

## Not addressed / limits
- Disease contrast has no specificity control: one cell type (excitatory 3+5, chosen for strong clock signal), no other cell type or non-clock gene set as companion, no other disease. Gene-level tests treat genes as independent units.
- Dispersion is measured among cells given the subject phase; errors in the CYCLOPS subject phases would inflate dispersion (our V14 simulation: bias only at low dispersion, see Ideas). Authors do not model it in the human data.
- Dispersion has to be read together with coupling and gene number (526 vs 464 genes, different filtering), so the AD-vs-CTL gap of 0.53 h is small against the 1-5 h window.
- No donor-level inference, no between-compartment (neuron vs glia) comparison; liver data lack biological replicates ("pseudo-replication"), authors call it hypothesis-generating.
- Needs fine time sampling or a prior ordering; hospital-hours samples would not give the 12 h variance rhythm.

## Link to us
- **Intro claim(s):** refs_circadian_desynchrony.md lists it as within-tissue single-cell desynchrony with loss in AD excitatory neurons: supported. The cell-to-cell variance readout is accurate but it needs time-resolved data, which the one-line summary hides.
- **Relation to clocksync:** complementary and partly competitor. ORPHEUS: cells within one cell type, needs 24 h coverage, no donor reference. clocksync: compartments within a donor, untimed. Both report AD neuron loss of coordination; ours (A141, results_preliminary.md) finds the AD contrasts in MSSM, RADC and SEA-AD noise-driven, not estimable at current design, so no agreement claim is possible.

## Ideas
- Already done: ORPHEUS bias under subject-phase error (V14, METHOD.md); `compute_orpheus_cells.py` in extended_data calibration reads ROSMAP-type sets.
- Cite as the within-compartment axis next to our between-compartment axis; state that the two measure different things (cell phase spread vs compartment coupling).
- Compare: ORPHEUS dispersion in RADC excitatory neurons with clocksync coupling c of the same compartment in the same donors (A141 sets); a convergent or divergent result would be informative for the Discussion.
