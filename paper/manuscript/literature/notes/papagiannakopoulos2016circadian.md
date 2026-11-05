# papagiannakopoulos2016circadian
**Ref:** Papagiannakopoulos T, Bauer MR, Davidson SM, Heimann M, Subbaraj L, Bhutkar A, Bartlebaugh J, Vander Heiden MG, Jacks T, Cell Metab 2016;24:324-331. DOI 10.1016/j.cmet.2016.07.001. Type: biology (mouse GEMM). Lab: Jacks.
**Source:** full text

## Question and design
Autochthonous lung adenocarcinoma, KrasLSL-G12D (K) and KrasLSL-G12D;p53flox/flox (KP) mice, plus KrasLA2/+; tumours induced by intratracheal Cre. Host clock perturbed in three ways: jet lag (8 h light advance every 2-3 days), germline Per2m/m or Bmal1-/-, and tumour-specific Bmal1 deletion (Bmal1flox/flox). Tumour clock not measured (no timed tumour expression); tumour and host not compared. **Sampling:** endpoint histology about 13 weeks after tumour initiation; no time-of-day series; TCGA lung adenocarcinoma expression untimed. n: KP LD 6 / jet-lag initiation 7 / progression 7 (histology); survival 15 / 13 / 13; KP Per2m/m 12 vs 12 (histology).

## Method (if any)
Histological tumour burden and grade, Kaplan-Meier, c-Myc and BrdU IHC, metabolite flux with [U-13C]glucose, CRISPR Per2 knockout lines, GSEA on tumours ranked by PER2 (top and bottom 25 %, n = 12).

## Key results
- Jet lag: "Animals placed in jet lag during tumor progression for 8 weeks post-infection had a significant increase in the lung tumor burden" (p < 0.05); survival decreased, median 9 days (p < 0.01). Starting a week before initiation gave "a modest but non-significant increase in tumor burden" but more grade 3-4 tumours (p < 0.05).
- Activity was disturbed ("we confirmed that jet lag disrupted the circadian behavior of animals", actograms), but no peripheral or tumour clock-gene rhythm was measured here; disruption of core genes in tissues is taken from refs.
- Germline Per2m/m: higher burden and grade in K and KP, survival decreased (KP p < 0.001); KrasLA2/+ with Per2m/m or Bmal1-/- more lung tumours (p < 0.05, p < 0.01) and worse survival; no effect on other sporadic Kras tumours.
- Cell-autonomous: "Loss of Bmal1 (Bmal1D/D) in tumors accelerated lung tumor progression in the K model" (n = 17 vs 14) but "there was no increase in tumor burden in KP tumors with conditional Bmal1 (Bmal1D/D) loss", p53-dependent.
- Per2-null KP cells grow faster in vitro; CRISPR Per2 knockout in isogenic C57BL/6 KP cells increases subcutaneous tumour formation in immunocompetent hosts.
- Mechanism: higher c-Myc by IHC in K;Per2m/m (p < 0.05), K;Bmal1D/D (p < 0.01), KP;Per2m/m (p < 0.0001); increased glucose and glutamine use in Per2m/m KP cells.
- Human: "We present evidence for decreased expression of most circadian genes in human lung adenocarcinoma patient tumors as compared to matched lung tissue" (TCGA, 57 pairs); CLOCK up.

## Not addressed / limits
- Does not separate host-systemic from tumour-intrinsic effects for jet lag and germline mutants; the authors note jet lag acts via "desynchronization of the circadian machinery" among many other routes (hormones, feeding, inflammation), none isolated.
- Tumour clock oscillation is never measured; low expression of clock genes in human tumours is an average, not a rhythm.
- Context-dependent: Bmal1 loss helps K but not KP tumours; Puram 2016 shows loss suppresses leukaemia.
- Small n in several arms (KP Per2m/m survival n = 5; Bmal1-/- n = 7).

## Link to us
- **Intro claim(s):** the refs_*.md row is supported (jet lag and Per2/Bmal1 loss shorten survival and promote growth) but the "tumour-intrinsic" reading is stronger than the data for jet lag; tumour-autonomous effect is shown only for Bmal1 deletion (K model) and Per2 loss in cells.
- **Relation to clocksync:** background / data we could not use (no expression time series); justifies measuring tumour and host clocks separately.

## Ideas
- Direction for human data: tumour clock-gene attenuation versus matched normal is already tested (TCGA pairs, A138/A85-A87 in paper/docs/results_preliminary.md); cite them with this paper.
- Missing reference: Puram 2016 (Cell, leukaemia stem cell clock), Altman 2015 (Mol Cell, MYC disrupts clock), Hadadi 2020 (chronic jet lag in a spontaneous mammary tumour model; scouted for A125 in paper/docs/provisory_plan_paper.md, not yet registered).
