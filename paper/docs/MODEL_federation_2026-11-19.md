# The body's clock as a federation: model statement (2026-11-19)

Take-home model of the paper. Supersedes `MODEL_hopf_2026-11-18.md`; the Hopf normal form survives as the cell-level
oscillator (section 2), not as the explanation of the organ-level data. The model is structural: each part is a claim
with its own registered test. A quantitative population fit of all parts at once was tried and is inadequate (A310).
Numbers are from the results entries cited (`results_preliminary.md`), mouse numbers from the literature cards
(`manuscript/literature/notes/`).

## 1. Statement

1. **Cells keep time autonomously.** The clock runs in each cell of the clock-carrying classes (stromal-vascular,
   epithelial, parenchymal, macrophages). Coupling between neighbouring cells is weak and permissive: it helps
   cells stay together but does not make the rhythm.
2. **Cell classes keep their relative loudness; composition moves gain, not phase.** Each cell class keeps the same
   gain relative to the organ's endothelium in every organ, and composition moves a sample's gain within a tissue, not
   its phase. Composition and conserved class gains do not explain the organ gain range: they spread organ gains
   1.3-fold where the organs spread 2.9-fold (A323). What sets the rest is open (biological or technical; check G1).
3. **An organ's timing is set by its drives.** Each organ follows a weighted mix of shared drives: a central drive
   while the brain lives, and behavioural drives (feeding, activity, light reaching the body by other routes). The
   weights are organ-specific. The organs' phase order is the order of these weighted mixes.
4. **The central drive sets the body's organ order.** Without it, organs follow the remaining shared drives in a
   different order (rearranged, not dissolved). With no shared drive at all, cells dephase and organs drift apart
   over days to weeks.
5. **Conflicting drives cost amplitude.** When the drives disagree, an organ moves by its weight on each drive and
   loses amplitude. Organs that follow the displaced drive least (mouse lung, heart, kidney) move 1.6-2.5 h and keep
   half to two thirds of their amplitude, more than a pure mix of drives predicts; the liver follows almost fully;
   adipose, with balanced weights, loses most (to a third or a tenth), at or beyond the mix (A322, A324).
6. **Each person adds a departure.** Around the shared state, each organ of each person carries a personal departure
   of nearly constant size, random in direction and not shared between organs. It is added variation, not the phase
   noise of an oscillator.

## 2. Formal version

For person p, organ o, cell class k:

    cell (class k):      dz/dt = (mu_k + i nu_k) z - |z|^2 z + sum_d w_{o,d} F_d(t) + eps * (local mean - z) + noise
                         mu_k > 0 (self-sustained), eps small (weak local coupling)
    organ state:         x_{p,o} = G_o exp(i phi_o) theta_p + delta_{p,o}
    gain:                G_o = D_o sum_k s_{o,k} a_k   (a_k conserved relative class gain, s_{o,k} composition,
                         D_o the organ gain not explained by them (A323); its nature is untested)
    lag:                 phi_o = arg( sum_d w_{o,d} exp(i Phi_d) )   (Phi_d: phase of drive d in internal time)
    amplitude under drive conflict:  R_o ~ | sum_d w_{o,d} exp(i Phi_d) | / sum_d w_{o,d}   (the mix of drives; lung, heart,
                         kidney stay above it, adipose at or below it, A322 / A324: not a bound)
    departure:           delta_{p,o} ~ CN(0, 2 psi_o), psi_o nearly constant across organs, independent between organs

- theta_p, the person's internal time, is what clocksync reads as the host. Coupling c_o = S_o / (S_o + 2 psi_o) with
  S_o = G_o^2 V follows from gain and departure: strong clocks look more coupled because the departure has a fixed
  size, not because they are stiffer oscillators.
- Removing a drive d sets w_{o,d} = 0. Because mu_k > 0, amplitude is held by the cells (saturation) and only phase
  moves: drives are not additive (A297).
- Brain death removes the central drive: phi_o is then set by the remaining drives with their own weights (new order);
  R_o stays near 1 while some shared drive remains; with none, cell phases spread at rate set by period dispersion
  and noise, and G_o falls with time (mouse: days to weeks).

## 3. Evidence for each part

| part | human (entries) | mouse (cards / entries) | status |
|---|---|---|---|
| 1. cells autonomous, weak local coupling | stromal-vascular, epithelial, parenchymal cells and macrophages carry the organ's time; tissue lymphocytes no readable clock (A272-A288); adjacency does not add clock-specific sharing (A313 FAIL) | hepatocyte clocks coherent without SCN (Sinturel 2021); hepatocytes shape other cells' outputs, not their clocks (Guan 2020); weak paracrine coupling (Finger 2021, TGF-beta; Noguchi 2013, density); dissociated peripheral cells not rescued by coupling (Liu 2007) | supported; local coupling not measurable in human (A313) |
| 2. gain = class x composition (x unexplained organ term) | class gain relative to endothelium conserved across organs, within a factor of about 1.4 (A293 H293.4); composition moves gain, not phase (within 0.13 h per SD, A295); gains 9.2-fold apart across tissues (shared variance S 85-fold; A294b); composition x class gains predict a 1.3-fold organ range against 2.9-fold observed over 8 tissues (A323 FAIL) | none at this resolution | class part supported; composition and class gains do not explain the organ gain range |
| 3. timing set by weighted drives | cell classes share the organ phase, about 0.6 h (A273, A279); night shift: muscle moves 2.7 h after a 12 h behavioural shift (A257) | feeding sets peripheral phases, weights differ by organ (Manella 2021; Izumo 2014: liver and kidney follow feeding, heart, lung, spleen less) | supported; one weight per organ not yet fitted blind (gap 5) |
| 4. central drive sets the order | order kept without alignment to the external day (slow deaths) and without meals (ventilated living-brain donors) (A293 Part 1); different in brain-dead organ donors, clock-specific (A317 H317.1, D_rho 0.621), after ischaemic time and RIN adjustment (A318, 0.584; registered re-tests, not blind), aliased with collection at procurement; shared clock between organs kept (A317 H317.2 FAIL); amplitude equal to equally ventilated living-brain donors (A315, clock-specific ratio 1.05) | SCN-clock knockout: order kept under light-dark (Husse 2014; Izumo 2014); under feeding, a reproducible new order (A319 H319.2, order share 1.09, 6/6 organs cluster); in constant dark, organs drift and dampen (Husse: by day 7; Izumo: > 30 days; A319 H319.1, no significant shared phase, p 0.070); SCN lesion: organs oscillate at lower amplitude (Tahara 2012) | supported in GTEx (single cohort) and mouse (explants, not blind); replication open (gap 1) |
| 5. conflict costs amplitude | night shift: muscle amplitude to about a quarter (A257 H257.3, R 0.233); forced desynchrony in blood (A262); mouse feeding inversion (two datasets): liver 0.44-0.71, lung / heart / kidney 0.50-0.65, adipose 0.11-0.32 of amplitude; a pure mix of drives is rejected in both (A322, A324) | the central drive slows feeding-driven shifts (Saini 2013: liver re-phases in 1-2 d without SCN, 6-7 h after 7 d with it); removing feeding keeps 0.72 of liver amplitude, inverting it moves the liver 11.7 h (A297) | supported |
| 6. personal departure | clock-specific, 2.6 times expression-matched genes, nearly constant size (A294b); not shared by meal-driven tissues (A302); not phase-dominated (A311 FAIL, median orientation share -0.028); persons differ in shared amplitude, every tissue follows it nearly in proportion (A312: rho -0.341, elasticities 0.83-1.77) | none | supported as added variation; its cause is open |

## 4. The paper's results in the model

| result | model reading |
|---|---|
| Lymphocytes carry no readable clock | not a clock-carrying class (a_k near 0); they share hormone and immediate-early states (A275), i.e. they feel the drives without a running loop (core loop open, gap 7) |
| Fetal organs: a driven state before an adult-size loop (A242-A269, A307, A309) | drives reach the organ before its cells run a loop (mu_k < 0 during development) |
| Chronic disease changes nothing (A164, A233, A235, A306) | no parameter moves |
| Psoriasis delays the lesion clock (A183, A308) | the lesion's drive mix or cell periods change (phi_o moves locally) |
| Tumours keep their patient's time and lose gain as they proliferate (A138, A304) | the tumour still follows the patient's drives; its gain falls (fewer clock-carrying cells or weaker cell clocks; not separable in bulk, A212, A305) |
| Age and sex leave couplings and offsets equal (A208, A231) | drives and weights conserved |

## 5. Tested limits (what the model is not)

- **Not a fitted population model.** The Ott-Antonsen population fit (A310) was inadequate: about a third of tissues
  inside their intervals, one noise intensity does not fit all tissues. The arithmetic bound leaves composition able
  to explain organ gain except in blood and testis.
- **Not oscillator stiffness.** A Hopf cycle above the bifurcation predicts tangential (phase) departures; observed
  departures are isotropic or slightly radial (A311 FAIL). Gains follow the shared amplitude nearly in proportion
  (A312): adult organ gains are at most weakly a matter of distance to the bifurcation.
- **No measurable local coupling in human.** Adjacent tissues share more state, but the gradient is the same in
  expression-matched genes (A313 FAIL).
- **Amplitude kept only in the short term.** Brain-dead donors are measured within hours to days, in ICU conditions;
  mouse tissues lose amplitude after days to weeks with no drive (Husse 2014, Izumo 2014, Tahara 2012).
- **Order kept between persons, not within one.** All statements are population-level estimates of within-person
  relations.

## 6. Predictions not yet tested

1. **Drift grows with time since loss of the central drive** (MODEL prediction 1; formerly synthesis gap 2). Organ disorder rises with time since brain death
   (GTEx has ventilation duration only; A132 tertiles in the drift direction, not significant).
2. **The unexplained organ gain.** Whether D_o is biological (drive reaching the organ, or the organ's cell clock
   strength) or technical (tissue-specific post-mortem attenuation in GTEx): check G1 (REVIEW_2026-11-19 section 5).
3. **One weight per organ.** A single b per organ predicts both the shift and the amplitude loss under any drive
   conflict (gap 5; Manella liver, lung and fat as the not-blind check; Saini 2013 as the mouse liver calibration).
4. **Within a brain-dead donor, cell types of an organ still share its phase** while the organ order differs from
   living-brain donors (Tabula Sapiens, gap 1).
5. **Tumours in single cells.** Malignant cells carry a weaker cell clock or fewer clock-carrying cells; not cells
   drifting out of phase.

## 7. Outside the model

The brain-body split (A263, A289: two domains, not clock-specific); the GTEx adrenal offset (read as an effect of
death, A296-A303); waveform and harmonics; the cause of the personal departure.
