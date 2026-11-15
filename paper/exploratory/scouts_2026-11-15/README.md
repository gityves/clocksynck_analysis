# Scouts for phase T (T0b, T0c), 2026-11-15

Roadmap items T0b (power of T1 to T3 at the program level) and T0c (calibration set). A scout reads designs, reference
data, noise and simulated data; it does not compute any statistic a later amendment registers on the real labels.
Gate G1: an item is registered only if its scout shows it can separate the times at the program level.

| item | scout | G1 | next |
|---|---|---|---|
| T0c calibration (mouse liver with and without a local clock) | `T0c_calibration/` | n/a (calibration, needed for G2 (ii)) | register (A260) |
| T1 night-shift muscle, programs | `T1_nightshift/` | **pass**: central / tissue clock / behaviour separable from 20 genes per program | register (A261) |
| T2 GTEx, body vs tissue vs day per gene | A259 (registered; `A259_tissue_deviation/`) | **pending A259**: the tissue term needs a person-specific deviation above noise | decide after A259 is scored |
| T3 repeat biopsies, rotation | `T3_repeat_biopsy/` | **fail as a separator of times** | not registered |
| T1b clock-autonomous core in night-shift muscle (follow-up of A261, 2026-11-15) | `T1b_clock_core/` | **fail: not estimable** (gate needs about 600 genes; a damped clock-following program is never told from central time) | not registered |
| T5 second misalignment design (G2 (iii)): blood, GSE48113 and GSE107537 / GSE122725 | `T5_blood_misalignment/` | **pass** for the genome-wide program (>= 0.965 correct at every attenuation); muscle or adipose: no misaligned design public | register (A262) |

## T0c: calibration set (`T0c_calibration/scout_t0c.py`, `output/`)

Manella 2021 (GSE159135): Alb-Cre control and hepatocyte Bmal1 knockout, ad libitum (AL) and daytime feeding (DF),
liver, WAT and lung, 12 ZT x 2 mice per arm (24 mice; the A111 assumption that a label is one mouse in the three tissue
files holds here too). Timed cosinor per arm and tissue (`output/design.csv`): rhythmic genes at q < 0.05 in the liver
1,595 (control AL), 99 (knockout AL), 1,518 (control DF), 588 (knockout DF); WAT 1,018 / 209 / 681 / 724; lung
1,178 / 1,035 / 55 / 93. Liver classes from the recorded times, non-clock genes (`output/classes.csv`): **local**
(rhythmic in control, lost in the knockout) 584, **systemic** (kept, amplitude ratio >= 0.5) 40, knockout-only 26.
Daytime feeding moves the clock genes of the control liver by -11.56 h, of WAT by -9.24 h and of lung by -2.31 h
(median per-gene shift; knockout WAT -9.57 h, lung -2.34 h; `output/shifts.csv`). Consequences for the calibration:
- the knockout gives a truth for "body vs local" (584 vs 40 genes): the body reader (untimed host loading, host from
  WAT and lung) can be scored on it;
- feeding moves WAT with the liver, so the host for a feeding contrast is the lung alone;
- nothing in the control arms dissociates the liver clock from the body time person by person (the liver is coupled at
  c 0.967, A129): the tissue-deviation reader cannot be calibrated on a natural deviation here, only described (A259
  estimator on the control liver).

## T1: night-shift muscle (`T1_nightshift/scout_t1.py`, `output/`)

Design GSE106800: 11 subjects with 07:00 and 19:00 biopsies in both conditions. Reference: timed muscle GSE108539 (10
subjects x 6 biopsies), per-gene cosinor with subject intercepts: 1,896 genes at p < 0.01, 1,597 of them on the array,
median amplitude 0.274 log2 (`output/noise_summary.csv`). Noise: residual vectors of the subject differences around the
condition means (median SD 0.317 aligned, 0.331 misaligned, against a median signal 2A of 0.548), resampled jointly over
genes. Hypotheses: central (s 0 h, r 1), tissue clock (2.68 h, 0.233: A257's muscle clock), behaviour (12 h, 1).
Estimator: per subject, regression of the program's differences on 2A (cos, sin)(w (19 - phi)); a hypothesis is
excluded when the Hotelling T2 of z_mis - r e^{iws} z_al over subjects rejects. 200 programs per cell, random or
phase-clustered genes, between-cohort phase error SD 1 h.

At level 0.01 (`output/separability_by_level.txt`; random and phase-clustered programs pooled, 5 to 100 genes): the truth is excluded in 0.028-0.045 (central),
0.050-0.092 (tissue clock), 0.018-0.045 (behaviour); the alternatives are excluded with power 0.90-1.00 at 10 genes
except tissue clock vs central (0.805 at 10, 0.915 at 20, 0.958 at 50). At 0.10 the truth is excluded too often
(0.08-0.25): the registered level must be 0.01. Most of the separation between central time and the tissue clock comes
from the excursion ratio (0.233 vs 1), not the angle. 39 Hallmark sets hold >= 10 reference-rhythmic genes, 28 hold
>= 20 (counted before the array overlap). **G1 pass** for programs of >= 20 genes on the array, at level 0.01.
Caveats for the registration: not blind (A124, A257 scored these arrays: clock genes only); the reference cohort ate
hourly while awake, so meal-driven rhythms are under-represented in its rhythmic set; a program that mixes times
excludes all three or none.

## T2: GTEx body vs tissue vs day (A259)

The two identification problems of the roadmap are one quantity: how much of a tissue's departure from its person's
body time is real and person-specific rather than measurement noise. A259 measures it with disjoint gene halves and an
unshrunk host (`clocksync.deviation`), so the host's measurement error cannot leak into the tissue term. If a tissue's
deviation is not detectable above its detection limit, a gene's loading on it cannot be told from a loading on the body
time and the T2 tissue term is not registered there; the two-way A91 (body vs day, rerun unchanged on current code)
then stands as the GTEx part of the partition, with companions added.

## T3: repeat biopsies (`T3_repeat_biopsy/`)

Designs (`design_inventory.md`, compiled by a Sonnet agent from sample sheets and the A252 / A256 design tables, checked
against the A256 registration): every same-day repeat-biopsy study on disk was scored in A252 or A256; arms with nothing
done between the biopsies are GSE40589 saline (adipose, 7 persons, 4.0 h) and GSE59088 control (muscle, 6 persons,
2.5 h): 13 persons, the same units that A256 left unscored (H256.1). Every other arm contains a meal, exercise, an
insulin clamp or surgery. Power (`scout_t3.py`, optimistic: independent gene noise): a program rotating at the muscle
reference amplitude is told from a flat one with power 0.80 (muscle arm) to 0.99 (adipose arm) at 10 genes
(`output/power.csv`). **G1 fail as a separator of times:** in a rest arm the body time and the external day advance by
the same hours, so rotation shows that a program follows time, not which time; in the intervention arms rotation and the
intervention's response are confounded, and the rest arms are too few to subtract it. T3 is not registered. What it
could add is a living replication of "rhythmic vs flat", which the timed muscle reference already gives.

## T1b: the clock-autonomous core of night-shift muscle (`T1b_clock_core/`, 2026-11-15)

Question (follow-up review of A261): do direct clock targets stay with the muscle clock (2.68 h, damped) while the rest
of the rhythmic transcriptome follows behaviour? Sources of an external core set (`clock_target_sources.md`, Sonnet
scout, no human muscle data opened): the rule would pick the tissue-matched knockout time series first (Dyar 2014,
GSE43071, muscle-specific Bmal1 KO, TA and soleus, 6 ZT x 3), then muscle BMAL1 / CLOCK ChIP (GSE143334), then the A260
liver class (clock-dependent in the hepatocyte knockout, tissue mismatch); none ships a ready gene list, and no
mouse-human orthology table is on disk (the project convention is the upper-cased symbol).
Power first, on matched stand-ins only (`scout_core_power.py`; the core genes are never read on the arrays): programs of
reference-rhythmic genes outside the core, the A261 rhythm gate (matched non-rhythmic q95) and the A261 Hotelling
classes, with each truth imposed on the stand-in's own misaligned residuals (`output/power.csv`, 300 draws per size):

| genes | 20 | 50 | 100 | 200 | 400 | 600 | 800 | 1,000 |
|---|---|---|---|---|---|---|---|---|
| gate passed | 0.000 | 0.000 | 0.003 | 0.033 | 0.157 | 0.827 | 0.930 | 0.970 |
| gated and "tissue clock" called, when true | 0 | 0 | 0 | 0.003 | 0.010 | 0 | 0 | 0 |
| gated and "behaviour" called, when true | 0 | 0 | 0.003 | 0.033 | 0.157 | 0.827 | 0.930 | 0.970 |

The A260 liver class (587 mouse genes) holds 65 reference-rhythmic genes on the array: gate 0.000 at that size
(`output/power_liver_local_A260.csv`). Under a damped clock-following truth the call is "unresolved" (central not
excluded) in 0.90 to 1.00 of gated stand-ins at every size: an excursion ratio of 0.233 on a program read at about a
fifth of its reference amplitude (aligned loading 0.222, A261) cannot be told from central time with 11 subjects.
**G1 fail: not estimable.** A core of the size the literature describes (10 to 14 % of rhythmic genes in liver; at most
a few hundred in muscle) stays far below the gate, and even a gated core could only show "with behaviour" or not, never
"with the clock". Not registered, not scored; no core gene set was read on GSE106800.

## T5: a second misalignment design (`T5_blood_misalignment/`, 2026-11-15)

Muscle or adipose with shifted meals or sleep: nothing on disk (`repeat_biopsy/`: every same-day design has meals,
exercise, clamps or surgery at fixed clock times; GSE118280 swaps meal composition, not meal times). Public (literature
agent, `T1b_clock_core/literature_core_clock_2026-11-15.md`): Lundell 2020 (GSE129843, muscle every 4 h,
time-restricted 10:00-18:00 vs extended 07:00-22:00 feeding) changes the feeding window, not its phase, and does not
move sleep, so clock and behaviour are not misaligned; Wilms 2019 (adipose, sleep loss) has no accession found;
Christou 2019 (GSE87761, adipose, constant routine) has no misalignment; Wehrens 2017 has no genome-wide data. **No
second non-blood misalignment design with public matrices.**
Blood (on disk): GSE48113 (Archer 2014, forced desynchrony: 22 / 21 subjects, 7 samples on day 1 and day 4; melatonin
at the same clock time on both days, sleep and meals displaced by about 12 h) and the dim-light arm of GSE122725
(= GSE107537, Kervezee 2018: 8 subjects, every 4 h at baseline and after 3 nights with sleep delayed by 10 h;
melatonin delayed about 1 to 2 h; measurement under constant posture with hourly snacks, so no meal response inside the
sampling window). A156 read GSE48113 for the clock and cue genes only; the genome-wide program of either dataset has not
been read. Reference for central-time blood rhythms: GSE39445 (Moller-Levet 2013, sleep-extension arm, extended
wakefulness, melatonin-referenced phase, linear hours-awake term), 8,240 of 15,632 genes at p < 0.01 (blood rhythms
include the diurnal cell composition, an expected exception for the clock readout, not for this question).
Power (`scout_blood.py`, simulation on the real designs, noise = within subject x condition SD with the rhythm left in,
`output/power.csv`): for the genome-wide program the three outcomes central / behaviour / lost are called correctly in
0.965 to 1.000 of programs in both designs at reference amplitude x 1, 0.5 and 0.25; at 20 genes and x 0.25 only 0.01
to 0.37. **G1 pass for the genome-wide program**; not for Hallmark-size programs (the A261 lesson). Caveats: the
sampling window of GSE48113 contains the displaced sleep and meals (masking possible, but a full cycle is sampled);
GSE122725's central time is known only at group level (melatonin delay 1 to 2 h, Kervezee 2019); blood has no core clock
for us, so the question is only which time the rhythmic program follows.
