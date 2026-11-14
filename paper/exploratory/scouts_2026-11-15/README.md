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
