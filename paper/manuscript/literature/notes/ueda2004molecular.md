# ueda2004molecular
**Ref:** Ueda et al., PNAS 101:11227-11232 (2004). DOI 10.1073/pnas.0401882101. Type: method. Lab: Hiroki R. Ueda (Yamanouchi / RIKEN).
**Source:** full text (main; supporting figures and tables not in the export)

## Question and design
Detect the body time (BT) of an individual from one genome-wide expression profile. Mouse liver (Balb/c, male), Affymetrix
U74A microarrays; also SCN and Drosophila head re-analysis of published series. **Sampling:** the timetable is built from a
timed full-cycle series (pooled livers of four mice, 12 time points every 4 h over 2 d, under LD and under DD). Test livers
are single time points from individual mice: ZT12 (n = 4), ZT6, ZT18, CT6, CT18 (n = 1 each), plus 3 C3H mice at ZT12,
3 wild type and 3 Clock/Clock at ZT8, 4 Clock/Clock at ZT12. All animals entrained to LD for 2 weeks before sampling.

## Method
Supervised, cohort-trained, single-tissue. (1) Select time-indicating genes from the timed series: best-cosine correlation > 0.8
in both LD and DD and CV > 0.15 (168 genes, 182 probes); each gets a molecular peak time. (2) For a new profile, z-score each
gene with the timetable mean and SD, correlate over genes the profile against 24 h cosines indexed by the gene peak times
(144 phases, 10 min steps); the best phase is the BT, the best correlation c a rhythm score. (3) Significance from 10 000
random profiles (P_r); rhythm disorder from 10 000 control profiles with 100 % noise (P_c). Needs timed training data of the same
tissue and species. A test sample needs no time coverage: one profile gives one BT.

## Key results
- Peak times span the day. "Importantly, the molecular peak times of time-indicating genes were distributed over 24 h (Fig. 1B)."
- Training-series self-test: BT error 0.0-1.3 h (mean 0.4 +- 0.3 h) (Fig. 1D).
- Independent mice, ZT12 (n = 4): BT 10.2-12.8; "These results suggest that BT can be accurately inferred from individual expression profiles with estimation errors from 0.3 to 1.8 h" (mean 1.0 +- 0.6 h).
- Other times and DD, 4 mice: errors 0.3-2.0 h (mean 1.1 +- 0.7 h). "These results suggest that BT can be accurately inferred from expression profiles at time points that were not used to construct the molecular timetable".
- Other strain (C3H, n = 3, ZT12): errors 0.3-0.7 h (0.6 +- 0.2 h).
- Rhythm disorder: Clock/Clock mice flagged at ZT12 and ZT8 (P < 0.0001 in all); wild type at ZT8 estimated BT 7.3, 6.7, 9.5. "To demonstrate the feasibility of expression-based diagnosis of circadian rhythm disorders, we sampled livers at" ZT12 from four mutants.
- Simulation: with 100 genes and 200 % noise, sensitivity and specificity 98.58 %, error 0.62 +- 0.48 h. "Even in the presence of 200% measurement noise, use of the molecular-timetable methods with 100 (or 150) time-indicating genes still had the capability of detecting BT with" high accuracy; 3 or 5 genes fail.
- Re-analysis: SCN timetable errors 0.3-1.8 h; Drosophila head (113 genes) 0-1.5 h, dClock mutant flagged.

## Not addressed / limits
- Mouse only, one tissue per timetable, tiny test set (8 wild-type livers for BT accuracy). Human not tested; "can be also applied to organisms other than the mouse" is argued from other species' series.
- Test mice are inbred, same age, LD-entrained, sampled at a few ZT/CT. No shift, disease or tumour.
- The rhythm-disorder call is a single threshold on c, without a coupling between tissues. No second compartment of the same animal is compared; the SCN timetable is built separately.
- Phase accuracy depends on timetable gene peak times being fixed across individuals.

## Link to us
- **Intro claim(s):** early supervised phase inference from a single sample (technical challenges): our refs entry (supervised, timed full-cycle training, single test sample needs no coverage) is supported.
- **Relation to clocksync:** background. The cosine-over-gene-peak-times score is an ancestor of CYCLOPS-type and gate scores; it assigns clock time and needs a timed reference for each tissue, whereas clocksync never assigns time and reads compartments against the donor's shared state.

## Ideas
- The P_r / P_c random-profile null is a model-free "clock-shaped" test on a single sample; compare with our gate (METHOD.md gate rows) as a simple baseline. Not done in registered_experiments.md (no hit for Ueda or molecular timetable).
- Missing reference check: Hughey 2016 ZeitZeiger (in bib), Braun 2018 TimeSignatR (in bib).
