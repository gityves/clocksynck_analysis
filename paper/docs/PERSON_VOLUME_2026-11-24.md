# The person's clock volume: what it is, what it is not (2026-11-24)

Planning note for `prompts/person_volume_desynchrony_2026-11-24.md` (task T3), after the assessment of that prompt and
the amendments A361-A363 and their gene-drop addendum. Numbers from `results_preliminary.md` entries cited. Proposed
wording only: the paper plan and the slides are not changed.

## 1. Question

Persons differ in how loudly all their organs express the shared clock time (the person volume). The prompt read a
quiet person as a candidate for body-wide desynchrony below the resolution (cells or conflicting cues), and proposed
heritability (MuTHER twins, T1) and exposures (GTEx, T2). The assessment found that T1 and T2 cannot separate that
reading from systemic damping of cell clocks, composition or the death process, that T1 is not decidable at MuTHER's
design (101 MZ and 154 DZ pairs; h^2 90 % range -1.3 to 2.3 at a true 0.5) and that MuTHER's centred radius can read
the shared biopsy time rather than loudness (synthetic leverage 0.18-0.29 with no amplitude variation at all). T1 and T2
were not registered; three discriminating tests were run instead.

## 2. Chain, with outcomes

| link | evidence | status |
|---|---|---|
| persons differ in a volume shared by their organs (deaths) | A312 leverage 0.75-0.86; A336 part 2 core 0.714 above matched q95 0.621 and responsive q95 0.692 | have (deaths) |
| every organ follows it nearly in proportion | A312, level calibrated post hoc (estimator reads about 1.1 under proportional scaling) | have |
| it scales the core clock genes alike | A363 part 2: E-box / PAR-bZip targets vs positive limb T -0.006 (recalibrated null -0.087 to 0.092) | have; a 1.5-fold cytokine-like distortion excluded, a 1.25-fold one not (power 0.65) |
| it falls with the glucocorticoid-target state | A360 (14 genes -0.235); A361 core -0.184 (90 % -0.297 to -0.063), beyond all 400 companions; gene-drop: the trio drop is ordinary (85th percentile) | have |
| inflammatory state | A363 part 1: core FAIL (p 0.30, power 0.76); 14 genes p 0.034; core value a typical 11-gene value (82nd percentile) | not decided |
| organ pattern of cue conflict | A362 FAIL (power 1.00): slope -0.09 against null q95 0.27; kidney follows most, lung least | not seen |
| living persons | A331 0.177 at the line on the core loop (0.107 < matched q95 0.130) and possibly biopsy timing; del Olmo 2022 and Saint-Antoine & Anafi 2025 (preprint) report person-level amplitude in timed skin | open: A364 registered |
| acute vs chronic glucocorticoid state | deaths only | open |
| genetic basis | not testable with public data at useful power | out of reach |
| cells desynchronised rather than damped | A334: public 10x cannot decide | out of reach |

## 3. Reading

The volume is a systemic, body-wide amplitude factor that every organ and every core clock gene follows, tied to the
glucocorticoid-target state. Its organ pattern is not the one cue conflict produces in mouse, and its gene pattern is
not a strong cytokine damping. What it physically is (uniform damping of cell clocks, cells dephasing, body-wide
composition) is not separable in bulk. In deaths, an agonal glucocorticoid surge and a chronic state are confounded.

## 4. Proposed wording

**Paper (results, figure 4 or 5 panel):**
> Each person has a clock volume that every organ follows and that scales the core clock genes alike. It is lower when
> the body's glucocorticoid-target state is higher, beyond any generic donor scale, and its organ pattern is not the one
> conflicting cues produce in mouse.

**Discussion:** one sentence that a lower volume could reflect damping or dephasing of cell clocks body-wide, which
bulk data cannot separate; "desynchrony" only as one of these possibilities. Remove "as conflicting cues produce".

**Slide 3.5 (proposed edits, not applied):**
- Bullet 2 (living women): move to the caveats until A364 is scored; MuTHER is at the line on the core loop and its
  radius may read the shared biopsy time. If A364 passes, replace it with the timed living result.
- Bullet 4: "The volume is lower in persons with a higher glucocorticoid-target state, in the core clock as well
  (standardised effect -0.18, 90 % -0.30 to -0.06; beyond all 400 companion gene sets), not explained by RIN, ischaemic
  time or time of death (A360, A361, registered, blind)."
- New bullet 5: "It is a body-wide scale, not a pattern: all core clock genes follow it alike (A363), and organs do
  not follow it by how much conflicting cues cost them in mouse (A362)."
- Panel C: add the core-loop GC coefficient beside the 14-gene one.
- Caveats: the gene-drop result (the smaller core coefficient is an ordinary three-gene drop); A363 part 1 at power 0.76.

## 5. Out of reach

Cells (single-cell clocks in deep data); a within-person repeat design (the same person sampled when quiet and loud);
heritability (MuTHER underpowered and timing-confounded; GTEx genotypes would need an access request and are far too
few); exposures with a useful count in GTEx (long-term steroids 5 of 244 fast deaths, depression 14).

## 6. Next

A364 (registered 2026-11-24): living timed skin, GSE326530: is the person's clock amplitude shared by two skin sites,
measured per volunteer from a full timed day (no biopsy-timing confound)?

## 7. Living link and data (added 2026-11-24, after A364)

- A364 (GSE326530) is not decidable: power 0.54 at the GTEx spread (20 persons, 4 noisy samples per site).
- Literature, carded (`manuscript/literature/notes/`): `delolmo2022interlayer` (11 adults, epidermis and dermis every 4 h:
  amplitude varies more across subjects than layers) and `saintantoine2025defining` (preprint: a subject's global
  amplitude is shared by epidermis and dermis, R^2 0.649 in del Olmo's data and 0.345 in Wu 2020; not strongly tied to
  core-clock amplitude). Wording: "a person-level rhythm amplitude shared by skin layers has been reported in timed
  living skin (del Olmo 2022; Saint-Antoine and Anafi 2025, preprint)".
- Data: del Olmo GSE205155 is public (raw Agilent TXT, 154 samples, 7 times per person and layer; not downloaded). A
  registered clock-plane test with companions on it would not be blind for the global amplitude (the preprint's R^2) but
  would be for the core-clock, companion-controlled version.
- Mouse GSE201207 (Wolff et al., Cell Rep 2023, DOI 10.1016/j.celrep.2022.111982): six tissues (hypothalamus, lung,
  heart, kidney, muscle, adrenal) x three ages (6, 18, 27 months) x 12 circadian times (every 4 h for 48 h), male, counts
  public. GEO design: "Two samples per age were collected and pooled for RNAseq analysis": each sample pools two mice,
  one pool per age and time; per-pool tissue identity is not stated in GEO (paper methods not read).

## 8. Living tests and the mouse ageing atlas (added 2026-11-24, after A365)

- A365 (del Olmo GSE205155, already on disk from archive A37) is not decidable either: power 0.51. The clock signal is
  above the noise, but 11 subjects and two layers of one punch (matched genes reach r 0.67 at q95) leave a GTEx-sized
  volume indistinguishable from the shared biopsy scale. Both public timed living designs are too small: the living link
  rests on the literature (section 7), with the caveat that it has no matched-gene control.
- Mouse GSE201207 (Wolff et al., Cell Rep 2023; methods read, PMC9929559): C57BL/6 Nia males of 6, 18 and 27 months,
  n = 24 per age, "2 animals for each age were collected at each time point", every 4 h for 48 h from CT18 in constant
  darkness (time in darkness before harvest not stated); poly-A RNA-seq, 2 x 100 bp, >= 40 M reads. GEO lists one
  sample per tissue x age x time with "Two samples per age were collected and pooled": the two animals of a time point
  are pooled, and whether one pool supplies all six tissues is not stated. Aged adrenal lacks the first 24 h. No
  amplitude results by age in the paper (rhythmicity counts only; Per1 / Per2 "dampened from Aged to Old across all
  peripheral tissues").
- What GSE201207 can test: not a person (animal) volume shared by organs (pools, identity across tissues unknown), but
  whether the clock amplitude of all peripheral organs falls together with age, in timed data (organ-shared amplitude
  as an age state). Group-level, 5 peripheral tissues x 3 ages x 12 times.

## 9. MuTHER resolved (added 2026-11-24, A366)

A366 FAIL (powered): MuTHER's states lie on a daytime arc (elongation 3.4, phase SD about 1.9 h); measured from the
ring's centre, the fat-skin leverage is 0.165, inside what shared biopsy time alone gives (timing-only q95 0.209, median
0.151), far below a GTEx-size (0.314) or half-size (0.262) volume; radii rise with the distance from the mean phase
(Spearman 0.21 / 0.19). A331 is withdrawn as living evidence; slide 3.5 bullet 2 and Figure 4f's living panel should go
(proposal, not applied). The person volume is shown in deaths; in living people it is reported only by others, without
a matched-gene control.

## 10. Matched-phase test (added 2026-11-24, A368)

A368 FAIL (powered 1.00 at GTEx size, 0.87 at half): at matched arc position, fat-skin amplitude correlation Q 0.228 is
above everything biopsy timing produces (q95 0.163) but below expression-matched random genes (q95 0.253), and at the
4th percentile of a half-GTEx-size volume. A GTEx-size person volume is not present in living MuTHER women. Either the
volume in deaths is larger than in life (agonal or post-mortem amplification: A367), or living and dead cohorts differ
in what they measure (one biopsy of two tissues vs 14 tissues per donor). The claim stays "in deaths".

## 11. STAGE (added 2026-11-24, A369)

A369 not decidable (power 0.04): 66 bypass patients, liver / muscle / fat from one operation, phase SD about 4.4 h,
noisy muscle and fat clocks. Synthetic check: no statistic decides at 66; about 600 persons would (STARNET). The public
living designs are exhausted for this question; STARNET (dbGaP, user application) is the remaining living test.
