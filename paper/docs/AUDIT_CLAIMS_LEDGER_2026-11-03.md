# Audit of the claims ledger (2026-11-03)

Independent audit of `CLAIMS_LEDGER_2026-11-03.md` (CL) against `VERDICT_TABLE.md` (VT), `results_preliminary.md` (RP),
`registered_experiments.md` (RE), the archive `scCircadianCancer/docs/RESULTS_preliminary.md` (AR) and the output CSVs.
Read only; no analysis run beyond reading CSVs. Paths are relative to `paper/exploratory/` unless stated. "L" = RP line.
About 150 ledger numbers were checked against CSVs (two Sonnet agents, key values re-read by the lead from the CSVs).

## 1. Summary (ranked by impact on the paper)

1. **The headline contradicts row 4.2.** Brain-dead organ donors lose the order of organ phase relations (rho 0.161 vs
   0.678, A132b post hoc, L258-264). By CL rule 3 this is body-wide internal desynchrony, which the headline says does not
   occur ("not as body-wide internal desynchrony"). RP still calls the brain-death reading post hoc after A136 (L2741).
2. **"Chronic disease" clause: A153 is omitted.** RP L1877-1880: "the paper cannot state 'coupling preserved in disease'
   as a tested claim from these cohorts". The clause rests on A164 (T2D, obesity, GTEx) plus 14 TOST rows (A180).
3. **"Development" clause: no human multi-organ support.** A181: fetal organs give "no evidence for or against a shared
   clock" (L2619) (omitted). A190's evidence was withdrawn post hoc: "not evidence of a coupled adrenal clock" (L2898),
   yet CL 3.2 cites it at status B. Only the n 16 mouse embryo (gate failed, p_acrophase 0.083) and one fetal organ remain.
4. **Row 6.2 "tumours stay entrained" overstates.** The same amendment has H138.1 FAIL (malignant host coupling 0.076,
   below threshold; L772-776). A149 single-cell pairing is not clock-specific (L1001-1006). The supported wording is a
   bulk, clock-specific tumour-normal correlation.
5. **Row 6.3 status A is generous for the proliferation arm.** Omitted: A178 V178.1/V178.2 FAIL (SCAN-B proliferation
   rests on ARNTL/NPAS2/PER3; gastric specificity rests on GSE15459 alone, L2510-2527). Also omitted: METABRIC
   proliferation H166.3 FAIL plus its erratum (L2245), A173 H173.2 FAIL, and the A170 ER− single-cell result, which is
   not clock-specific.
6. **Row 6.1 and the "local" claim: A115 is omitted.** Uninvolved psoriatic skin already runs late against healthy skin,
   −1.81 h [−3.81, −0.69], p 0.0084 (descriptive, L741). That is a person-level component, not a lesion-local break.
   NR1D1 carries about a third of the delay (L1347).
7. **Rows 4.1-4.2 omit the drop in coupling and the registered reading.** Excess C falls from fast deaths by 0.159
   (ventilated) and 0.214 (Hardy 4) (L105-106). The registered reading of A109 is "share a donor state rather than
   running clocks" (L126-128). A109b Tabula Sapiens did not replicate (L154-163).
8. **Gate false positives are omitted next to the rows that rest on gate passes.** A145 H145.3: in no-clock simulations
   the gate passes 4 of 30 and false-claims 0.127 vs CoFE 0.06 (L898-899). A154 V154.1 fails narrowly. Rows 1.3, 3.1 and
   7.2 rely on gate passes.
   **Update 2026-11-07:** the claim rule (gate + matched pseudo-clock companion) makes 0 of 150 false claims in the same
   no-clock arm (A209, A209b; gate alone 26 and 13 of 150) and keeps 86 of 90 true couplings >= 0.5. State the gate-alone
   rate next to the companion rate; rows without a companion (A174) are checked against surgery windows.
9. **Row 7.1 says the reverse of its source, and mislabels A124.** "Blood time signal follows cell composition" vs RP:
   composition removes 6 % and the explanation is "only weakly supported" (L3194-3195). A124 (night shift) is a muscle
   study (GSE106800, 11 subjects, L1263), not blood: its null is a power limit (MDE 7.3 h), not a blood limit.
   (Added 2026-11-03 after the first version.)
10. **Row 6.4: a contrary descriptive result is omitted and the test is uninformative.** NeoPalAna baseline → C1D1
    (anastrozole) coherence rose +0.316, p 0.037 (L3165). H196.2 is uninformative (no cross-sectional relation within
    ER+). The "[−0.10, +0.14]" interval is not a computed CI.
11. **Number errors.** Row 6.1 "A54 discovery −1.43 h" is the post hoc 15-series value; A54 discovery is −1.50 h. Row
    1.4: baboon 0.455 h → 0.452 h, and R 0.987 comes from a different readout. Row 6.5: A192 p 0.12 → 0.11. Row 8.1: seed
    concordance 0.43 (whole day) / 0.45 (working hours).
12. **Status higher than evidence.** Row 1.3 is at A on seen data (V25.5 first look; A200 on the same 42 patients).
    A136 H136.3 is described as a registered confirmation but uses the same donors as A132b. Row 3.2 (A190) should be C.
    Row 6.3 gastric should be B.
13. **(Resolved for 5.1 by A202, 2026-11-04: T2D interval re-read as [−0.037, +0.076]; ledger row updated.)** **Rows 1.5 and 5.1 use bootstrap intervals**, while "Do not claim" excludes bootstrap CIs for coupling (A154,
    A162, A163). The equivalence of A164 rests on a donor-bootstrap 90 % interval.

## 2. Coverage

Used = cited as evidence; DNC = "Do not claim"; omitted. "Matters" says whether an omission supports (S), qualifies (Q)
or contradicts (X) a ledger row, or would be a claim of its own (N).

### 2a. Registered tests from A109 (VT rows)

| result | what it shows | in ledger | matters? |
|---|---|---|---|
| A109 H109.1, sec, -R | inter-organ C above pseudo and responsive nulls in all death groups | used 4.1 | excess falls fast→vent by 0.159, fast→Hardy 4 by 0.214 (L105): Q 4.1 |
| A109 H109.2 | vent offsets not preserved (rho 0.261); pulled to 0/12 h | used 4.2 | registered reading "donor state rather than running clocks" (L126) not stated: Q 4.1 |
| A109 H109.2 sec Hardy 4 | offsets preserved rho 0.784 | used 4.2 | (secondary, should be labelled so in the sentence) |
| A109b C, offsets (Tabula Sapiens) | no replication (+0.009, MDE 0.104); offsets uninformative | omitted | Q 4.1 / 2.1: the only independent replication failed (underpowered) |
| A109c H109c.1, H109c.2 | Hardy 4 wall time R 0.097 (80 % MDE 0.208); residual rho 0.750 | used 4.2 | "lose wall time" needs the power limit: entrainment ≤ 0.15 not excluded (L184-187) |
| H109.1-R / H138.4-R | responsive-gene companions pass | used 4.1, 6.2 | — |
| A110 H110.1 | infeasible | omitted | no |
| A110 H110.2a | strong clocks carry coupling (rho 0.498) | omitted | S 2.1 |
| A110 H110.2b | retained coupling vs stromal fraction, opposite sign | omitted | minor Q (stromal backbone reading) |
| A110 H110.2c, H110.4 | descriptive; hypothalamus +0.087 vs other regions −0.121 | omitted | H110.4 contrasts with A143 H143.3 FAIL; no claim |
| A110 H110.3 | temperature excluded (0.017 < 0.05) | used 4.3 | — |
| A110 H110.5 | feeding: 2 pairs, weak exclusion | used 2.3 limit | — |
| A132 H132.1, H132.2 | brain-dead vs ventilated not separated (MDE 0.54 h); both pulled to 0/12 h | omitted | Q 4.2: the registered tests of brain death failed; only the post hoc split separates |
| A132 D132.3 | dropped | omitted | no |
| A133 V133.1 | circular-part estimator fails | omitted | Q 4.2: post hoc sims say the brain-dead loss is "not a larger displacement of the same kind" (L315-322) |
| A147 V147.1/2 | order readout validated, shrink ratio biased; overall FAIL | omitted | Q 4.2 (pull-vs-order readout is descriptive) |
| A111 H111.1 Xin / Manella | feeding shift tracked (FL 0.988, min p) / fail by design | used 1.5 (Xin) | Manella fail by design should be named with Xin |
| A111 H111.2 | daytime feeding does not decouple the liver (WAT, lung drop as much) | omitted | Q 1.5: a known perturbation the readout did not localise |
| A111/A129 H111.3 | Bmal1 KO liver d* −0.287, specific | used 1.5 | — |
| A112 H112.1 | wrong sign | omitted | no |
| A113 H113.1 adipose / muscle | pass / fail as registered | used 1.4 | — |
| A130 H130.1 | internal minus wall time vs age −0.339 h/decade, p 0.054 | omitted | Q 2.4 (age phase drift not excluded) |
| A131 | not scored | omitted | no |
| A123 H123.1-3 | GTEx working-hours map vs truth (i) | used 1.1 | — |
| A123b H123.4 | orderers fail in working hours | used 1.2 | — |
| A123c H123.5 / H123.6 | offsets 0.679 FAIL; amplitudes 0.965 PASS | used 1.1 | H123.6 (amplitudes pass vs disjoint truth) not cited: S 1.1 |
| A134 H134.1, H134.2 | copy number small effect, CI includes 0; MYC/hypoxia not supported | omitted | Q 6.3 "composition" objection |
| A138 H138.1 | malignant host coupling FAIL (0.076; design change) | omitted | **X 6.2** |
| A138 H138.2/3/6/7 | not scored (gate) | omitted | no |
| A138 H138.4 | tumour-normal pairing 0.365, clock-specific | used 6.2 | — |
| A138 H138.5 | single-cell malignant vs normal not detected | omitted | Q 6.2 |
| A138 H138.8 | purity: no dependence detected (+ erratum L942) | omitted | Q 6.2 |
| A139 V139.1/2 | per-gene refit and flags fail; conserved waveform stays an assumption | omitted | Q section 1 (method assumption) |
| A140 Q1 gate, H140.1-3 | CYCLOPS fails simulator gate (0.060; core seed 0.058) | used 8.2 | core-seed result should be cited |
| A140 Q2 (i), P140.1-3 | immune gate fails, not estimable | omitted | no |
| A115 H115.1 (D115.1) | uninvolved psoriatic skin late vs healthy, −1.81 h | omitted | **X/Q 6.1 and headline "local"** |
| A116 | withdrawn | omitted | no |
| A135 | no brain reference | omitted | Q brain-disease statements (anchor dependence) |
| A142 H142.1 / H142.2 | brain split; digestive departure on r* only | used 2.2, 2.3 | — |
| A143 H143.1 / H143.2 / H143.3 | quality-adjusted 0.627 / batch excess 0.039 in brain and 0.040 in periphery (p 0.0035) / hypothalamus FAIL | used 2.2 | "no batch effect" overstated: a non-brain-specific batch effect exists |
| A144 H144.1 | cross-half rho 0.690 | used 2.1 | A144 is the post hoc direction of failed A104 (VT 221): label |
| A145 H145.1-3 | no registered advantage; **gate false-claim 0.127 vs CoFE 0.06** | used 1.2 limit | **Q 1.2-1.3, 3.1, 7.2** (gate false positives) |
| A146 H146.1-5 | recorded-time pass, internal-time fails, pan-cancer fails; GC check p 0.054 | used 6.6 | H146.4 fails in predicted direction: Q 6.6 (GC not excluded) |
| A121 H121.1/2 | onset within late window not shown; phase uninformative | used 3.1 (mis-cited) | CL attributes "earlier windows show no clock-shaped state" to A121; that is A108 descriptive gates |
| A122 H122.1/2 | immune cue-following not cue-specific | DNC | — |
| A148 | stopped at acceptance | omitted | no |
| A149 H149.1-3 | single-cell pairing not clock-specific; stroma pairing without control | omitted | **X 6.2** ("entrained") |
| A150 H150.1 | atlas coupling clock-specific overall; myeloid/lymphocyte/mural not | used 2.1 | family results needed for any cell-type claim |
| A151 H151.1/2 | fibroblast pairing 0.860; specificity uninterpretable | omitted | no (FM uses it; ledger rightly does not) |
| A152 H152.1 | residency does not carry lymphocyte coupling | omitted | no |
| A124, A156 | night shift, forced desynchrony: not detected (MDE 7.3, 13.1 h) | used 7.1, DNC | — |
| A154 V154.1-3 | gate slightly anti-conservative; c coverage fails, c biased low | used DNC | V154.1 Q 1.3, 3.1 |
| A155 | LOO anchor readout | omitted | no |
| A157 V157.1 | JS on 4-h grid fails literally | omitted | no |
| A162, A163 | coverage fails; reliability correction fails | DNC | Q 1.5, 5.1 (bootstrap intervals used there) |
| A153 H153.1-3 | no disease primary estimable (MDE > 0.2 everywhere) | omitted | **Q headline "chronic disease", row 5** |
| A159 H159.1-2, D159.1 | not scorable; 3 of 31 compartments clock-specific | omitted | Q 2.1 |
| A164 H164.1-3 | T2D, obesity equivalent within ±0.15 | used 5.1 | vent obese pancreas +0.363 (descriptive, L1910) unmentioned |
| A165 H165.1/2 | 0 of 89; endothelial z2 +0.409 [0.027, 0.792] against predicted sign; myocarditis EC +1.14, ILD EC −1.68 h at raw p 0.001 | used 5.2 | Q 5.2 (leads of change omitted) |
| A167 H167.1/2 | liver fibrosis 0 of 3, hepatocyte d* −0.209 (p 0.034, Holm 0.102) leans negative | omitted | Q 5 (a disease lead in the "loss" direction) |
| A126 H126.0-3 | menstrual natural experiment fails | omitted | Q 6.3 (no physiological-proliferation analogue) |
| A127 H127.1/2, H127.3 | SCAN-B proliferation (narrowly specific) and ER− | used 6.3 | — |
| A161 H161.1 | endometrium opposite sign | omitted | Q 6.3 |
| A166 H166.1/2/3 | METABRIC ER− replicates; **proliferation H166.3 FAIL** | used 6.3 (ER− only) | Q 6.3 (erratum L2245 says the proliferation loss rests on A97, A127) |
| A168 H168.1-3 | AI does not reduce coherence | used 6.4 | — |
| A169 H169.1/2 | no progressive loss with duration | omitted | S 6.4 |
| A170 H170.1 | ER− malignant single-cell loss not clock-specific | omitted | Q 6.3 |
| A171 H171.1-4 | LUAD pass (GSE72094 not clock-specific at 200 sets); HCC fail | used 6.3 | "not HCC" = not detected (MDE 0.32-0.36): rule 1 |
| A172 H172.1/2 | not pan-cancer (MDE 0.087) | used 6.3 | — |
| A173 H173.1-4 | gastric 2 of 3 pass; **H173.2 FAIL**; ADC vs SQC p 0.058 | used 6.3 | Q 6.3 ("not squamous" is not shown) |
| A174 V174.1-3 | heart AM/PM (archive result known) | used 1.4 | label confirmation on seen data |
| A175 H175.1/2 | TCGA pairing drop clock-specific (seen); ER not estimable | used 6.3 limit | — |
| A176 H176.1 / H176.2 | swamped 7 of 7 / disordered not informative (design flaw) | used 6.3 | "not disordered" not supported (H176.2 uninformative) |
| A178 V178.1/2 | **robustness FAIL**: proliferation rests on ARNTL/NPAS2/PER3; gastric specificity on GSE15459 alone | omitted | **Q 6.3** |
| A179 H179.1 | not cell-cycle composition (0 of 5) | used 6.3 objection | tests cell-cycle conditioning, not cell-type composition |
| A180 age / disease TOST | age not equivalent; 7/8 phase, 7/10 amplitude | used 2.4, 5.2 | non-equivalent rows (skin_ad fibroblast −1.12 h; myocarditis and DCM EC amplitude up) omitted: Q 5.2 |
| A181 V181.1, H181.2 | fetal gate uninformative by design | omitted | **Q headline "development"** |
| A182 H182.1-5 | blind types: none passes, not estimable | omitted | Q 6.3 ("not squamous": not estimable) |
| A183 H183.1/2 | psoriasis delay replicates in size, not clock-specific; not specific vs eczema | used 6.1 | — |
| A184 H184.1 | cord × maternal not resolved | omitted | no (development limit) |
| A185 H185.1 | spread-matched pairing drop not resolved (p 0.064) | omitted | Q 6.3 (pairing loss partly amplitude) |
| A136 H136.1-3 | GR gate fails (not robust: p 0.021 without steroid users); H136.3 0.458 | used 4.2 | label: same donors as A132b |
| A186 H186.1 | CPTAC replicates pairing drop, not clock-specific | used 6.3 limit | — |
| A187, A187b, A188, A195 | not estimable / post hoc null | DNC | — |
| A189 H189.1-3 | fetal adrenal coherent (Holm 0.0495); others uninformative; lung gate not clock-specific | used 3.2 | — |
| A190 H190.1 + post hoc | gate passes but confounded; "not evidence of a coupled adrenal clock" | used 3.2 (B) | **X 3.2 status** |
| A191-A194, A198 | prognosis: no information beyond proliferation (A192); not prognostic (A193, A194); no interaction (A198) | used 6.5, DNC | — |
| A196 H196.1/2 | not estimable; calibration uninformative; NeoPalAna BL→C1D1 +0.316 (p 0.037) | used 6.4 | **Q 6.4** |
| A197 H197.1-4 | dedifferentiation fails | DNC | — |
| A199 H199.1/2 | CYCLOPS magnitude barely tracks clock strength | used 8.1 | — |
| A200 H200.1/2 | breast surgery clock-specific (seen data) | used 1.3 | label seen data |
| A201 H201.0/1 | aorta shift; difference not estimable | used 7.2 | — |

### 2b. Unregistered entries of RP (post hoc, descriptive, integrity, errata)

| result (L) | what it shows | in ledger | matters? |
|---|---|---|---|
| A19 errata (L10), A19 rerun (L1992) | H19.1 is own noise (ψ) over 83 compartments incl. brain; without brain −0.009 (SE 0.084); H19.2 amplitude uninformative | used 2.4 | **row 2.4 misdescribes the readout** (see 5) |
| A20/A23 errata (L30, L2591) | neuron-loss in brain disease no longer holds | omitted | no (not claimed) |
| WP0 item 1/3 (L75) | figure uses contradict archive verdicts (A104, A54 held-out, embryo "from E15") | omitted | Q 2.1 (A144 rests on A104 direction), 6.1, 3.1 |
| WP0 item 5 (L166), D4 correction (L2559) | CHIRAL age damping vanishes with death-type matching (CPM p 0.024; TMM p 0.092) | used 2.4 objection as "not re-verified" | stale: it was re-verified, with the TMM caveat |
| WP0 item 8 (L286), A141 (L529), WP0 Q3 (L559) | brain host anchor-dependent; 4 of 12 brain contrasts anchor-dependent | omitted | Q any brain-disease statement; not 2.2 (r* is anchor-free) |
| A132b + descriptive (L258, L267) | post hoc split; ventilated pull with order kept, organ donors reorder | used 4.2 | **X headline** (see 7) |
| A109c D109c.3-4 COVID (L195) | no 0/12 h push; could not have shown it | used 4.2 limit | — |
| A144 part 2 (L661) | inconclusive | omitted | no |
| A138 power step (L721), H138.1 design change (L2327) | registered rule anti-conservative; rule changed before data | omitted | Q 6.2 |
| A149 note 5 / A151 note 2 (L1066), Phase 5 bounded (L1514) | bounded rescores; A149 still not clock-specific | omitted (except A142 erratum) | Q 6.2 |
| A123c post hoc (L1167, L1204) | estimate split-half 0.873; truth halves 0.544; Whole Blood drives the gap | used 1.1 | "corrected about 0.84" uses an RP-text reliability 0.93; the CSV's truth-noise-only value is 0.815, RP L1213 gives 0.73-0.81 |
| A123b reflection (L991) | not a reflection artefact | used 1.2 | — |
| Psoriasis robustness (L1332) | 15-series core −1.426 h; 0 of 3×200 pseudo sets; NR1D1 a third | used 6.1 (as "A54") | **mislabelled** (see 3) |
| A127 post hoc (L1828) | ER− holds within PAM50 | used 6.3 indirectly | — |
| A160 (L1482) | keratinocyte delay not testable | omitted | Q 6.1 (cell of origin unknown) |
| A187b (L2787) | post hoc null | omitted | no |
| A201 post hoc (L3134) | baseline offset +0.18 h | used 7.2 | — |
| A192 post hoc RORC (L2972) | single-gene residual | used 6.5 (RORC exception) | — |
| Blood post hoc (L3182) | no blood core clock; composition removes 6 % | used 7.1 | **X 7.1 wording** |
| A181 gate sensitivity statement (L2624) | — | omitted | Q 3 |
| A147 n_ref 244 (L1250) | order criterion missed (0.755) | omitted | Q 4.2 |

### 2c. Archive results cited by the figure plan or FINDINGS_MAP

| archive | archive verdict | in ledger | matters? |
|---|---|---|---|
| A3 | primary met (−0.532, p 0.022) | used 2.4 | — |
| A19 | not supported (rerun) | used 2.4 | see 2b |
| A21 R21.4 (with V17.2 FAIL) | colon shift recovered blind (−8.95 vs −9.03 h) | omitted (FM uses it) | S 1.5 (a blind known shift; anchors set by hand) |
| A25 V25.x, T1-T3 | V25.1 pass (M3), T1, T2 pass, T3 not supported | omitted | S/Q section 1 |
| V25.5 | post hoc first look | used 1.3 | per-gene correlations only in AR L1513 (no CSV) |
| A36 baboon | JS degenerate; offsets 0.452 h (current), R 0.972 for offsets | used 1.4 | number mix (see 3) |
| A42/A103 | brain r*, sheep pass, baboon readout fails; "message 3 does not replicate" (sheep 0.33) | used 2.2 | — |
| A54 | discovery −1.50 h; held-out −0.93 h, H54.1 and H54.3 fail | used 6.1 (wrong value) | held-out failures omitted: Q 6.1 |
| A56/A57/A59/A74/A78/A79/A81b/A82 (plan Fig 5) | A56 dose and treatment specificity fail; A79 "specific to psoriasis" | omitted | A79 vs A183 H183.2 (matched 20th pct): conflict to resolve before any psoriasis-specific text |
| A85 | pass after correction; not cue-specific (A122) | DNC | — |
| A86/A87 (plan Fig 5) | malignant coherence 0.23 vs 0.78, alignment 0.27 vs 0.64 | omitted | Q 6.2 (low malignant alignment to host) |
| A91 | rhythmic genes follow internal time | omitted | S (no row) |
| A92/A93/A95/A33 (plan Fig 3) | no decay; Hardy 4 C 0.41 vs 0.60; liver and muscle lose coupling | omitted | **Q 4.1-4.2** (coupling loss without cues) |
| A94 | H94.1-3 FAIL; small phase advance 0.2-0.5 h/decade | omitted | Q 2.4 |
| A104 | H104.1 not supported, H104.2 not replicated | omitted | A144 is its post hoc direction: label 2.1 |
| A108 | H108.1 met narrowly (R 0.516 vs bar 0.5); leave-one-day-out needs all 4 days; gate p_acrophase 0.083 | used 3.1 | — |
| A5/A17/A98 fibrotic lung T | holds, conditional on procurement | DNC | — |
| A46/A47/A96/A99, A43-A45, A52/A53, A100-A107b | mostly nulls or uninformative | omitted | A43-A45 liver-wide loss with procurement confound: Q 5 |

## 3. Number mismatches

All other numbers match their CSV to the stated digits.

| ledger row | ledger value | source value | file |
|---|---|---|---|
| 6.1 | A54 discovery −1.43 h | A54 discovery −1.499 h (9 series); −1.426 h is the post hoc 15-series A78/A81b core | `paper/results/figure_5/output/psoriasis/a54_discovery_pooled.csv` (offset_h); `psoriasis_robustness_2026-10-28/output/reproduction.csv` |
| 6.1 | "0 of 600 pseudo sets" vs A54 | 0 of 3 × 200, compared with the −1.426 h core, not with A54 | `psoriasis_robustness_2026-10-28/output/pseudo_sets_summary.csv` |
| 1.4 | baboon error 0.455 h, R 0.987 | offsets readout: error 0.4525 h, R 0.972; R 0.987 is the host-phase readout (error 0.392 h) | `paper/results/figure_1/output/species/a36_readouts.csv` rows 0, 3 |
| 6.5 | A192 0.88 (p 0.12) "with G2M" | 0.8825 [0.756, 1.031], p 0.1146, adjusted for G2M, ESR1 and ERBB2; adjusted score at the 4.1th pct of matched sets | `C8b_clock_prognosis_replication/output/a192_tests.csv` (H192.2b) |
| 8.1 | seed concordance 0.43 | 0.425 (whole day), 0.447 (08-18 h) | `A199_cyclops_working_hours/output/per_tissue.csv` |
| 6.4 | "about [−0.10, +0.14]" | not stored; equals 0.0193 ± 1.96 × sd_T_alone 0.0619 (a pooled permutation-null SD); with the row's sd_T 0.0804: about [−0.14, +0.18] | `A196_proliferation_arrest/output/calibration_posthoc.csv` |
| 1.1 | corrected about 0.84 | RP text only (reliability 0.93 not in a CSV); CSV disattenuation for truth noise 0.815 | `A123_working_hours/output/a123c/diag/split_half.csv` |
| 1.3 | ARNTL +0.69, NR1D1 −0.63, NR1D2 −0.67, PER3 −0.71 | unverified: archive doc only (AR L1513), no CSV | — |
| 7.1 | "A124 night shift" listed as a blood study | A124 is skeletal muscle (GSE106800), MDE 7.3 h | RP L1263; `A124_night_shift/output/mde.csv` |
| 8.1 | tumour seed concordance ≥ 0.976 | unverified beyond RP L3075 | — |

## 4. Label and status problems

- **1.3 (A):** V25.5 is a post hoc first look; A200 is "seen data: a robustness test, not a replication" (L3104). With
  one cohort and dependence on NR1D1, B is the ceiling until a fresh cohort.
- **1.4:** A174 is a "confirmation on current code; archive result known" (L2400), i.e. seen data. Its truth is a 4-h
  window per patient, not a recorded time.
- **1.5:** "Brain departure" and "psoriasis delay" are findings, not positive controls with a known answer. The known
  positives are Bmal1 KO, Xin feeding (minimum attainable p) and the A21 colon shift (archive; anchors set by hand).
- **2.1:** A144 is the post hoc reversal of A104, whose registered prediction (rho < 0) failed (VT 221; L83).
- **2.2:** "A143 no batch effect" should read "no brain-specific batch effect". A batch excess of 0.040 is present in
  periphery (p 0.0035, L690-694).
- **2.4:** the objection cell says the death-type reanalysis is "not re-verified here". It was re-verified (L166-178,
  L2559), with the TMM caveat (p 0.092).
- **3.1 (A/B):** the gate failed (p_acrophase 0.083, passed False) and the bar R > 0.5 was met at 0.516. Leave-one-day-out
  loses significance in 2 of 4 (AR L3035-3045). Archive, n 16: B. "Earlier windows show no clock-shaped state" is a
  failed-gate negative (rule: uninformative) and is A108 descriptive, not A121.
- **3.2:** A190 → C, or drop it (post hoc contamination, L2887-2900). A189 → B (one organ, 7 fetuses, Holm 0.0495, "fragile",
  L2859).
- **4.2:** A132b is post hoc. A136 H136.3 re-tests the same 79 + 79 donors without PER1/PER2, so it is a registered
  sensitivity on seen data, not a confirmation (RP: "still post hoc", L2741). Hardy 4 preservation is a registered
  secondary. Status: B for Hardy 4; C for the brain-death split.
- **5.1 (A):** the intervals are donor-bootstrap 90 % intervals of d*. Coverage of bootstrap intervals for c fails (A154,
  A162, A163); for d* coverage was not validated. Either state the intervals as descriptive and add the permutation p and
  MDE (0.083, 0.130), or show d* coverage.
- **6.2 (A):** A only for "bulk tumour-normal clock-state correlation, clock-specific (TCGA)". "Entrained" is not supported
  (H138.1 FAIL, A149 not clock-specific).
- **6.3:** breast ER− A (robust in A178). Breast proliferation B (A178 V178.1 FAIL; METABRIC H166.3 FAIL). LUAD A
  (GSE68465; GSE72094 clock-specific only at 1,000 sets). Gastric B (one clock-specific cohort; H173.2 FAIL). Pairing
  B (TCGA seen; CPTAC not clock-specific).
- **6.4 (B):** C. H196.1 not estimable, H196.2 uninformative, A168 secondary only; one contrary descriptive rise.
- **6.5:** label which cohorts were seen for outcome (METABRIC in A191; SCAN-B and METABRIC in A193, A194, A198), and that
  A194's "aligned" readout was chosen after seeing its distribution (before outcomes).
- **7.1:** "post hoc 2026-11-03" was suggested by the user's prior observation (L3184). The ledger's wording reverses the
  reading.
- **8.2:** the A140 gate is a failed simulator gate whose registered consequence is "H140.1-3 not scored". The CYCLOPS
  statement is a descriptive reading of the gate (concordance 0.967, L604), not the number 0.060.

## 5. Wording problems and proposed text

| row | problem | proposed sentence |
|---|---|---|
| 0 | contradicts 4.2 (brain death = body-wide reordering); "development" and "chronic disease" clauses lack A/B support (A181, A153) | "Across adult organs, internal time is shared and survives the loss of daily cues and common metabolic disease within stated bounds; where we detect clock breaks, they are confined to one organ, lesion or tumour type, with one body-wide exception, brain death." (see 7) |
| 1.1 | fine; add the amplitude pass vs the disjoint truth (0.965) and quote 0.73-0.84 as the descriptive corrected range | "…recovers the organ phase map (r 0.945 vs whole-day timing; 0.679 vs a donor-disjoint truth, below the registered 0.7, with that truth's split-half reliability 0.54) and the relative amplitudes (0.965 vs the disjoint truth)." |
| 1.3 | status, scope | "In one cohort of 42 breast resections, internal time follows the recorded surgery hour, clock-specifically (Spearman 0.558; all 200 matched gene sets below 0.215), in a registered rerun on data seen before." |
| 1.5 | brain and psoriasis are not positive controls | list Bmal1 KO (A), Xin feeding (B, minimum p), A21 colon shift (archive, B); move brain and psoriasis out |
| 2.1 | A150 is cell types against their tissue host, not organs; per-compartment specificity weak (A159 3 of 31) | "Within tissues, stromal-vascular and epithelial cell types share the host clock state beyond matched non-clock genes (mean c 0.393 vs 0.242); myeloid, lymphocyte and mural couplings are not distinguishable from matched genes." |
| 2.2 | batch wording | "…A143: quality adjustment leaves r* 0.627; no brain-specific batch effect above about 0.07." |
| 2.4 | A19 H19.1 is own noise (ψ) across all compartments, not amplitude; "no peripheral change" scope | "Prefrontal clock amplitude is lower above 60 (archive). In peripheral compartments we detect no age effect on clock noise (−0.009 per decade, SE 0.084) or relative amplitude; effects up to about 0.15 per decade are not excluded." |
| 3.1 | causal "organised"; failed gates read as absence | "In the late mouse embryo (E15-E18.75, 16 embryos), clock-gene states across cell types follow the time of day of collection (R 0.52, p 0.016; matched gene sets 0.19-0.33); whether the clock is fetal or maternally driven, and when it starts, is not resolved." |
| 3.2 | A190 withdrawn | "The human fetal adrenal at 13-17 weeks carries clock-gene phase relations beyond matched gene sets (R 0.734, 7 fetuses, Holm 0.0495); other fetal organs could not be tested at this design." |
| 4.1 | omits drop and registered reading | "Without daily cues, organs keep a clock-gene-shaped shared state above matched and responsive gene sets, reduced from fast deaths (excess −0.16 to −0.21), and mostly along a non-temporal, glucocorticoid-aligned axis." |
| 4.2 | power of the wall-time null; post hoc split; brain death is desynchrony by rule 3 | "After slow deaths (Hardy 4), wall time is not detected (R 0.097; R ≥ 0.21 excluded at 80 % power) while the order of organ phase relations is kept (rho 0.784, registered secondary). In a post hoc split, brain-dead organ donors lose this order (0.161 vs 0.678): the one body-wide internal desynchrony we observe." |
| 5.1 | "inter-organ coupling" vs d* (compartment-to-host); bootstrap intervals | "In matched GTEx donors, T2D and obesity change the coupling of metabolic organs to the donor clock by less than ±0.15 (d* +0.019 [−0.021, 0.060]; +0.004 [−0.081, 0.089], donor bootstrap)." |
| 5.2 | omits non-equivalent rows | "Where testable (8 phase and 10 amplitude rows of 32 each), most disease rows are equivalent within ±1.5 h and 1.5-fold; one phase row (atopic non-lesional fibroblasts, −1.12 h) and three amplitude rows (heart endothelium up in DCM and myocarditis; atopic endothelium down) are not." |
| 6.1 | A54 number; "fresh"; omits A115 and held-out failures | "Psoriatic lesions run their clock genes about 1.4 h late against the same person's uninvolved skin (A54 discovery −1.50 h; held-out −0.93 h; blind replication −1.45 h, Holm 0.0002, at the specificity line); uninvolved skin itself appears late against healthy skin (descriptive, −1.8 h)." |
| 6.2 | "entrained" | "Tumour clock-gene states correlate with the same patient's normal tissue beyond other patients and matched genes (TCGA bulk, 0.365 vs 0.122; pseudo q95 0.261); coupling of malignant cells to the host is not detected (below about 0.2)." |
| 6.3 | "not squamous, not HCC" = not detected; "not disordered" not tested; proliferation arm weaker | "…ER-negative breast tumours (robust to every single gene) and, less robustly, proliferative breast, lung and gastric adenocarcinomas lose clock coherence; not detected pan-cancer (MDE 0.087), in HCC (MDE 0.32-0.36) or in squamous cohorts (not estimable)." Drop "not disordered"; keep "swamped". |
| 6.4 | not supported | "Coherence did not measurably rise during two weeks of endocrine or CDK4/6 treatment (pooled change +0.019); a rise of the size of the cross-sectional proliferation contrast (~0.3) is excluded, smaller rises are not, and one descriptive series rose (+0.32)." (C) |
| 6.5 | "reproduce in direction" holds only for Cadenas genes; TimeTeller and CYCLOPS are not proliferation proxies | "The published clock-gene association (Cadenas) reproduces in direction and is largely carried by proliferation (adjusted HR 0.88 [0.76, 1.03]; RORC excepted); TimeTeller and CYCLOPS clock-strength scores are not prognostic (HR 0.966, 0.920; intervals exclude 1.06 and 1.11)." |
| 6.6 | fine (B); add GC check fails in predicted direction (p 0.054) | — |
| 7.1 | reverses source | "Human misalignment studies in blood cannot address tissue clocks: blood carries no coherent core-clock time signal (clock genes at the 52nd gene percentile; phase concordance 0.43 vs 0.85-0.94 in solid tissues); what little time signal it has leans on composition markers (post hoc)." |
| 7.2 | fine; state the gate did not fail in the Bmal1 KO (clock share halves), so the gate is not a clean negative control there (L3131) | — |
| 8.1 | seed concordance range | "…fits unstable across seeds (minimum Spearman 0.43-0.45)…" |

Undersold:
- **1.1:** H123.6 (amplitudes 0.965 against the disjoint truth) passes but is not cited. It is the cleanest registered
  pass of the method on an independent truth.
- **1.2:** TimeTeller, which is supervised, reaches 0.754 in working hours and the timed oracle 0.738; clocksync 0.945
  (A123b, L979). Stronger than comparing only with untimed tools.
- **2.2:** H142.1 also holds without PER1/PER2 (0.585) and on the uncorrected r (0.575 vs reference q05 0.802, p 0.001;
  L2318). The brain split does not depend on the reliability correction.

## 6. Missing objections

| section | objection | answered by recorded results? |
|---|---|---|
| 0 | "Brain-dead donors show body-wide desynchrony, so your headline is false." | Not answered; the only answer is to include it (section 7). |
| 1 | "Your gate passes shared non-clock structure (A145 H145.3 0.127 false-claim vs CoFE 0.06; A154 V154.1)." | Partly: A154 N1 map (factor raises pass rate to ≤ 0.10-0.14); matched pseudo-clock companions for most rows. Not for A174, V25.5 gates without companions. **Update 2026-11-07:** answered for the claim rule: gate + companion 0 of 150 false claims (A209, A209b); V25.5 now has its companion (A200); A174 is checked against surgery windows. |
| 1 | "A conserved waveform is assumed across tissues and people." | Not answered: A139 failed (L852-854). |
| 1 | "Post-mortem GTEx is not physiology (agonal state, ischaemia)." | Partly: A143 quality adjustment, A109 procurement checks, living cohorts (A113, A174, A200). |
| 2 | "Coupling c is biased low and depends on host spread; comparisons are not like-for-like." | Partly: A154/A162 document the bias; d* used for contrasts; raw c should not be compared. |
| 2 | "The brain departure is a sampling or anchor artefact." | r* is anchor-free (A142); WP0-8 shows the brain host is anchor-dependent for contrasts, which a reviewer may extend to 2.2. |
| 3 | "Embryo signal is maternal entrainment, not a fetal clock." | Open (stated). |
| 4 | "Hardy 4 and ventilated deaths are not loss of daily cues but ICU care with drugs, lights and feeding." | Partly: A110 temperature, A136 GR (gate fails, not robust). Feeding not tested well. |
| 4 | "Offset preservation from the same GTEx tissues is circular (fast-death reference vs subsets of the same tissue panel)." | Partly: A109c residual rho 0.750; donor groups are disjoint. |
| 5 | "Equivalence in GTEx (post-mortem, mostly treated T2D) says nothing about clock disruption in untreated metabolic disease." | Not answered (no living timed metabolic cohort). |
| 5 | "Single-cell disease atlases cannot test anything (A153)." | Answered by stating it; it limits the headline. |
| 6 | "Coherence loss in proliferative tumours is purity or stromal composition." | Partly: A179 is cell-cycle conditioning; A134 copy number small; A138 H138.8 (purity) failed as a test. No direct purity-stratified coherence test cited. |
| 6 | "Psoriatic delay is a keratinocyte-differentiation (transit) effect, not a clock break." | Not answered (A160 not testable, A148 stopped; A115 shows co-occurrence with differentiation, g 0.47). |
| 6 | "Prognosis cohorts were seen; multiple readouts were tried." | Partly: A192 used fresh cohorts; others not. |
| 7 | "Aorta cell types start in phase, so the test cannot show desynchrony." | Stated in the limit (post hoc +0.18 h). |
| 8 | "Your CYCLOPS refit is not their pipeline (seed list, normal eigenspace)." | Stated in 6.5; 8.1 relies on GTEx fits with 58-98 samples, which the authors would call too small. |

## 7. Headline assessment

- **Supported by A/B rows:** adult organs share an internal time readable in hospital hours (1.1-1.3, 2.1); the brain
  runs apart (2.2); metabolic disease leaves coupling within ±0.15 (5.1, single matched GTEx design); local breaks in the
  psoriatic lesion (6.1, not clock-specific in replication) and in ER−/proliferative adenocarcinomas (6.3).
- **Not supported as worded:**
  - "Development": A189 is one organ (B), A108 is n 16 with a failed gate (B), A190 is withdrawn and A181 is
    uninformative. "From late gestation" would be a hypothesis, not a result.
  - "Ageing": no A/B row shows a shared time preserved with age. Row 2.4 is about amplitude and noise, and age is not
    equivalent (A180).
  - "Loss of daily cues": holds for Hardy 4 (registered secondary); the registered primary failed in ventilated donors.
  - "Chronic disease": only T2D and obesity (A164); A153 says disease preservation cannot be a tested claim from the
    atlases.
  - "Not as body-wide internal desynchrony": contradicted by brain-dead donors (4.2), a post hoc split of a registered
    null (A132).
- **Rows it actually depends on:** 1.1, 1.3, 2.1, 2.2, 4.2 (Hardy 4 part), 5.1, 6.1, 6.3. Of these, only 1.1 (partly),
  2.2, 5.1 and the ER− part of 6.3 are A after this audit.
- **The single missing result that would most change it:** a multi-organ deliberate misalignment with timed or
  recoverable baseline offsets (a shift-work or jet-lag series with several organs). Preserved offsets would turn the
  headline from "no desynchrony where nobody expected it" into a test against the setting where the field expects it.
  A break would make brain death and shifts the two body-wide exceptions. The second most influential result is a
  fresh, registered replication of the brain-death reordering in an independent organ-donor cohort, since it is the one
  observation that currently contradicts the headline.
