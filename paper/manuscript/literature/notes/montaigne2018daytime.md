# montaigne2018daytime
**Ref:** Montaigne D et al. (Staels lab), Lancet 391:59-69, 2018. DOI 10.1016/S0140-6736(17)32132-3. Type: clinical (with human tissue transcriptomics). Lab: Staels / Lefebvre (Lille).
**Source:** full text

## Question and design
Does time of day of aortic valve replacement (AVR) change perioperative myocardial ischaemia-reperfusion injury, and by which mechanism? Single centre (Lille University Hospital), severe aortic stenosis, LVEF > 50 %. (i) Prospective observational cohort, propensity-matched: n = 596 (298 morning, 298 afternoon), 500 days follow-up. (ii) Randomised study NCT02812901: 88 patients, 44 morning / 44 afternoon, troponin T AUC. (iii) Right atrial appendage biopsies, ex vivo hypoxia-reoxygenation of trabeculae (n = 14 morning, 16 afternoon) and Agilent microarray (18 morning vs 11 afternoon patients, plus qPCR n = 28-30 per group). (iv) Mouse Langendorff, Rev-Erba KO or SR8278 at ZT8. **Sampling:** timed but coarse: two slots (morning vs afternoon, same surgical team, same day, i.e. working hours only). Exact clock hours of the slots are not given in the main text (appendix not available); one biopsy per patient, taken during cardiopulmonary bypass preparation. Each patient contributes a single time point.

## Method (if any)
Propensity score matching, Cox model, randomisation, differential expression (fold change > 1.2, p < 0.05, FDR), KEGG enrichment. No phase inference.

## Key results
- Cohort: major adverse cardiac events lower in the afternoon: "afternoon surgery group than in the morning group: hazard ratio 0·50 (95% CI" 0.32-0.77, p = 0.0021 (Summary).
- Randomised: troponin T AUC 179 ng/L (afternoon) vs 225 ng/L (morning); "estimated ratio of geometric means for afternoon to morning of 0·79 [95% CI 0·68–0·93; p=0·0045]" (Fig 2A).
- Same team did both slots: "the same surgical team did both morning and afternoon" surgery on the same day.
- Atrial biopsies: contraction recovery after hypoxia-reoxygenation better in afternoon samples (Fig 2B).
- Transcriptome: "showed that the expression of 287 genes was" regulated by time of day, the KEGG theme circadian rhythm is enriched; REV-ERBa high in the morning, BMAL1/ARNTL anti-phasic, with BMAL1 "prominently expressed during the afternoon when" REV-ERBa is low (Fig 2E).
- Mouse: Rev-Erba deletion or antagonist protects at the sleep-to-wake transition via CDKN1a/p21.

## Not addressed / limits
- Two time slots only; no full 24 h series; human phase of the clock never estimated, and morning vs afternoon contrast confounds clock phase with anaesthesia, fasting, surgery duration, sleep state.
- Inter-individual phase (chronotype) ignored; "Moreover, the mouse is a nocturnal species."; "remains correlative. Moreover, time-of-the-day probably" interferes with other systemic processes.
- Only AVR, one centre, atrial tissue (not ventricle); no cancer, no other organs.
- Sample sizes of the transcriptome are small (18 vs 11); differential expression between two groups, not rhythm fits.

## Link to us
- **Intro claim(s):** our current claim (refs_*.md): afternoon AVR associated with fewer perioperative myocardial injury events than morning (propensity-matched cohort plus trial; REV-ERB link). Supported (HR 0.50, randomised troponin ratio 0.79). Not a hospital-hours method paper: it shows that surgery-time clock phase matters and that human tissue sampled in the clinic carries a detectable clock signal (REV-ERBa vs BMAL1 antiphase between 2 daytime slots).
- **Relation to clocksync:** background; also a data / design reference (hospital-hours samples with a ~4-6 h clock-time contrast; one organ).

## Ideas
- Direct test of our hospital-hours claim: in any cohort with morning and afternoon surgery biopsies and clock genes, estimate offsets arg z per donor and see whether the morning-afternoon difference appears; the paper shows antiphase REV-ERBa / BMAL1 between slots, so a coarse expected effect size exists.
- Sensitivity check for the method: shift simulated sample times by a morning-vs-afternoon contrast and verify that donor-referenced phase offsets cancel it (constant offset), as claimed.
