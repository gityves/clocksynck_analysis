# matsuo2003control
**Ref:** Matsuo, Yamaguchi, Mitsui, Emi, Shimoda, Okamura, Science 2003 (302:255-259). DOI not in text (not confirmed). Type: biology. Lab: Okamura (Kyoto).
**Source:** full text

## Question and design
Mouse (wild type, Cry1/Cry2-deficient, Clock/Clock) liver regeneration after two-thirds partial hepatectomy (PH) done at ZT8 or ZT0 (12:12 light-dark). Readouts: BrdU (S phase), mitotic index (3000 hepatocytes/animal), Cdc2 and Wee1 kinase activity, northern blots and a 68-gene cell-cycle microarray, Per2 immunohistochemistry on sections. **Sampling:** timed, hours after PH (to ~72 h), fixed tissue from different animals; "single cell" means immunostained hepatocytes, not live imaging. Liver-weight n = 4-6 per group.

## Method (if any)
None. Perturbation by surgery time, genotype and a retrorsine mitotic block.

## Key results
- Mitosis timing depends on the time of PH, S phase does not. "only a few cells entered the M phase within 44 hours, and a mitotic peak was reached 48 hours after PH." (PH/ZT0 versus massive entry within 40 h for PH/ZT8.)
- Cdc2 follows: "occurred 40 and 48 hours after PH, respectively, corresponding to the observed mitotic peaks".
- Of 68 cell-cycle genes only cyclin B1, cdc2 and wee1 differ markedly: "showed remarkably different expression profiles"; mRNA peaks delayed 8-12 h after PH/ZT0.
- Cry-deficient mice lose the time-of-operation effect and impair S-to-M progression: "Cry-deficient mice did not exhibit a mitotic wave that depended on the time of operation"; mitotic maximum below 4%.
- wee1 is a clock target: "wee1 mRNA levels were high at both ZT0 and ZT12" (Cry-deficient), low in Clock/Clock; "Three E-box elements were found within 1.2 kb of the mouse wee1 gene" (CLOCK-BMAL1 reporter activation, suppressed by PER and CRY).
- Clock runs in dividing cells: "more than 70% of hepatocytes" cycled on day 2 while clock genes kept oscillating; "Therefore, the circadian clockwork oscillated accurately in a single proliferating cell." With retrorsine (mitosis below 0.2%): "PER2 expression was still observed in almost all hepatocytes (99.1% of 2000 cells)" at ZT20.

## Not addressed / limits
- Shows the clock does not need division to oscillate (PER2 nuclear timing); does not test whether division resets or shifts the clock. Peak timing of one clock readout at coarse resolution.
- Mouse liver regeneration, a perturbed non-steady state; no live single-cell traces; the wee1 link is mechanistic only for G2/M entry.
- "or neural timing cues may also participate in cell-cycle regulation (12)": systemic cues not excluded.

## Link to us
- **Intro claim(s):** refs_cancer_phase_inference (clock controls cell-cycle genes; in single cells the clock runs independently of the cell cycle, i.e. unidirectional clock-to-cell-cycle control) is partly overstated: independence is shown for PER2 timing, not unidirectionality; Bieler 2014 and Droin 2019 later found the cell cycle shifts the clock in NIH3T3 and U2OS cells (cell cycle to clock dominant). Wording: clock gates G2/M via Wee1 in liver, while the clock itself persists in dividing cells.
- **Relation to clocksync:** background (clock to cell-cycle direction).

## Ideas
- In tumour data, test whether WEE1 and CCNB1 behave as clock-locked in healthy tissue and lose that in tumours (check whether done: grep A97/A127 uses G2M set as proliferation, not WEE1-clock coupling).
