# Results section outline (2026-11-23, built slide by slide with the user)

Prompt: `docs/prompts/presentation_results_2026-11-23.md`. Framing: "loudness, not time" (`docs/ABSTRACT_DRAFT_2026-11-22.md`,
`docs/PAPER_PLAN_2026-11-19.md` section 11). Slides kept from `results.md` are referenced, not repeated. Tags: blind / not
blind / descriptive / post hoc. Every number below was read from the output table named in its **Source** line.

---

# 1. Validation: the readout sees time, and sees it move

| slide | status | figure |
|---|---|---|
| 1.1 GTEx organ map from working-hours deaths | kept (`results.md` R1.1) | `figures/gtex_organ_map/` |
| 1.2 Living patients: breast, heart, adipose | kept as is (`results.md` R1.3; user, 2026-11-23) | `figures/living_patients/` |
| 1.3 Baboon atlas, 63 tissue offsets | kept (`results.md` R1.4; user, 2026-11-23) | `figures/baboon_offsets/` |
| 1.4 Timed living skin and muscle, and hospital hours | **new** (below) | `figures/living_timed/` |
| 1.5 Positive controls: time moves | to revise (next) | `figures/positive_controls/` |

## 1.4 In timed living skin and muscle, the untimed readouts match the timed ones, and survive hospital hours

**Bullets**
1. Offset between compartments of one person: skin dermis minus epidermis, clocksync -1.58 h (90 % -1.81 to -1.34; 148
   untimed women) against timed -1.98 h (-3.15 to -0.95; 20 men of the same study) (A251 H251.4, pass). Exposed minus
   protected skin site: 1.33 h against 1.57 h timed (A342).
2. Pre-mRNA lead over mRNA: skin protected 1.19 h against 1.18 h timed, exposed 1.41 against 1.48 (A342 H342.1-2, pass);
   muscle 1.42 against 1.75 (A337 Perrin H337.1-2, pass on 14 genes).
3. Hospital hours (same skin, 06, 12, 18 h only): the lead stays within the bound at both sites (untimed minus truth
   protected +0.10 h, -0.44 to 0.63; exposed -0.09, -0.57 to 0.53; H344.2 pass); the offset keeps its value (1.37 against
   1.57 h) but its interval is too wide to decide (difference -0.20, -1.68 to 0.88; H344.1 undecidable). With two
   daytime points (06, 12 h) the protected lead drifts +0.85 h (0.14 to 1.64; H344.3 fail).
4. Hospital hours in muscle (08, 12, 16 h only): the lead is recovered, untimed minus truth -0.01 h (-0.50 to 0.63;
   A353 H353.1 pass); also on 08-20 h (-0.35, -0.57 to -0.10) and on two points, 08, 12 h (+0.16, -0.55 to 0.89; H353.2
   pass).

**Figure:** `figures/living_timed/living_timed.png` (`scripts/living_timed.py`; wide). A. offset between compartments,
clocksync (blue circles) against timed (black diamonds). B. pre-mRNA lead, same encoding. Filled: 14 clock genes; open:
core loop. C. skin (A344) and muscle (A353) cut to subsets of clock times (full day, three day times, two day times):
untimed minus the full-cycle timed truth, skin lead at each site, muscle lead and the skin offset; grey bands, registered bounds (lead +-1.5 h, offset +-1 h). 90 % intervals throughout.

**Tags:** A251 not blind (registered, one cohort); A342 registered, scored on direct introns after the erratum; A337
Perrin blind for the unspliced counts; A344 and A353 not blind for the full cycle, blind for every subset value.

**Caveats (speaker notes)**
- The exposed-protected pairing is not clock-specific against companions (H342.3 fails on that condition only;
  coherence 0.467 against matched q95 0.549): the offsets agree, the pairing could be shared visit-level variation.
- Dermis-epidermis compares different people: untimed women against timed men.
- Core loop: muscle lead under-read by about 1.3 h (1.63 against 2.97 h; H337.2 core FAIL). Skin hospital hours on the
  core: pass on 06-18 h (leads -0.13, -0.48 h); on 06-12 h the exposed lead also fails (-1.16 h).
- Reading each visit as a time of day degrades as the arc shrinks (median error 2.03 h full day, 2.59 h 06-18 h, 3.27 h
  06-12 h; night pair 3.79 h; descriptive): the method reports relationships, not sample times.
- The H344.1 interval is wide largely because the timed truth itself is (full set -1.00 to 0.51).
- Muscle passes on two daytime points where skin does not: six full-day times and a well-entrained protocol (hourly
  isocaloric meals, controlled light), the most favourable design. Its per-biopsy time error also grows as the arc
  shrinks (1.42 h full day, 2.51 h on 08-16 h, 4.00 h on 08-12 h; A353, descriptive). Muscle has one compartment: no
  offset there.

**Source:** `exploratory/A251_skin_sites/output/{score_offsets,score_descriptive}.csv`;
`exploratory/A342_skin_timed/output/score_direct/{pairing,leads}.csv`;
`exploratory/A337_timed_lead/output/perrin/verdicts.csv`; `exploratory/A342_skin_timed/output/score_a344/subsets.csv`;
`exploratory/A337_timed_lead/output/a353/subsets.csv`; `docs/results_preliminary.md` A251, A337 Perrin, A342, A344,
A353.
