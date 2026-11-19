# Registry of amendments A109 onward (2026-11-19)

219 rows in `registry_2026-11-19.tsv` (A109 to A327, lettered variants and the A214 scout included). Built by a
Sonnet subagent; lead corrections below.

## By verdict

| verdict | n |
|---|---|
| PASS | 31 |
| partial | 74 |
| FAIL | 75 |
| not estimable | 15 |
| descriptive | 9 |
| gate failed | 6 |
| not scored | 7 |
| not decidable | 2 |

## By blind status

"(inventory)" marks rows whose entry states no blind status, filled from the `blind` column of
`inventory_2026-11-18.tsv` (yes / partly / no).

| blind | n |
|---|---|
| not blind | 66 |
| blind (inventory) | 64 |
| partly blind (inventory) | 37 |
| not stated | 35 |
| blind | 6 |
| descriptive | 5 |
| partly blind | 5 |
| registered re-test | 1 |

## By figure (an ID cited in several figures counts in each)

| figure | n |
|---|---|
| Fig1 | 5 |
| Fig1 Supp | 7 |
| Fig2 | 10 |
| Fig2 Supp | 6 |
| Fig3 | 11 |
| Fig3 Supp | 8 |
| Fig4 | 10 |
| Fig4 Supp | 7 |
| Fig5 | 10 |
| Fig5 Supp | 6 |
| Fig6 | 14 |
| Fig6 Supp | 7 |
| Fig7 | 4 |
| Supplement | 5 |
| out | 129 |

## Verdict x blind (inventory fallback merged)

| verdict | blind | descriptive | not blind | not stated | partly blind | registered re-test | total |
|---|---|---|---|---|---|---|---|
| PASS | 13 | 0 | 8 | 4 | 5 | 1 | 31 |
| partial | 19 | 0 | 26 | 11 | 18 | 0 | 74 |
| FAIL | 25 | 0 | 24 | 11 | 15 | 0 | 75 |
| not estimable | 7 | 0 | 2 | 3 | 3 | 0 | 15 |
| descriptive | 2 | 5 | 2 | 0 | 0 | 0 | 9 |
| gate failed | 2 | 0 | 1 | 3 | 0 | 0 | 6 |
| not scored | 2 | 0 | 2 | 2 | 1 | 0 | 7 |
| not decidable | 0 | 0 | 1 | 1 | 0 | 0 | 2 |
| total | 70 | 5 | 66 | 35 | 42 | 1 | 219 |

## Lead corrections (2026-11-19, spot-check of 24 rows)

Dates and figure assignments matched the logs and PAPER_PLAN section 5 in all 24 rows checked
(A109, A138, A164, A183, A208, A231, A242, A244, A254, A257, A264, A266, A291, A293, A304, A306, A312, A315, A317,
A321, A323, A324, A326, A327). Three labels were wrong and are corrected in the table:
- A321 blind (the RP heading states "BLIND for the atlases"), not "not blind";
- A324 partly blind (blind for heart, kidney and VAT; the Xin liver was seen);
- A254 partial (H254.1-2 FAIL, H254.3 clock-specific at the line), not FAIL.
The blind status of 101 rows without a statement in the entry was filled from the inventory; 35 remain "not stated".
A208 and A231 sit in PAPER_PLAN section 3 (the map) but in no figure of section 5, hence "out"; they enter Figure 5a
through the two-axis table.

## Method note

Registration date is the date of the first commit that added the amendment heading to `registered_experiments.md`; scoring date is the first commit that added a heading naming the ID to `results_preliminary.md` (blank when none). Both were read in one pass over `git log -p` of each file, so multi-ID headings (A111/A112/A129, A113/A123/A130/A131, A215/A216, A217/A218, A237-A239, A157/A158) give every member the same date. The repository starts on 2026-10-23, so no date needed the heading fallback, except A294b (renumbered on 2026-11-18 from a second A294; both dates are the heading dates, 2026-11-17, marked `*`) and A214 (a scout with no registration heading; registered date blank). Scoring dates for A111, A112, A129, A132b, A133 and A159, whose results are under non-ID headings, were found with `git log -S` on the heading text. Blind status is read from the entry (heading, registration text, results heading): "NOT BLIND", "partly blind", "blind" or "registered re-test"; where the entry says none of these it is "not stated" (or "descriptive" for descriptive verdicts). Verdicts reuse the 2026-11-18 inventory (composite inventory rows mapped to their members) and were read directly from `results_preliminary.md` for the IDs the inventory lacks (A111, A112, A129, A116, A131, A132b, A133, A159, A290, A298, A316-A327) and for A144, A263, A271, A286, A305. Figures follow section 5 of `PAPER_PLAN_2026-11-19.md`; ranges count for every ID inside them.
