# Scout: GTEx donors with a cancer history (2026-11-12)

Question: can the Masri 2016 claim (a distant tumour rewires circadian output of other organs without moving their
core clock) be tested in GTEx donors with cancer, with the A233 / A235 matched designs? Aggregates only (controlled
access); no clock statistic computed. Source: `common.gtex.phenotypes()` (pht002742), 980 subjects.

| variable | meaning | cases (all Hardy classes) | fast (Hardy 1-2) | ventilated (Hardy 0) |
|---|---|---|---|---|
| MHCANCERC | current cancer diagnosis | 7 | < 5 | 5 |
| MHCANCER5 | cancer diagnosis within 5 years | 20 | 7 | 9 |
| MHCANCERNM | history of non-metastatic cancer | 69 | 24 | 25 |

Reading: not estimable for an active tumour. A233 needed >= 10 matched cases per tissue and cohort; current cancer
has 7 donors in total and the 5-year variable 7 fast and 9 ventilated donors before any tissue filter. Only the
history variable reaches the A233 rule (24 fast, 25 ventilated), and it is a past diagnosis, mostly without a tumour
at death: it tests "ever had cancer", not Masri's tumour-bearing state. GTEx eligibility also excludes metastatic
disease and recent chemotherapy. Alternatives (not scouted): the tumour-free organs of cancer patients in public
single-cell or bulk cohorts (e.g. normal liver resected with a metastasis) against healthy donors, which is a
between-cohort contrast.
