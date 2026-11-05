# aiello2020circadian
**Ref:** Aiello et al., Sci Adv 6:eaaz4530 (2020). DOI 10.1126/sciadv.aaz4530. Type: biology (mouse). Lab: Golombek / Paladino (Quilmes), Finkielstein (Virginia Tech).
**Source:** full text

## Question and design
Does chronic jet lag (CJL; 6-h advance of LD12:12 every 2 days) remodel the tumour-immune microenvironment and speed up B16F0
melanoma (non-metastatic) in C57BL/6J mice? Tumour growth, latency and survival n = 14 per arm; tumour and spleen
macrophages by flow at ZT3, ZT9, ZT15, ZT21; cytokines (ELISA) at ZT6 and ZT18; liver clock and cell-cycle genes by qPCR at
8 times (ZT0-ZT21). n = 3-4 per condition and time point. 3 weeks of CJL before tumour challenge.
**Sampling:** timed, full 24 h but sparse (4 times for flow, 8 for liver), tumour-free and tumour-bearing mice in LD and CJL.

## Method (if any)
Cosinor (24 h) and Kruskal-Wallis per tissue; day (ZT3+ZT9) vs night (ZT15+ZT21) grouping for immune cells.

## Key results
- Growth: CJL tumours larger (day 33: 1607.5 vs 432.7 mm3, P < 0.01), survival P = 0.0153, latency 15 vs 25.8 days (P = 0.037).
- Tumour macrophages: under LD "M1 and M2 macrophage levels were in antiphase, with an increase in the level of M1 cells during the night"; in CJL "the levels of M2 macrophages lost their day-night differences and remained constant throughout the whole cycle" (M1 pattern inverted, LD Night vs CJL Night P = 0.0004).
- Tumour cytokines IL-1, IL-6 daily cycle in LD, flat in CJL; TNF-alpha antiphase between LD and CJL.
- Host liver clock: LD Bmal1 cosinor P = 0.006, acrophase 23.59 h; Cry1 P = 0.024, acrophase 21.23 h; under CJL cosinor not significant, "CJL might cause a phase advance in the rhythmic expression of clock genes" (peak at ZT3) and "the liver circadian clock is altered but not completely abrogated".
- Tumour clock: "The expression of Bmal1 and Cry1 genes did not show a rhythmic pattern under either an LD or CJL schedule" (fig. S3).
- Specificity: the tumour does not change the liver clock, "the presence of the tumor does not play a substantial role in such alterations" (it changes cell-cycle genes), and in spleen "the sole presence of the tumor was sufficient to influence the polarization of macrophages under LD conditions compared to sham animals".

## Not addressed / limits
- Tiny n per time point (3-4), cosinor on 8 points, one tumour model; the arrhythmic tumour clock is unresolved (asynchronous cells vs real loss), as the authors note.
- Host liver and tumour are compared only qualitatively (rhythmic vs not); no phase or amplitude estimate with confidence for the tumour, no test of a tumour-host phase relation.
- CJL changes light schedule for host and tumour alike; no tumour-specific clock manipulation.

## Link to us
- **Intro claim(s):** CJL speeds tumour growth and removes or inverts day-night immune patterns; our refs row is SUPPORTED. It also supports the claim that tumours may be arrhythmic while the host liver keeps (damped) rhythms.
- **Relation to clocksync:** background; a timed tumour plus host-organ design that could serve as a mouse positive control for tumour-vs-host coupling.

## Ideas
- Data: ask for the 8-time liver/tumour qPCR to test whether a covariance-based coupling c (no times) recovers the LD vs CJL and tumour vs liver differences seen by cosinor (not done: no mouse tumour-host test in registered experiments).
- Cite the host-clock-reprogramming-by-tumour work they build on (Masri 2016 lung adenocarcinoma/liver; Hojo 2017 breast cancer/liver); check whether in our literature set.
