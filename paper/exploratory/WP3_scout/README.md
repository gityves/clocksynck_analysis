# WP3 scout (A111, A112, A129)

Data scout and builders only, written before any hypothesis statistic of A111, A112 or A129 is computed (2026-10-24).
`scout_wp3.py` writes `output/scout_wp3.csv` (git-ignored): one row per dataset x arm with mice per time point,
tissues, mice carrying all / at least two tissues, count type, clock genes present (of 14) and with median count >= 5.

Builders (`paper/common/datasets/<name>/build.py`, run `paper/common/run.sh <builder>` to write the standard table;
shared id -> symbol helper `paper/common/mouse_ids.py`, A48 orthologue convention = upper-cased symbol through
`clocksync.genes.canonical`): `manella2021`, `xin2021`, `greenwell2019` (+ `counts()` for the A112 gene sets),
`jouffe2022_bmal1ko`, `vipko_gse296764`. Raw files: `$SCCC_DATA/<name>/` with `PROVENANCE.md`; Xin 2021 is the
existing `$SCCC_DATA/mouse_ageing/Xin2021/` (downloaded 2026-10-16).

| dataset | accession | design | count type | tissues | same mouse across tissues | use |
| --- | --- | --- | --- | --- | --- | --- |
| Manella 2021 | GSE159135 | male C57BL/6, 3-4 months; Alb-Cre (control) or Alb-Cre Bmal1 fl/fl (BLKO, hepatocyte Bmal1 knockout); AL or light-phase feeding (DF) 30 days; ZT0-22 every 2 h; 2 mice per ZT per arm; 4 arms x 24 mice | raw 3' MARS-seq reads, ~11.6 k pre-filtered genes | liver, WAT, lung | labels `<ZT><A/B>` match across tissue files (assumed same mouse, not stated) | A111 primary, H111.3, A129 |
| Xin 2021 | GSE150380 / 150381 / 151221 / 151228 | female C57BL/6J, 9 weeks; 1 week DRF or NRF; ZT0-24 every 4 h (ZT0 and ZT24 different mice); n = 4 | TPM only | liver, heart, kidney, VAT | replicate `_k` matched across series (assumed) | A111 replication |
| Greenwell 2019 | GSE118967 | male C57BL/6J; ad libitum, night-restricted, arrhythmic, arrhythmic-unadapted; ZT2-22 every 4 h; n = 3 | raw 3' counts + TPM, 51.9 k genes | liver | single tissue | A112 |
| Jouffe 2022 | GSE150102 | whole-body Bmal1 KO x ob/ob, 20 weeks, ZT12 only | raw | liver, muscle, eWAT (1-3 per mouse) | yes (animal id) | A129 candidate: global KO, one ZT |
| VIP-Cre Bmal1 KO | GSE296764 | WT vs Bmal1-VIPKO, AL or TRF, CT0 and CT12, 3-4 mice per cell | raw | SCN, liver, heart | yes (`S<n>`) | A129 candidate: central, two times |

Facts for H111.3 (Manella): the mutant is a hepatocyte-specific Bmal1 knockout (`Alb-Cre+ Bmal1 fl/fl`, "BLKO");
controls are `Alb-Cre+` mice without the floxed allele. Floxed
exon deletion leaves the Arntl transcript detected in 3' data (liver ARNTL median count AL 35 control vs 39.5 BLKO), so
Arntl abundance is not a knockout readout. Only liver, WAT and lung were sequenced (the 6 tissues were qPCR only).

A129 candidates (public multi-tissue sets of tissue-specific clock loss or reconstitution):

| set | accession | status |
| --- | --- | --- |
| Koronowski 2019 Cell, liver-RE (Bmal1 stop-FL) | GSE117134 (WT, KO, Liver-RE, ZT0-20, n = 3, female); GSE158600 (+ AL / night feeding, male) | liver only: no second tissue, fails |
| Welz 2019 Cell, epidermis-RE | GSE115104 (LD), GSE114943 (DD): WT, KO, RE, 6 times, n = 4 | epidermis only: fails |
| Manella 2021 BLKO | GSE159135 | passes: 24 mice per genotype x feeding with liver + WAT + lung; only 2 intact companion tissues |
| Jouffe 2022 | GSE150102 | whole-body KO, one ZT, <= 6 mice with >= 2 tissues per arm: descriptive only |
| VIP-Cre Bmal1 KO | GSE296764 | central clock; 12-14 mice per genotype; 2 times only |
| GSE288746-8 (liver clock and WAT; WT, LKO, Liver-RE) | WAT only in these series | liver not public here: fails |
| tissue-specific KO, liver only (GSE171432, GSE95357, GSE115264, GSE190221, GSE284601) | liver only | fails |

No set with a tissue lacking its own clock in more than 3 tissues was found; A129 can use Manella BLKO.
