# Scout: post-mortem brain (GSE45642, GSE71620) and living skin (GSE139300 + GSE139301), 2026-11-12

No clock statistic was computed. Script: `scripts/brain_skin_scout.py` (+ `scripts/brain_skin_pivot.py`), tables in `output/`.
Data root `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/`. Per-donor covariates in `output/*_samples.csv` (GSE71620 TOD is public GEO data).

## 1. Accessions and files

| dataset | accession(s) | status | folder | size |
|---|---|---|---|---|
| Li 2013 brain, 6 regions | GSE45642 (title "[control set]") | correct accession; already on disk since 2026-10-11 (not re-downloaded) | `brain_tod/GSE45642/` | 132 MB (matrix 43 MB, platform soft 91 MB) |
| Chen 2016 PFC | GSE71620 | correct; already on disk since 2026-10-11 | `brain_tod/GSE71620/` | 58 MB |
| Wu 2020 skin, 154 women | **GSE139301 (epidermis) + GSE139300 (dermis)**, both subseries of SuperSeries GSE139305, platform GPL13667 (HG-U219 GeneTitan) | new, downloaded 2026-11-12 | `skin_wu2020_gse139301/` (PROVENANCE.md with sha256) | 423 MB (two matrices 141 MB each, GPL13667 table 162 MB) |

Skin accession correction: the candidate in the review notes (Kimball 2018 / P&G ageing series) is not the one. The 154-woman three-site design with both layers is GSE139300 (dermis) and GSE139301 (epidermis), submitted by Procter and Gamble (Mills, Oblong; Sherrill) in 2019; GSE112660 (Wu 2018 PNAS, 20 men + 219 women, forearm epidermis only, 298 samples) is a separate series, not downloaded. Existing PROVENANCE.md of the brain folders carries URL, date and sha256; only the skin folder was written today. Total new downloads 0.42 GB.

## 2. Problems that matter (read first)

1. **GSE45642 contains controls only, and no covariates.** 670 chips, 86 control subjects. The 34 MDD patients are NOT in GEO (title "[control set]"; no companion MDD series found among adjacent accessions or by GEO search). Age, sex, time of death, pH, PMI, RIN, diagnosis details exist only in the PNAS SI (Table S4 per third-party description); SI download was blocked (403) again today. Without it, GSE45642 has no time axis and no MDD arm.
2. **GSE45642 regions were hybridised at different assay sites.** Region x assay site (chips): amygdala D30/I32; anterior cingulate I64/M60; cerebellum D21/I42/M16; DLPFC D84/I9/M68; hippocampus D47/M61; NAcc D46/M90. 78 of 86 subjects have chips from at least 2 assay sites. Between-region contrasts within a person are partly confounded with site and with chip type (U133A 393 vs U133Plus2 277, only the U133A probes used). 261 of 365 subject-region cells have more than one chip (replicates, to be averaged or modelled).
3. **Skin: every (layer x site) compartment was hybridised in its own batch.** Women's six compartments sit in six separate GSS batches (epidermis cheek GSS1339, buttock GSS1376, forearm GSS1403; dermis cheek GSS1340, buttock GSS1377, forearm GSS1404; `output/wu2020_women_batch_by_layer.csv`). Compartment and batch are fully confounded, so any between-compartment offset in mean level is untestable; only within-compartment across-woman covariance (each compartment against the same women's other compartments) is clean. Dermis and epidermis matrices are separate RMA runs with different global levels (matrix median 4.37 dermis vs 2.98 epidermis log2).
4. **Skin: no age in GEO** (only tissue, collection time, sex). Age 20-74 is from the paper; the Genome Med supplement could not be fetched (site returned a redirect/stub). Age would have to come from Wu 2020 additional files or the authors. Women have no collection time (NA, as the paper says); the timed arm is 20 men (hhmm 0, 600, 1200, 1800), forearm, both layers.
5. **Skin subject keys.** Women: subject number in the title (101-321) identical across all six compartments, so linking is direct (154 women, 137 with all six compartments). Men: epidermis numbers 101-120, dermis 1001-1020; the mapping 1000+n <-> n is plausible but unverified (not stated in GEO). Women and men numbers overlap numerically (17 shared integers), so the sex arm must be part of the key. Forearm epidermis is labelled "arm" in the epidermis series, "forearm" in the dermis series.
6. **GSE71620 (Chen) is not controls-versus-MDD and has two compartments only (BA11, BA47)**, no donor id in GEO. 64 of 210 subjects have no TOD (unwitnessed deaths).

## 3. Design tables

### GSE45642 (Li 2013), controls only, HG-U133A / U133Plus2 chips, custom CDF Entrez probesets
Subject key: GEO characteristic "subject id" (integer, shared across regions); checked that it recurs across regions (365 subject-region cells for 86 subjects).

| region | subjects | chips |
|---|---|---|
| amygdala | 32 | 62 |
| anterior cingulate | 70 | 124 |
| cerebellum | 51 | 79 |
| DLPFC | 83 | 161 |
| hippocampus | 63 | 108 |
| nucleus accumbens | 66 | 136 |

Regions per subject: 1 region 8, 2 regions 9, 3 regions 5, 4 regions 24, 5 regions 12, **6 regions 28 subjects**. Group x region: all control (no MDD in GEO). Chips per subject-region cell: 1 chip 104, 2 chips 227, 3 chips 28, 4 chips 2, 5 chips 4. Subjects with chips from >1 assay site: 78; per subject-region cell 2 sites in 250 cells and 3 sites in 8 cells (replicate chips of the same cell are split across sites). Time of death, age, sex: not in GEO (empty columns in `GSE45642_samples.csv`).
Files: `output/GSE45642_samples.csv`, `_design.csv`, `_subject_by_region.csv`, `_region_by_assay_site.csv`.

### GSE71620 (Chen 2016), controls only, Human Gene 1.1 ST (GPL11532)
Subject key: derived. No donor id in GEO; the tuple (age, sex, race, PMI, pH, TOD) is identical in a BA11 and a BA47 chip and yields 210 distinct subjects, each with exactly one BA11 and one BA47 chip (210 with both, 0 with one region).

| item | value |
|---|---|
| chips / subjects | 420 / 210 (BA11 210, BA47 210) |
| subjects with TOD | 146 (64 without) |
| sex | 166 M, 44 F |
| age min / median / max | 16 / 52 / 96 |
| age bins (subjects) | <=29: 26; 30-39: 17; 40-49: 41; 50-59: 73; 60-69: 37; 70-79: 10; 80+: 6 |

TOD: decimal hours, ZT with ZT0 = sunrise (per Chen 2016; GEO does not state the convention); range -5.96 to 17.40, stored mod 24 in `tod_zt_0_24`. Histogram, 146 subjects, 2 h bins of ZT (only group = control): 0-2: 2; 2-4: 9; 4-6: 17; 6-8: 21; 8-10: 20; 10-12: 17; 12-14: 14; 14-16: 16; 16-18: 10; 18-20: 10; 20-22: 6; 22-24: 4 (the late-night/early-morning bins are thin). Also in the sheet: pH, PMI, RIN, race. No batch or site fields exist in GEO.
Files: `output/GSE71620_samples.csv`, `_design.csv`, `_age_bins.csv`, `_tod_hist_2h.csv`.

### Wu 2020 skin (GSE139301 epidermis, GSE139300 dermis), HG-U219 GeneTitan, RMA
Women (untimed): 154 women, 906 samples, 6 compartments per woman = 3 sites (forearm, buttock, cheek) x 2 layers; 137 women have all six, 16 have five, 1 has four.

| compartment (women) | forearm | buttock | cheek |
|---|---|---|---|
| epidermis | 152 | 151 | 151 |
| dermis | 149 | 151 | 152 |

All 154 women are present in both layer series; no woman has a duplicate sample in a compartment. Subject key: subject number from the title `GSSxxxx_<number>_F_<site letter><layer letter>` (A arm, B buttock, C cheek; E epidermis, D dermis; letter-to-site confirmed against the source field). Group: healthy (all). Age: not in GEO. Collection time: NA for all women. Men (timed, forearm, 20 men, four biopsies each, hhmm clock time 0/600/1200/1800 as given by GEO, whether sunrise-adjusted is not stated): epidermis 79 samples, dermis 79 samples (38 men have 4 per layer, 2 have 3); collection time counts 0: 38, 600: 40, 1200: 40, 1800: 40.
Files: `output/wu2020_GSE139300_GSE139301_samples.csv` (sample, subject, layer, site, batch, sex arm, collection time), `wu2020_design.csv`, `wu2020_women_subject_by_compartment.csv`, `wu2020_women_batch_by_layer.csv`, `wu2020_women_batch_by_site.csv`.

## 4. Normalisation state

| dataset | matrix | notes |
|---|---|---|
| GSE45642 | log2 RMA, custom CDF hgu133ahsentrezg, 11,979 Entrez-gene probesets; range 3.5-14.9 (median 6.55) | one global RMA on all chips; regions and sites not corrected; no batch correction |
| GSE71620 | log2 RMA, quantile-normalised (Expression Console 1.2.1.20), 33,297 transcript clusters; range 0.72-13.2 (median 4.58) | one run for both regions, no batch field |
| GSE139301 / GSE139300 | log2 RMA, 49,386 probesets each; medians 2.98 (epidermis) and 4.37 (dermis) | separate RMA per layer series, no batch correction, six batches confounded with compartment |

## 5. Clock-gene probes and detection (percentile of the compartment's median probe level among all probes' medians in that compartment)

Probes per gene: GSE45642 one per gene, **CIART absent** (not on HG-U133A); RORA present. GSE71620 one transcript cluster per gene, all 14 present (PER3 probe 7897378 shared with UTS2, NR1D1 probe 8014956 shared with THRA). GPL13667 (U219, 2010 annotation): several probe sets per gene (ARNTL 3, CLOCK 2, NPAS2 4, PER1 3, PER2 5, PER3 3, CRY1 1, CRY2 2, NR1D1 1, NR1D2 4, DBP 2, RORA 3, TEF 2, CIART 1 via alias C1orf51), some of them _x_/_s_ non-specific sets; the percentile range over probe sets is shown. Percentiles are not comparable across platforms (probe universes differ: 11,979 expressed-gene probesets vs 33,297 vs 49,386).

**GSE45642 (percentile)**

| gene | probes | amygdala | anterior cingulate | cerebellum | dorsolateral prefrontal cortex | hippocampus | nucleus accumbens |
|---|---|---|---|---|---|---|---|
| ARNTL | 1 | 45.0 | 60.0 | 66.0 | 64.0 | 47.0 | 55.0 |
| CLOCK | 1 | 85.0 | 84.0 | 93.0 | 86.0 | 86.0 | 81.0 |
| NPAS2 | 1 | 46.0 | 50.0 | 23.0 | 46.0 | 45.0 | 57.0 |
| PER1 | 1 | 86.0 | 86.0 | 87.0 | 87.0 | 89.0 | 90.0 |
| PER2 | 1 | 40.0 | 37.0 | 35.0 | 37.0 | 35.0 | 41.0 |
| PER3 | 1 | 68.0 | 58.0 | 86.0 | 55.0 | 62.0 | 57.0 |
| CRY1 | 1 | 31.0 | 33.0 | 82.0 | 23.0 | 40.0 | 41.0 |
| CRY2 | 1 | 79.0 | 80.0 | 87.0 | 83.0 | 79.0 | 83.0 |
| NR1D1 | 1 | 89.0 | 91.0 | 94.0 | 93.0 | 87.0 | 90.0 |
| NR1D2 | 1 | 56.0 | 73.0 | 85.0 | 82.0 | 76.0 | 77.0 |
| DBP | 1 | 73.0 | 70.0 | 84.0 | 70.0 | 62.0 | 83.0 |
| RORA | 1 | 24.0 | 31.0 | 48.0 | 14.0 | 13.0 | 14.0 |
| TEF | 1 | 23.0 | 27.0 | 30.0 | 31.0 | 31.0 | 39.0 |
| CIART | 0 | absent | absent | absent | absent | absent | absent |

**GSE71620 (percentile)**

| gene | probes | BA11 | BA47 |
|---|---|---|---|
| ARNTL | 1 | 77 | 78 |
| CLOCK | 1 | 79 | 80 |
| NPAS2 | 1 | 93 | 93 |
| PER1 | 1 | 74 | 73 |
| PER2 | 1 | 83 | 82 |
| PER3 | 1 | 84 | 84 |
| CRY1 | 1 | 61 | 60 |
| CRY2 | 1 | 75 | 76 |
| NR1D1 | 1 | 93 | 93 |
| NR1D2 | 1 | 80 | 81 |
| DBP | 1 | 73 | 72 |
| RORA | 1 | 75 | 76 |
| TEF | 1 | 82 | 82 |
| CIART | 1 | 67 | 65 |

**GSE139301 epidermis (percentile; range over probe sets)**

| gene | probes | epidermis_buttock | epidermis_cheek | epidermis_forearm |
|---|---|---|---|---|
| ARNTL | 3 | 55-83 | 55-87 | 55-85 |
| CLOCK | 2 | 34-78 | 38-78 | 38-78 |
| NPAS2 | 4 | 17-69 | 9-69 | 17-71 |
| PER1 | 3 | 20-81 | 20-67 | 16-74 |
| PER2 | 5 | 4-91 | 4-88 | 3-90 |
| PER3 | 3 | 40-94 | 42-93 | 46-95 |
| CRY1 | 1 | 91 | 91 | 91 |
| CRY2 | 2 | 47-80 | 43-71 | 45-77 |
| NR1D1 | 1 | 68 | 56 | 64 |
| NR1D2 | 4 | 56-92 | 53-91 | 55-92 |
| DBP | 2 | 68-71 | 62-65 | 64-67 |
| RORA | 3 | 97-99 | 95-99 | 97-99 |
| TEF | 2 | 44-90 | 43-87 | 43-89 |
| CIART | 1 | 85 | 73 | 81 |

**GSE139300 dermis (percentile; range over probe sets)**

| gene | probes | dermis_buttock | dermis_cheek | dermis_forearm |
|---|---|---|---|---|
| ARNTL | 3 | 34-69 | 37-75 | 39-71 |
| CLOCK | 2 | 33-73 | 35-73 | 35-73 |
| NPAS2 | 4 | 26-59 | 24-59 | 28-62 |
| PER1 | 3 | 20-86 | 13-75 | 22-79 |
| PER2 | 5 | 13-92 | 9-90 | 13-90 |
| PER3 | 3 | 40-96 | 36-95 | 43-95 |
| CRY1 | 1 | 90 | 89 | 87 |
| CRY2 | 2 | 35-75 | 34-72 | 39-77 |
| NR1D1 | 1 | 65 | 54 | 66 |
| NR1D2 | 4 | 58-93 | 54-92 | 56-93 |
| DBP | 2 | 59-62 | 56-58 | 58-61 |
| RORA | 3 | 64-96 | 61-95 | 62-96 |
| TEF | 2 | 30-88 | 30-86 | 35-88 |
| CIART | 1 | 70 | 64 | 67 |

Full per-probe table: `output/brain_skin_clock_probe_detection.csv` (median log2 and percentile per probe and compartment).

## 6. Within-person comparability summary

- GSE45642: subject IDs consistent across regions (good); 28 subjects have all six regions, 64 have four or more regions (24+12+28); but regions and replicate chips are spread over three assay sites and two chip types, no TOD in GEO, no MDD in GEO.
- GSE71620: pairing BA11/BA47 recoverable for all 210 subjects, two compartments only, TOD for 146 subjects, adjacent PFC areas (weak between-organ contrast).
- Skin: best design (154 women, 137 complete with six compartments, three sites x two layers, untimed, healthy), but compartment-batch confounding, no age in GEO, no dermis/epidermis common normalisation; timed arm limited to 20 men on forearm.
