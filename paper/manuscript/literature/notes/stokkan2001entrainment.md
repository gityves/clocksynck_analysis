# stokkan2001entrainment
**Ref:** Stokkan, Yamazaki, Tei, Sakaki, Menaker, Science 2001. DOI 10.1126/science.1060081. Type: biology. Lab: Menaker.
**Source:** full text (the file starts with the end of an unrelated paper; the Stokkan paper follows)

## Question and design
Does a feeding cycle entrain peripheral clocks independently of the SCN and the light cycle? **Rats** (Per1-luc
transgenic, young animals) on 12:12 LD, food available 4 h or 8 h per day in the light phase; after 2, 7 or 19 days of
restricted feeding (RF) the liver, lung and SCN were explanted and luminescence was recorded in vitro. Per time point and
tissue only a few animals (e.g. four liver explants on day 2; sample sizes in the Fig. 3 legend). Controls: ad libitum;
corticosterone injections (5 mg i.p., 7 days, ZT5). **Sampling:** timed; one explant per animal, phase read from the first
full peak 12-36 h after explant (animals killed at a fixed clock time).

## Method (if any)
Tissue phase = peak of the in vitro luciferase rhythm 12-36 h after explantation, taken to reflect the phase in vivo at removal.

## Key results
- SCN is not shifted by RF: "the phase of the SCN rhythm was unaffected and remained phase-locked to the light cycle, even after 19 days of RF." (Results)
- Liver 4 h RF: average shift of 10 h on day 2 (n = 4, more variable), 12 h after 7 or 19 days (Fig. 3A).
- Lung 4 h RF: on day 2 two of four explants were arrhythmic; 6 h shift after 7 days, no further shift after 19 days.
- 8 h RF (less stress): liver shifted 12 h after 7 days, lung unshifted; fasting 1-2 days after RF keeps the shifted liver phase (entrainment, not masking).
- Corticosterone injections (15-fold rise, 581 +/- 174 ng/ml) do not shift liver or lung phase (Fig. 3C): not the mediator.
- Abstract: "restricted feeding rapidly entrained the liver, shifting its rhythm by 10 hours within 2 days." They propose that peripheral oscillators may be coupled to the SCN mainly through rhythmic behaviour such as feeding.

## Not addressed / limits
- Rats, not mice (our refs entry says mice). Per1 promoter reporter in explants, not endogenous gene expression in vivo; a handful of animals per cell.
- Only liver, lung and SCN; no inter-organ phase map within one animal under RF.
- Signal responsible (insulin, temperature, food intake) not identified; the lung shift is ascribed possibly to activity.
- Corticosterone tested at only one time (ZT5), where dexamethasone also has no effect in balsalobre2000resetting.

## Link to us
- **Intro claim(s):** feeding is a zeitgeber that can decouple liver from SCN (Internal desynchrony). refs entry says
  mice with ~12 h shift (Per1-luc reporter): **CLAIM WRONG on species** (rat); the 12 h shift after >= 7 days is right,
  10 h is the day-2 value.
- **Relation to clocksync:** background. Meal timing and fasting in hospitalised patients are a plausible source of
  organ-specific phase offsets.

## Ideas
- Their shift in liver but not SCN is the case a host-referenced offset (arg z of liver versus brain-like host) should detect: use as a simulation benchmark of a 10-12 h offset.
- Mouse counterpart with several tissues in one experiment: damiola2000restricted. GTEx phenotype proxies for feeding
  (ventilation, fasting before death) could be checked for offsets (aggregates only).
