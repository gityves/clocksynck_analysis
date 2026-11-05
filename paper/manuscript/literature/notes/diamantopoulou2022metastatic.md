# diamantopoulou2022metastatic
**Ref:** Diamantopoulou et al., Nature 607:156-162 (2022). DOI 10.1038/s41586-022-04875-y. Type: clinical + biology (patients, mouse). Lab: Aceto (ETH Zurich) with Basel clinicians.
**Source:** full text

## Question and design
When do circulating tumour cells (CTC) enter blood? 30 hospitalised women with breast cancer (21 early stage, 9 stage IV), "no treatment or were temporarily off-treatment", blood at 10:00 and 04:00 of the same day; four mouse models (NSG-CDX-BR16, NSG-LM2, NSG-4T1, BALB/c-4T1) with terminal blood every 4 h over 24 h; scRNA-seq of CTCs at ZT4 vs ZT16 and of patient CTCs (10:00 vs 04:00); tail-vein injection of CTCs sorted by phase; jet lag, melatonin, luzindole, dexamethasone, testosterone, insulin; Bmal1-/- E0771.lmb mice. No adjacent normal.
**Sampling:** timed. Patients: two clock times (10:00, 04:00; paired within patient). Mice: 24 h at 4 h intervals (ZT0, 4, 12, 16 plus others).

## Method (if any)
None (CTC capture by microfluidics, scRNA-seq, GSEA of rest vs active signatures).

## Key results
- Patients: "we found most CTCs (78.3%) in samples obtained at night time during the rest phase, including single CTCs, CTC clusters and CTC–white blood cell (WBC) clusters" (Fig 1a).
- Mice: 87-99% of CTC events in rest phase across the four models; ZT4 vs ZT16 "with a 6- to 88-fold increase for single CTCs" (Fig 1c); same in the tumour draining vessel.
- Jet lag: "jet-lag induction leads to a 38- to 282-fold decrease in single CTCs" with unchanged primary tumour size; Bmal1-/- host "mice generally fail to generate CTCs despite identical tumour size and timing of sample collection" (Ext Data Fig 4).
- Metastatic potential: "CTCs obtained during ZT4 exhibit an extraordinary metastasis-forming capacity compared to CTCs that are obtained during ZT16" (Fig 2e); clusters stronger than single CTCs.
- Transcriptome: mitosis programmes up at rest, ribosome biogenesis/translation up in active phase, also in patients; Ki67 higher at rest in primary tumours (Fig 3).
- Tumour clocks: "we could not detect differential expression in CTCs or primary tumour cells" of clock genes between ZT4 and ZT16, whereas rhythmicity was confirmed in non-neoplastic tissues (Ext Data Fig 8). Systemic hormones (melatonin, testosterone, glucocorticoid, insulin) drive timing.

## Not addressed / limits
- Patient samples have two time points; whether clock genes are rhythmic in patient tumours is not measured; the clock-gene null in tumour cells rests on ZT4 vs ZT16 (a 12 h contrast cannot show damped rhythm with a different phase).
- Metastasis assay is mouse only (CTC transfer, bioluminescence); patient claim is CTC counts and gene expression, not metastatic ability.
- Sleep in patients not monitored; CTC capture efficiency and clearance by phase tested only in mice. n per mouse time point small (n = 3-6).

## Link to us
- **Intro claim(s):** refs line 84 (CTC intravasation peaks during sleep; rest-phase CTCs far more metastatic): supported for counts in patients and mice; the higher metastatic ability is demonstrated in mice.
- **Relation to clocksync:** background and motivation. Host rhythms (hormones) rather than tumour-intrinsic clocks drive tumour timing, in line with host-referenced readout of compartments; hospital-hours sampling at a single clock time would miss it.

## Ideas
- Tumour-cell clock absence plus host-driven proliferation timing parallels A138/A149 (malignant clock tracks the stromal host, not clock-specific beyond tumour-normal bulk match; A149).
- Patient CTC scRNA at 10:00 vs 04:00 is a paired time contrast; request its CTC scRNA data (GEO in paper) to test whether a host-referenced readout shows phase differences.
