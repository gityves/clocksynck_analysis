# shaashua2020stromal
**Ref:** Shaashua et al., Front Cell Dev Biol 8:587697 (2020). DOI 10.3389/fcell.2020.587697. Type: biology (mouse). Lab: Scherz-Shouval (Weizmann).
**Source:** full text

## Question and design
Is the clock in the non-malignant host (tumour microenvironment, TME) needed for tumour growth? Wild-type, Per1-/-, Per2-/- and
Per1-/- Per2-/- C57BL/6 mice receive wild-type cancer cells: MC38 colon via portal vein (liver metastasis, 21 days; n = 9-16 per
genotype), E0771 breast in mammary fat pad (day 32; n = 11-13), MC38 in colon wall (day 14; n = 10-14). Bulk liver RNA-seq 1 week after
MC38 or PBS (n = 4 and 2 per genotype). In vitro Per2::luc rhythms of NIH3T3 vs MC38, LLC, E0771.
**Sampling:** terminal single time points; no circadian time reported; in vitro luminescence time courses after dexamethasone.

## Method (if any)
Genetic loss of Per genes in the whole host (not cell-type specific); tumour burden and histology (Sirius red, alpha-SMA CAFs);
CIBERSORT with single-cell reference for liver composition.

## Key results
- In vitro: "all three cancer cell lines were less rhythmic than NIH3T3"; periods LLC 28.3 h, E0771 27.2 h, MC38 26.3 h vs fibroblasts, and Per1/Per2 knockdown impairs fibroblast rhythm.
- Per1-/- Per2-/- hosts have fewer metastases and lower liver/body weight; Per2-/- alone "exhibited significantly lower metastatic burden than WT mice" while "Loss of Per1, on the other hand, did not affect metastasis" (Fig. 3).
- Primary tumours: Per2-/- hosts "exhibited significantly smaller E0771 tumors compared to WT mice", "a similar trend was seen for the MC38 colon model" (Fig. 5).
- T cells are not the mediator: Per1-/- Per2-/- remain resistant after CD4/CD8 depletion.
- Liver RNA-seq: 1222 differentially expressed genes between genotypes and conditions; Per2-/- livers upregulate autophagy, circadian regulation and anti-proliferation genes after cancer-cell injection; hepatocytes ">95%" of cells.
- Interpretation: "thus implying that the effect of stromal Per2 knock-out on metastatic colonization may be clock-independent".

## Not addressed / limits
- Constitutive whole-body knockouts, not stroma-specific; "stromal" is the entire host, including liver parenchyma and the central clock.
- No clock readout (phase, amplitude, coupling) in host or tumour in vivo; whether the knockout disrupts tumour-host synchrony is untested.
- Per1 vs Per2 difference argues the Per2 effect may not be a rhythm effect; no tumour-vs-adjacent-normal comparison; no time of sampling.
- Opposite to the usual tumour-suppressor reading of Per2 (cell-autonomous), so direction depends on compartment.

## Link to us
- **Intro claim(s):** host clock acts on tumours, host and tumour clocks can have opposite roles. Our refs row ("stromal (not tumour) Per2 required for tumour initiation and metastatic colonisation") is SUPPORTED but should say "host Per2 (global knockout)" and note the authors' clock-independence caveat.
- **Relation to clocksync:** background; the authors pose "synchronization, or lack of synchronization, between the TME and the cancer cells" as a hypothesis, which is what a host-referenced readout tests.

## Ideas
- Test in human data whether stromal-vascular coupling to the host is retained in tumours: done for fibroblasts and endothelium (A151 post hoc, tumour stroma vs same-patient normal stroma clock-specific; A138.4 TCGA pairing). Not done: a stroma-specific Per2 conditional knockout design.
- Missing reference: Fuhr 2018/2019 (CAFs improve cancer-cell rhythmicity in co-culture), cited here.
