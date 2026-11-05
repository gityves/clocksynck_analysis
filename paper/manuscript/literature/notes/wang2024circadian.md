# wang2024circadian
**Ref:** Wang et al., Cell 187:2690-2702 (2024). DOI 10.1016/j.cell.2024.04.015. Type: biology + translational (mouse, human melanoma). Lab: Scheiermann (Geneva), Pittet.
**Source:** full text

## Question and design
Does the tumour microenvironment oscillate over the day and does it dictate immunotherapy efficacy? Mouse B16-F10-OVA (s.c.) and spontaneous Tyr::CreERT2;Pten;BrafCA melanoma, DoHH2 lymphoma
xenografts with human CD19 CAR T cells, MC38 colon; tumours harvested at 4 zeitgeber times (ZT1, ZT7, ZT13, ZT19, n = 11-15 per time), scRNA-seq of CD45+ cells at the same 4 times, DD and 12 h-inverted light, acute jet lag, endothelial Bmal1 deletion (Cdh5-CreERT2), T-cell Bmal1 deletion. Human: melanoma resection
specimens at different clock times (CD4 n = 18, CD8 n = 19 patients, IHC), TCGA SKCM (103 primary samples), public skin microarray.
**Sampling:** timed. Mouse 4 times per 24 h (6 h spacing); human surgical hours only, with recorded time of removal for IHC.

## Method (if any)
ZeitZeiger (supervised time prediction), trained on mouse tumour and skin scRNA/bulk (57 conserved genes); applied to untimed human samples. Assumes conserved acrophases between mouse and human, and a fixed 24 h clock in the sample type.

## Key results
- Mouse TILs oscillate: total leukocytes peak at ZT13; "Constant darkness did not alter the observed time-of-day differences", light inversion inverts them, jet lag abolishes them (Fig 1).
- Clock genes rhythmic in tumour leukocytes: "Circadian clock genes were expressed rhythmically in all leukocyte subsets with Per1 exhibiting the most wide-spread oscillation, peaking at ZT13" (Fig 3A); each subset has distinct oscillating genes.
- Host endothelial clock gates infiltration: "time-of-day differences in ICAM-1 expression in tumor endothelial cells were completely abrogated" in Bmal1DEC mice, and infiltration stops being time-dependent (Fig 2).
- CAR T cells: "infusion of CAR T cells at ZT13 showing a reduction in tumor volume, while infusion at ZT1 had no effect on tumor control" (Fig 2I); OT-I therapy better at ZT13, lost in Bmal1DEC (Fig 2H).
- PD-1: "Pdcd1 (encoding for PD-1) was rhythmic in this CD8 T cell cluster, peaking at ZT1 and exhibiting a trough at ZT13"; anti-PD-1 works at ZT13 not ZT1, lost in T-cell Bmal1 deletion (Fig 5).
- Human: CD8/CD4 counts peak "in the early afternoon" (IHC, cosinor, Fig 6A, B); ZeitZeiger predicted time "with an average error of 1.5 h" in held-out mouse data and "we could adequately predict the time of day in these samples with an error of <4 h" in healthy human skin; "we next assigned time-of-day information to 103 primary melanoma samples from the TCGA dataset", exhausted/non-exhausted CD8 ratio higher anti-tumour in the morning than the afternoon: "with a higher anti-tumor ratio in the morning compared with the afternoon" (24 vs 42 patients, Fig 6F).

## Not addressed / limits
- Tumour-cell clocks are not the subject: oscillation is in host-derived leukocytes and endothelium; tumour-cell rhythms not characterised. No tumour vs adjacent-normal comparison; no specificity control for the time-stamped TCGA result beyond skin validation.
- Human time stamps for TCGA are predicted by a mouse-trained model validated only in healthy skin; error in tumours unknown, and a tumour with a damped clock would be mis-timed. "We focused on melanoma in this study".
- Mouse optimum ZT13 (active phase start); human optimum assumed morning from phase inversion, not measured. Sample times in humans are surgery hours only.

## Link to us
- **Intro claim(s):** refs line 82 (CD8 T-cell infiltration and function oscillate via leukocyte clocks and endothelial-clock-dependent trafficking; ICB and CAR-T efficacy improve with time-of-day adjustment): supported (mouse; human is IHC counts and predicted times).
- **Relation to clocksync:** background and comparator (time stamping from untimed tumour biopsies with ZeitZeiger); supports a host-referenced view (host endothelium/immune clocks drive the tumour-resident rhythm, stromal-vascular host anchor).

## Ideas
- Their 57-gene ZeitZeiger as a comparator against clocksync on the same human tumour sets: A145 compares two-step orderers (CHIRAL/CYCLOPS/CoFE) with clocksync; ZeitZeiger on TCGA/skin not found in the docs (not done).
- Immune and endothelial compartments of tumours vs stromal host: tumour T/NK and myeloid readouts done in A138 (P4-P6, immune) and A140 (CYCLOPS 2.0, immune compartments); endothelial-specific host readout not found.
