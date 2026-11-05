# anafi2017cyclops
**Ref:** Anafi, Francey, Hogenesch, Kim, PNAS 114:5312 (2017). DOI 10.1073/pnas.1619320114. Type: method + data.
Lab: Ron Anafi (Penn) / Hogenesch.
**Source:** full text (main; SI not in export)

## Question and design
Order untimed human samples on the circadian cycle (CYCLOPS: circular-node autoencoder on eigengenes of seed genes =
human homologs of genes cycling in >75 % of mouse tissues). Validation: mouse time courses (liver etc.), 146 human
prefrontal cortex autopsy samples with time of death. Applications: lung biopsies (Groningen, Laval), liver: 249
non-cancerous biopsies (243 tumour-adjacent "normal margin"), HCC biopsies; STZ chronotherapy in mice.
**Sampling:** biopsies, surgical working hours; autopsy brain with TOD; mouse full 24 h.

## Method
One phase per sample from an ellipse in eigengene space (85 % variance); cohort-level; relative phase anchored by
the PAR bZip acrophase set to π. Quality: Metsmooth, Staterr. Assumes a common phase relationship of seed genes
across subjects. HCC was ordered by projecting onto the eigenvectors of normal samples.

## Key results
- Brain TOD: "provides an excellent estimate of time of death (TOD) (ρc = 0.68, ηc = 0.55, median absolute error = 1.69 h)".
- Lung, two sites: acrophases consistent, "median absolute discrepancy = 0.32 radians ∼1.2 h".
- Inferred sample phases not uniform: "temporal reconstruction with CYCLOPS did not uniformly distribute samples across the circadian cycle".
- Hospital hours stated explicitly: "Biopsies are obtained during surgical working hours (∼6:00 AM to 6:00 PM)."; night phases attributed to shift workers ("∼15–20% of the population are shift workers").
- Desynchrony as hypothesis: "resulting in intertissue (36) or intratissue (30) desynchrony" (36 = Archer 2014).
- HCC: no quality ordering from mouse seeds; after projection onto NC eigenvectors most core clock genes cycle except PER1, CRY1; "we estimate that ∼15% of the transcripts that cycled in NC samples lost rhythmic expression in HCC".
- **Full-cycle requirement (key quote):** "It requires data from the entire periodic cycle to form an ellipse. Biopsies are almost exclusively obtained during the day." and "Our experience suggests that >250 samples are required to order biopsy samples".
- Phase is relative: "Additional information is needed to assign a circadian time to any particular CYCLOPS phase."
- Statistics on inferred order are optimistic: "standard statistical significance tests tend to be too" (liberal).

## Not addressed / limits
- Shift-worker explanation of night phases is an assumption, not tested; an alternative is that hospital-hours data are forced onto a full ellipse.
- One phase per sample; tissues ordered separately, so no within-person comparison of tissues.
- Tumour ordering needed the normal-tissue projection: tumour phases are not estimated independently.

## Link to us
- **Intro claims:** (i) orderers need the whole cycle: this is THE citation, verbatim; (ii) biopsies come from working hours (6 am-6 pm); (iii) >250 samples; (iv) tumours keep a clock with altered output (HCC). Our refs table says no explicit full-cycle requirement was quoted: WRONG, update it.
- **Relation:** competitor (orderer); already in A145.

## Ideas
- Quantify in our data the paradox they explain by shift work: run CYCLOPS/CHIRAL on a hospital-hours subset of GTEx (A123 window) and show the inferred phases still span the circle (spurious full-cycle ordering).
- Use their HCC data (Lamb 2011, GSE25097? check) as a tumour + adjacent normal pair cohort for clocksync.
