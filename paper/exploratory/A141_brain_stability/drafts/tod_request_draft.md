# DRAFT (not sent; PI decides): request for time-of-death fields, brain cohorts (A141)

Why: A141 found 4 brain disease contrasts whose dependence on the host definition exceeds sampling noise. They are
PD_mssm, VaD_mssm (PsychAD MSSM), ASCVD_hbcc (PsychAD HBCC) and the PD atlas. With the time of death, each compartment
can be read against clock time directly. That shows which compartments follow time, and so which host definition is
right. Needed: time of death (hour; date or season optional) per donor. Aggregates only would be released.

Targets:
- PsychAD (MSSM, HBCC): Synapse syn52160016 (controlled; data-use request through the PsychAD / NIMH data access process).
- PD atlas (Prashant et al. 2024): AMP PD Knowledge Platform (data use agreement; clinical tables).
- Optional for coverage: ROSMAP clinical (AD Knowledge Portal), for RADC.

Draft message (to the data access committee / corresponding authors):
> We are using the PsychAD single-nucleus data to study the internal coordination of circadian clocks across brain cell
> types without sample times. To validate the method's readouts against clock time, we request the time of death (hour)
> for the donors already in our approved dataset. We will report aggregate statistics only and will not redistribute
> donor-level values.

## Update 2026-11-16 (cell-type case; still a DRAFT, the PI decides)
**What exists (checked at the source, 2026-11-16).**
- ROSMAP single-nucleus prefrontal cortex, 427 participants (Mathys et al., Cell 2023; >2.3 M nuclei): raw data and
  metadata on Synapse (syn52293417; participant metadata syn3157322 or through the RADC Research Resource Sharing Hub),
  Data Use Certificate required (https://compbio.mit.edu/ad_aging_brain/).
- Time of death in ROSMAP: recorded for ROS and MAP donors (all are brain donors) and used by Lim et al., Nat Commun
  2017 (DLPFC bulk RNA-seq; time of death relative to sunrise and local clock time). It is a RADC variable, requested
  through the RADC hub, not part of the open Synapse metadata (to be confirmed with RADC).
- Our `psychad_radc` dataset is the RADC (ROSMAP) part of PsychAD; whether its donors overlap the Mathys 2023
  participants is to be checked from the IDs once access exists.
**What the times would test (added scientific case).** The anchor-free readout shows that the cell classes of the
human cortex share a clock-gene state (A274, H274.3; pending the calibration of Q36), and the order of cortical cell
classes disagrees in sign between PsychAD and SEA-AD (neuron-vascular offset +1.3 h against -1.5 h). With the hour of
death per donor: (1) the A273 calibration in humans (untimed pairwise offsets between cell classes against the timed
offsets of the same donors; so far mouse colon only); (2) the tie-break of the cortex offset sign; (3) whether the
shared cortical state reads the time of death (as the clock host does in GTEx).
**Draft paragraph to add to the message:**
> We additionally request the hour of death (local clock time) for the ROSMAP participants of the single-nucleus
> prefrontal cortex dataset (Mathys et al. 2023), through the RADC Research Resource Sharing Hub. The values would be
> used only to validate untimed estimates of clock phase between cell types against time of death; only aggregate
> statistics would be reported.
