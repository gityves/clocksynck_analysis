# feillet2014phase
**Ref:** Feillet, Krusche, Tamanini, Janssens, Downey, Martin, Teboul, Saito, Levi, Bretschneider, van der Horst, Delaunay, Rand, PNAS 2014 (111:9828-9833). DOI 10.1073/pnas.1320474111. Type: biology + model. Lab: Delaunay (Nice), Rand (Warwick), van der Horst (Rotterdam). NOT a Naef-lab paper: Naef is not an author (text has no occurrence).
**Source:** full text

## Question and design
NIH3T3 fibroblasts with a Rev-erbα::Venus clock reporter and a FUCCI cell-cycle reporter, single live cells, tracked 72 h with lineage trees; FBS 10, 15, 20% to change cell-cycle length; unsynchronised cells versus a dexamethasone pulse. **Sampling:** live imaging, continuous, 72 h. Cell numbers are in the SI (not in the text read).

## Method (if any)
Joint phase of clock and cell cycle on a torus; stochastic model where "the cell cycle progresses at a speed that depends upon the phase of the clock in that it is slowed down" in a region of the torus (clock acts on cell cycle); winding ratios p:q.

## Key results
- Unsynchronised cells lock 1:1 with a common period that scales with serum: 10% FBS clock 21.9 +/- 1.1 h and cell cycle 21.3 +/- 1.3 h; 15% FBS 19.4 +/- 0.5 and 18.6 +/- 0.6 h. "We show that in unsynchronized cells the cell cycle and circadian clock robustly phase lock each other in a 1:1 fashion"
- Dexamethasone reveals coexisting states. "there are distinct coexisting states with a significantly higher period clock"; in 20% FBS: "had median periods for the clock and cell cycle of 27 h and 17 h, respectively (Fig. 4C), roughly a 3:2 period ratio", other group "The cells in the second group locked 1:1".
- 10% FBS after dexamethasone: clock 24.2 +/- 0.5 h, cell cycle 20.1 +/- 0.94 h, about 5:4 (modelled), no multi-peak division pattern; the three division gates of Nagoshi 2004 are explained by 3:2 locking, not gating. "The gating hypothesis is also not consistent with the dynamics seen in our unsynchronized cells."
- Multistability: dexamethasone "knocks the state out of the low-period 1:1 locked state" into other attractors.
- Reverse direction suggested: "suggest that the coupling also operates in the opposite direction" (both periods fall together as FBS rises).

## Not addressed / limits
- Modes observed: 1:1, about 3:2, about 5:4. A 2:1 state is not reported. Higher p:q "decreases very rapidly" in robustness.
- Cause of reverse coupling "is not currently supported by known mechanistic links"; no genetic perturbation of the clock.
- One cell line, in vitro; dexamethasone is a strong non-physiological perturbation; no human cells, no tissue.

## Link to us
- **Intro claim(s):** refs_cancer_phase_inference (1:1, 3:2, 2:1 attractors, tuned by period mismatch and coupling; Naef lab) is partly wrong: modes are 1:1, 3:2 and about 5:4 (2:1 is only predicted by Droin 2019 and not seen); tuning is by serum (cell-cycle length) and dexamethasone, not directly by coupling strength; and the group is not Naef's. Supports that clock and cell cycle lock in single cells.
- **Relation to clocksync:** background; explains that proliferative states can change the clock period (21.9 versus 24.2 h) and hence relative phase.

## Ideas
- Cite Feillet and Bieler together as the two independent demonstrations of 1:1 locking; state the disagreement on direction (clock-to-cell cycle in the Feillet model; cell cycle to clock in Bieler and Droin).
- Missing reference: Nagoshi 2004 (already in notes as nagoshi2004circadian).
