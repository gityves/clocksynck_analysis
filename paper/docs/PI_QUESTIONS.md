# Questions for the PI (open first, answered below)

Each entry: date, question, what is on hold until it is answered, the evidence (with `src` pointers in
`results_preliminary.md`). Work that does not depend on an answer continues. When answered, move the entry to
"Answered" with the answer and its date.

## Open

### Q1 (2026-10-23): Figure 3 framing after A109
H109.1 passes (organs share a clock-shaped state in every death group); H109.2 fails in ventilated donors (phase
relations collapse to 0 / 12 h; about half of the shared covariance is a displacement) but passes in Hardy 4
(secondary). Registered reading: "organs share a donor state", WP2 out of the main text. Proposed alternative: coupled
running clocks without cues, shown in Hardy 4 / Hardy 3, with ventilated donors as a shared displacement overriding
the clock relations. PI comment (2026-10-23): check first that Hardy 4 clocks are decorrelated from wall time, and
whether the 0 / 12 h push is seen in COVID ventilated patients -> registered as A109c (dab55ab), being scored.
- **A109c H109c.1-2 (2026-10-23):** no detectable wall-time alignment in Hardy 4 (R 0.097, p 0.56; MDE 0.145-0.208);
  offsets preserved after removing wall time (rho 0.750, p 0.0002), also without ventilator cases and in inpatients
  only. COVID (D109c.3-4, descriptive): Delorey organs uninformative (too few shared donors); Melms lung shows no
  push in COVID vs controls. The circularity ratio is confounded by host-phase concentration (sampling window), so the
  offset test (H109.2) is the decisive readout; the ventilated collapse to 0 / 12 h rests on it (p_closer 0.001).
- **Question:** after A109c, which framing (registered default, or the Hardy 4 framing as a registered amendment
  labelled post hoc relative to A109)?
- **On hold:** WP2 (A110); edits of the Figure 3 table and narrative.
- Evidence: results_preliminary.md, sections A109, A109b (and A109c when scored).

### Q2 (2026-10-23): figure-table numbers used against their archive verdict (WP0 item 1)
- A104 rho +0.73 quoted as support; H104.1 predicted rho < 0 (archive: not supported; H104.2 not replicated).
- A54 psoriasis "replicated in discovery and held-out": the delay replicates (H54.2) but H54.1 and H54.3 fail.
- Mouse embryo "from E15": H108.1 passes in the late window but the window gate does not (p_acrophase 0.08).
- Quoted values that differ from their table: benchmark offset r 0.934 / 0.900 (quoted 0.91 / 0.88); baboon offsets
  0.455 h (0.5 h); healthy atlases 22 targets / 17 studies ("about 30").
- **Question:** remove, restate or keep each claim (with the verdict stated)?
- **On hold:** edits of the Figure 1, 2 and 5 tables.

### Q3 (2026-10-23): panels not on the anchored A17+ estimator (WP0 item 4)
13 of 30 figure-table rows flagged (`paper/exploratory/WP0_integrity/output/estimator_consistency.csv`). By design
(the text must name the readout): r* (A103, A42), A104, A85 index, A92/A93 coherence C, A33 tissue host. Needing a
decision: Fig 1 A21 colon (an anchor decouples, V17.2 fails); Fig 2 GTEx organ departures (pre-A17 timed) and
neuron/retina panels (pre-A17 free fit); Fig 4 A20 "damping" (withdrawn in A22) and COVID A105 (2 anchors, no A22
confirmation); Fig 5 fibrotic-lung T cells (A5 free-fit 0.45).
- **Question:** drop, recompute on the anchored estimator, or keep with the readout named, per panel?
- **On hold:** recomputation of those panels.

### Q4 (2026-10-23): clocksync change
`clocksync.genes.gene_columns` returns -1 silently for a missing symbol (indexing the last column; the A85 bug).
Proposed: `missing="ignore"|"warn"|"raise"` and an all-absent guard
(`paper/exploratory/WP0_integrity/output/gene_audit/proposed_fix.diff`, with a test). clocksync holds your uncommitted
work, so it was not applied.
- **Question:** apply it (own commit, only these files)?
- **On hold:** nothing (no current result affected).

## Answered
