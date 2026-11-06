# A212: coherence mechanisms (plasmode)

Registered in `paper/docs/registered_experiments.md` (A212, commit e419f38). Run: `python scripts/compute_a212.py --stage grid`, then
`--stage clouds`, then `--stage tests` (about 1 min in total with 4 jobs).

Choices left open by the registration:
- The clock part uses the gene-centred waveform design (clocksync `plane`, the projector the A176 statistics use), so C = Y P;
  the registration wrote beta = Y X (X'X)^-1 without saying whether X is gene-centred.
- Swamping loadings: SD = s x median per-gene SD of E (one scalar per base), s on linspace(0, 2, 10).
- Matching: scan the grid from weakest to strongest; the first value whose mean dR is at or below the target, interpolated
  with the previous one. Mechanisms that never reach the target run at the grid value with the lowest mean dR.
- W draws one shift vector per replicate (shared by all case tumours of all strata); H draws one per tumour.

Property found while running (not a test): dR is not monotone in the strength of D, W and H. With strong damping (rho <= 0.4) or
large shifts (>= 4-5 h) coherence rises again, because the leave-one-out phases then come from the off-plane rest, which is
orthogonal to the plane by construction. The matched strengths are on the monotone branch (rho 0.83, shifts 1.8-2.5 h).

Outputs: `grid.csv`, `matched.csv`, `clouds.csv`, `tests.csv`, `auc_all.csv`, `observed.csv`.
