# Preliminary results (append only)

Results of the amendments registered in `registered_experiments.md` (A109 onwards), each with its verdict and a
`<!-- src: ... -->` pointer to the output table it was read from. Results A1-A108c are in the archive
(`$SCCC_ARCHIVE/docs/RESULTS_preliminary.md`, frozen); corrections to them go in the next section, dated, with the
archive section they correct.

## Errata to the archive

### 2026-10-23, WP0 item 7: A22 V17.2 coherence values (corrects archive RESULTS "V17.2, A21 colon", line 1262)

The legacy V17.2 coherences were computed on a pseudobulk table whose empty `stratum` column was read as NaN, so pair
coherence was not centred within stratum (legacy e_control 0.16-0.32; delta -0.11, SE 0.07; case excess 0.044). The
port (`results/extended_data/scripts/simulations/compute_a22_colon.py`, colon as one centred stratum) gives e_control
0.45-0.55 over the four anchors (lymphatic EC 0.448, endothelial 0.517, smooth muscle 0.544, fibroblast 0.546) and
delta -0.17 to -0.40; lymphatic EC: delta -0.399 (SE 0.164), e_case 0.048. The verdict is unchanged: no anchor is
dropped, lymphatic EC is not dropped, V17.2 fails as in the archive. Other legacy coherence runs on tables with an empty
stratum may carry the same artefact (to be checked where cited).
<!-- src: paper/results/extended_data/output/simulations/a22_colon_confirm.csv (archive-port snapshot, sha256 e00e28cc...) -->
