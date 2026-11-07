# IBD paired inflamed / non-inflamed scout, update 2026-11-07 (metadata only)

Base: /home/maxine/Documents/paychere/scCircadianCancer/legacy/desync/validation/inflammation_scout_2026-10-16.md (not redone).
Legend: [V] verified this session (GEO text API / paper), [P] paper, [U] unverified. Nothing downloaded except GEO sample metadata text.

## Already seen in this repo
- TAURUS (GSE282122 / Zenodo 13768607; on disk in ibd/TAURUS): analysed already (PI_QUESTIONS Q19, A45): 29 patients paired, 23-24 with both families; gut stroma coupling 0.80-0.91 in both states, phases unchanged. A within-patient epithelium-vs-stroma offset contrast was NOT registered, so a TAURUS readout of the new hypothesis would be non-blind (secondary only).
- Bulk IBD pairs (A55) on disk: ibd_bulk GSE107499, GSE179285, GSE193677, GSE208303, GSE243625, GSE66407, GSE95095.
- Smillie 2019 (SCP259): never acquired (needs the PI's portal login; Q19 option i; PI 2026-10-31: revisit only if no bulk design reaches MDE <= 0.15 on d*).
- No hit in results_preliminary / registered_experiments for GSE190564, GSE250498, GSE266546, GSE214695 beyond the 2026-10-16 scout. Only "ibd_crohn" bulk contrasts.

## New findings this session
1. GSE266616 identified: paper = Zhou, Xie, Basu et al., bioRxiv 10.1101/2025.03.26.645552 / PMID 39837850 (CD on biologics, "key mediators of inflammation"). [P] 29 CD + 17 non-IBD (GEO says 32 CD + 17 controls, 85 samples listed, "79 biopsies"); 48 TI + 37 AC biopsies; epithelial and lamina propria fractions separated enzymatically then POOLED for scRNA-seq => whole-mucosa-like, epithelium AND stroma in one library; 10x 3' v3.1. Data: GEO raw tars open (3.7 GB); no annotation deposited (FASTQ/matrices "on request" in paper). [V] from GEO metadata (patient, tissue, sample region): per-sample regions: inf 17, adjacent 19, non 16, control 33. CD patients with inflamed AND (adjacent or non) at the SAME site: **14** [V]; with inflamed + strictly "non": 1 [V]. So the contrast is inflamed vs adjacent-to-inflamed (adjacent is the usual "uninvolved" comparator; whether it counts as non-inflamed is a design decision). Treatment (anti-TNF etc.) annotated. Collection time: not recorded [U]. Median UMI / cells per sample: not checked [U]; needs own annotation of fibroblast/endothelial. Only candidate besides TAURUS with open counts, whole mucosa, paired same-site.
2. GSE206299 etc. (irColitis, Thomas 2024 Nat Med, PMID 38724705): GEO tissue_cells has only 1 uninflamed sample patient (SIC_43 left colon uninflamed) [V], CD45+ enriched/sorted libraries alongside unsorted, deposited h5ad are lineage subsets (B, CD4, CD8, myeloid) only, raw in dbGaP. NOT usable.
3. GSE253720 (Kattah/Oh, steroid-refractory CPI colitis; PMID 38642938): HC, UC, CPI colitis, pooled/demuxed, one biopsy per patient, h5ad open but no paired uninflamed [U]; raw missing for privacy. NOT usable for pairs.
4. GSE288517: NEW bulk. Laser-capture-microdissected epithelium and lamina propria from paired inflamed and uninflamed biopsies, biologic-naive UC, pre anti-TNF. Epithelium: 10 uninflamed, 13 inflamed; LP: 10 uninflamed, 12 inflamed (about 10 patients with both states) [V]. But raw/human counts not submitted to GEO [V statement]; processed matrix availability [U]. A natural fit for the epithelium-vs-stroma design (two compartments per sample, same patient and biopsy) but n about 10 and counts possibly unavailable. Secondary.
5. Others met in search, not usable or already known: IBDverse (unpaired), Dryad 10.5061/dryad.rfj6q57dw (bulk RNA-seq + proteomics TI inflamed vs non-inflamed paired, bulk, whole tissue only; 1 compartment => out per estimable-bulk-designs rule), Frontiers Immunol 2026 (Eshghi = SCP3755, already listed), PSC-IBD vs UC medRxiv 10.1101/2024.07.23.24310848 (47 colonic biopsies, active vs non-inflamed; could not fetch, 403; accession and pairing [U], worth a manual check), Nottingham IBD stem/TA preprint Dec 2025 (epithelial focus, [U]).
Honest limit: I could not query SCP, CELLxGENE Census or dbGaP interactively this session; no 2025-26 open whole-biopsy paired cohort beyond GSE266616 surfaced. Search engine coverage of GEO is weak, so absence is not proof.

## Ranking (design rules: >= 6 paired patients, whole biopsy epi + stroma, counts public, median UMI >= 1,500, >= 50 cells/compartment)
| Rank | Cohort | Paired patients | Counts | Chemistry | Whole biopsy | Status for us |
|---|---|---|---|---|---|---|
| 1 | TAURUS GSE282122 | 28-29 (18 baseline; 20 same-site) [V earlier] | open, on disk | 10x v3/3.1 | yes | eligible, already seen, non-blind |
| 2 | Smillie SCP259 | ~18 UC [P] | SCP login | 10x v1/v2 | Epi + LP fractions of same biopsy | eligible, unseen by us; needs PI login |
| 3 | Eshghi SCP3755 | ~18 UC [P] | SCP login | BD Rhapsody | Epi + LP fractions | eligible, separate stratum, unseen |
| 4 | GSE266616 | 14 CD, same site, inflamed vs adjacent [V] | open | 10x v3.1 | epi + LP pooled | probably eligible, unseen; UMI and annotation to be checked; adjacency caveat |
| 5 | Martin GSE134809 | 11 ileal CD [V] | open | 10x v1/v2 | LP, stroma sparse | marginal (5-8 pass) |
| 6 | Kong SCP1884 | 5 TI + 3 colon [V] | SCP login | 10x | Epi + LP | too few |
| 7 | GSE331003 | ~9-13 CD resections [V] | open | 10x CITE | unclear stroma | small, own annotation |
| 8 | GSE288517 (bulk LCM) | ~10 UC | counts uncertain | bulk | two compartments | secondary |
Not usable: GSE190564, GSE206299, GSE253720, IBDverse, GSE266546, GSE214695, GSE250498 (no pairs / sorted / stroma depleted / too small).

## Pooled estimate of eligible paired patients
- Not yet seen by us (candidate blind/primary set): Smillie ~18 + Eshghi ~18 + GSE266616 14 = ~50 (SCP pair counts [P]/[U]; ~35 if SCP data fail checks). Open without login: only GSE266616 (14).
- Including TAURUS 28-29 (seen): ~78-80. Adding Martin/GSE331003 marginal: ~85-90 at best.
- Practical blind-from-scratch estimate with open data only: 14 patients (GSE266616), adequate vs the >= 6 rule but with adjacent-region caveat, unknown UMI and no annotation. Best design: Smillie (adjacent biopsies, same procedure) once the PI obtains the login.

## Bulk, brief (secondary)
On disk already (A55). New: GSE288517 (LCM epithelium + LP paired, ~10 UC, raw not on GEO); Dryad rfj6q57dw (TI paired, one compartment). Bulk needs >= 2 compartments per patient; only the LCM set qualifies.
Collection time: not recorded in any cohort checked [U]; the design relies on constant-offset cancellation.
