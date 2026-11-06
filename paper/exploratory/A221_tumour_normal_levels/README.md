# A221: clock-gene levels, tumour vs adjacent normal vs GTEx

Registration: `paper/docs/registered_experiments.md` (A221). Code: `scripts/compute_a221.py`.

    RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A221_tumour_normal_levels/scripts/compute_a221.py --stage power --jobs 4
    ... --stage tests   (registered; run only after the power note is dated and committed)
    ... --smoke --out DIR   (KICH + UCEC, 20 sets)

## Data
GTEx recount3 G029 for the matched tissues (`$SCCC_DATA/gtex_recount3/`, PROVENANCE.md, MANIFEST.tsv; 0.54 GB raw, 1.2 GB
with the cache `cache_a221/`). TCGA recount3 as on disk. 13 projects: BRCA COAD ESCA KICH KIRC KIRP LIHC LUAD LUSC PRAD STAD THCA UCEC.
Time of death (labelled h5ad `tod`, by subject) and Hardy class are held in memory; only aggregate counts and statistics are written.

## Power stage (what is computed)
eligibility.csv (pairs, tumours, GTEx samples, samples with time of death, daytime 08-18 h count); 200 random sets per project;
null draws of L and Delta from the random sets only; power.csv (null SD, MDE = 2.49 x SD). No clock-gene tumour-vs-normal
difference is computed in this stage.

## Choices left open by the registration
- Counts: recount3 coverage sums, CPM = coverage / summed protein-coding coverage x 1e6 (AUC scaling cancels), log2(CPM + 1); duplicated symbols summed.
- Tumours = one primary tumour per patient (lowest barcode, all histologies); normals = solid-tissue normals (as `tcga_paired`); pairs = patients with both. LFC_gtex uses all such tumours.
- GTEx tissue by SMTSD (colon transverse, oesophagus mucosa, kidney cortex); the three kidney projects (and the two lung projects) share one GTEx sample set, so projects are not independent in H221.2/3.
- Candidate random genes: protein-coding, mean CPM >= 1 over all tumours and normals, outside the clock set; deciles of the mean log2 CPM over those samples; each clock gene is replaced by a distinct gene of its decile (decile counts of the clock set reproduced).
- Null draw of set k: its mean over projects minus the mean of the 200 sets in each project (the reference of L). Delta (tests) is read against the raw random-set values, two-sided around their mean; null SD is identical. Permutation p = (1 + count) / (1 + n).
- H221.3 daytime = time of death in [8, 18) h; offset over donors with a time of death; Welch p in the descriptive: tumour vs normal groups (paired patients) and tumours vs GTEx.
- Hardy descriptive: classes with >= 10 samples in the tissue; purity: ABSOLUTE `tcga_purity_absolute/mastercalls.txt`, Spearman with mean clock log2 CPM within tumours.
- Tests stage is implemented but untested on real data (helpers checked on synthetic data only).

## Runtime
Power stage 36 s, peak 2.9 GB (4 workers).
