"""A126 step 2: assemble donor x compartment pseudobulks (panel genes) of Kumar (rebuilt by broad_cell_type), Reed (existing HBCA pseudobulk,
Census terms) and Murrow RM, and the donor tables. Output OUT/pb_all.parquet (study, donor, comp, n_cells, depth, genes), OUT/donors_all.csv."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd, anndata as ad
import a126_lib as L
OUT, D = L.OUT, L.DATA
panel = (OUT / "panel_genes.txt").read_text().split()
COMPS = ["hr", "lumprog", "basal", "fibro", "endo"]
rows, meta = [], []
# Kumar
c = pd.read_parquet(OUT / "kumar_pb_counts.parquet"); g = pd.read_csv(OUT / "kumar_pb_groups.csv", index_col=0)
c = c.rename(columns={"BMAL1": "ARNTL"}) if "ARNTL" not in c.columns else c  # Census symbol
mp = {"lumhr": "hr", "lumsec": "lumprog", "basal": "basal", "fibroblasts": "fibro", "vascular": "endo"}
g = g[g.broad.isin(mp)]
k = pd.concat([g[["donor", "n_cells", "depth"]].assign(study="kumar", comp=g.broad.map(mp)), c.loc[g.index]], axis=1)
rows.append(k)
kd = pd.read_csv(D / "atlas_wg/breast_kumar/donors.csv")
meta.append(pd.DataFrame(dict(study="kumar", donor=kd.donor_id, premeno=kd.donor_menopausal_status.eq("premenopausal"), post=kd.donor_menopausal_status.eq("postmenopausal"),
    age=np.nan, bmi=kd.donor_BMI_at_collection, parity=np.nan, brca=np.nan, hc=np.nan, site=kd.sample_source, stage=kd.development_stage)))
# Reed
a = ad.read_h5ad(D / "atlas_wg/breast_hbca/pseudobulk.h5ad")
sym = a.var.feature_name.astype(str).to_numpy()
pos = {}
for i, s in enumerate(sym): pos.setdefault(L.canonical(s), i)
cols = [(gname, pos[L.canonical(gname)]) for gname in panel if L.canonical(gname) in pos]
mpr = {"luminal hormone-sensing cell of mammary gland": "hr", "luminal adaptive secretory precursor cell of mammary gland": "lumprog",
       "basal-myoepithelial cell of mammary gland": "basal", "fibroblast of mammary gland": "fibro", "blood vessel endothelial cell": "endo"}
sel = a.obs.cell_type.isin(mpr).to_numpy()
X = a.X[sel]; X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
r = pd.DataFrame(X[:, [j for _, j in cols]], columns=[n for n, _ in cols], index=a.obs_names[sel])
o = a.obs[sel]
r = pd.concat([pd.DataFrame(dict(donor=o.donor_id.astype(str), n_cells=o.n_cells, depth=o.depth, study="reed", comp=o.cell_type.map(mpr))), r], axis=1)
assert not r.duplicated(["donor", "comp"]).any()
rows.append(r)
rd = pd.read_csv(D / "atlas_wg/breast_hbca/donors.csv")
meta.append(pd.DataFrame(dict(study="reed", donor=rd.donor_id, premeno=rd.donor_age <= 45, post=rd.donor_age >= 55, age=rd.donor_age, bmi=rd.BMI,
    parity=pd.to_numeric(rd.parity, errors="coerce"), brca=rd.brca_status, hc=np.nan, site=np.nan, stage=np.nan)))
# Murrow
mc = pd.read_parquet(OUT / "murrow_pb_counts.parquet"); mg = pd.read_csv(OUT / "murrow_pb_groups.csv", index_col=0)
mpm = {"HRpos_Luminal": "hr", "Secretory_Luminal": "lumprog", "Basal": "basal", "Fibroblast": "fibro", "Vascular_Endothelial": "endo"}
mg = mg[mg.type.isin(mpm)]
rows.append(pd.concat([mg[["donor", "n_cells", "depth"]].assign(study="murrow", comp=mg.type.map(mpm)), mc.loc[mg.index]], axis=1))
md = mg.groupby("donor").first()
meta.append(pd.DataFrame(dict(study="murrow", donor=md.index, premeno=True, post=False, age=md.age.to_numpy(), bmi=md.bmi.to_numpy(), parity=md.parity.to_numpy(),
    brca=np.nan, hc=md.hc.to_numpy(), site=np.nan, stage=md.obesity.astype(str).to_numpy())))
pb = pd.concat(rows, ignore_index=True)
common = [x for x in panel if all(x in df.columns for df in rows)]
pb = pb[["study", "donor", "comp", "n_cells", "depth"] + common]
pb["donor"] = pb.donor.astype(str); pb.to_parquet(OUT / "pb_all.parquet")
dm = pd.concat(meta, ignore_index=True); dm["donor"] = dm.donor.astype(str); dm.to_csv(OUT / "donors_all.csv", index=False)
print(pb.groupby(["study", "comp"]).size().unstack(), "\ncommon panel genes", len(common), "clock/luteal present", all(x in common for x in L.GENES + L.LUTEAL))
