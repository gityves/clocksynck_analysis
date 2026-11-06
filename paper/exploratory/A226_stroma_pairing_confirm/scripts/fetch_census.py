"""A226: download the source h5ad of Census datasets (metadata-selected in scout/census_paired_candidates.csv). Parse only."""
import sys, hashlib
import cellxgene_census as cc
from pathlib import Path
root = Path(sys.argv[1])
for pre in sys.argv[2:]:
    with cc.open_soma(census_version="2025-11-08") as C:
        ds = C["census_info"]["datasets"].read().concat().to_pandas()
    row = ds[ds.dataset_id.str.startswith(pre)].iloc[0]
    d = root / ("cellxgene_" + pre); d.mkdir(parents=True, exist_ok=True)
    f = d / f"{pre}.h5ad"
    cc.download_source_h5ad(row.dataset_id, to_path=str(f), census_version="2025-11-08")
    h = hashlib.sha256(f.read_bytes()).hexdigest()
    (d / "PROVENANCE.md").write_text(f"# Provenance\n- URL: CELLxGENE Census 2025-11-08 download_source_h5ad, dataset {row.dataset_id} ({row.collection_name}; {row.dataset_title}; doi {row.collection_doi})\n- Downloaded: 2026-11-06\n- sha256: {h} ({f.name}, {f.stat().st_size/1e6:.0f} MB)\n- Parse only (python -I), not executed. A226.\n")
    print(pre, f.stat().st_size / 1e6, flush=True)
