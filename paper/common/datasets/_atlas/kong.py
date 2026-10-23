"""Kong et al. 2023 gut atlas, normal terminal ileum and colon: shared builder of `kong_ileum` and `kong_colon`.

Legacy: desync/18_gut_kong/scripts/kong.py (`SEGMENTS`, `table`, `run`); A17 targets `kong_ileum` (A16 primary) and
`kong_colon` (A16 secondary), healthy fits.
Source: six Census datasets (epithelium, immune and stromal fractions of ileum and colon), raw inputs
$SCCC_DATA/atlas/{kong_ti_epi,ileum_kong,kong_ti_str} and {kong_co_epi,kong_co_imm,kong_co_str}.

Choices: normal groups of the segment's sample tissues; classes `COARSE`; compartment `<segment>:<class>`; all files
share the donor key `kong:<donor_id>`; the majority assay is taken per file; pseudobulks >= 20 cells; classes with >= 15
donors. All donors are `healthy`.
"""

from __future__ import annotations

import pandas as pd

import census as C
import legacy as L
import organ as O

SEGMENTS = {  # segment: (datasets, sample tissues)
    "ileum": (("kong_ti_epi", "ileum_kong", "kong_ti_str"), ("ileal epithelium", "ileum lamina propria")),
    "colon": (("kong_co_epi", "kong_co_imm", "kong_co_str"),
              ("colonic epithelium", "lamina propria of mucosa of colon")),
}


def specs(seg: str):
    dss, tissues = SEGMENTS[seg]
    keep = lambda o: o.tissue.astype(str).isin(tissues).to_numpy() & (o.disease.astype(str) == "normal").to_numpy()  # noqa: E731,E501
    return [O.Spec(C.ATLAS / d, "kong", seg, O.COARSE, keep=keep) for d in dss]


def _donors(seg: str) -> pd.DataFrame:
    tabs = [C.census_donors(s.path, "kong") for s in specs(seg)]
    return pd.concat(tabs).loc[lambda t: ~t.index.duplicated()]


def build(name: str, donors=None, genes="clock"):
    seg = name.removeprefix("kong_")
    F = O.concat_groups(specs(seg), donors, genes)
    pb = O.table(F, F.resection.unique(), 20, min_donors=None if donors is not None else 15)
    return O.finish(pb, _donors(seg), {"healthy": O.flag(pb.resection.unique())})


def from_legacy(name: str):
    d = _donors(name.removeprefix("kong_"))
    return L.from_pickles([name], covariates=C.donors_frame(d, d.index))
