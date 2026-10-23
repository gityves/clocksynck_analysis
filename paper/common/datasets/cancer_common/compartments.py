"""Cells -> compartments from marker genes, for the single-cell cancer cohorts.

A cell gets the compartment with the highest mean log-normalised marker expression when it clearly wins (score above
`min_score` and second best below `margin` x best); otherwise `unassigned`. Markers are canonical and tissue-agnostic.
Malignant cells are not a compartment here: in tumour samples they sit inside `epithelial` until `cnv.study_calls`
splits them. Ported unchanged from legacy/src/sccc/prep/compartments.py (the Cell Ontology rules there are unused).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from clocksync.genes import gene_columns


# --- From marker genes --------------------------------------------------------------------------------------

MARKERS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT5", "KRT14", "CLDN4"],
    "fibroblast": ["DCN", "LUM", "COL1A1", "COL1A2", "PDGFRA", "FBLN1", "COL3A1"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "EGFL7", "FLT1"],
    "mural": ["RGS5", "PDGFRB", "MYH11", "ACTA2", "TAGLN", "KCNJ8"],
    "T_NK": ["CD3E", "CD3D", "CD2", "NKG7", "GNLY", "IL7R"],
    "B_plasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGKC"],
    "myeloid": ["CD68", "LYZ", "CD14", "C1QA", "C1QB", "AIF1", "FCER1G"],
    "mast": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"],
}


def marker_compartments(X, var_names, total_umi, markers=MARKERS, min_score: float = 0.3,
                        margin: float = 0.7) -> np.ndarray:
    names = list(markers)
    X = sp.csr_matrix(X)
    scale = 1e4 / np.maximum(np.asarray(total_umi, float), 1.0)
    cols = []
    for n in names:
        c = gene_columns(var_names, markers[n])
        c = c[c >= 0]
        M = X[:, c].multiply(scale[:, None]).tocsr()
        M.data = np.log1p(M.data)
        cols.append(np.asarray(M.mean(1)).ravel())
    S = np.column_stack(cols)
    order = np.argsort(-S, 1)
    r = np.arange(len(S))
    best, second = S[r, order[:, 0]], S[r, order[:, 1]]
    return np.where((best > min_score) & (second < margin * best), np.array(names)[order[:, 0]], "unassigned")
