"""Design of A266 (fetal adrenal cortex: entrained loop or driven response?), shared by the feasibility and the scoring.

Rows and covariates: the `cortex` readout of A244 (`design_a264`: cohort indicator, common quadratic age). Waveforms:
  CORE   the 14-gene GTEx waveform without the clock genes with a documented direct acute response to glucocorticoid or
         cAMP signalling: PER1 (GRE and CRE), PER2 (GRE) and NR1D1 (glucocorticoid repression); 11 genes.
Responsive genes of H266.2: `clocksync._core.controls.RESPONSIVE_GENES` (glucocorticoid targets and immediate-early genes,
the pool of the responsive companion) with a mean count >= MIN_COUNT in the fetal cortex pseudobulks and present in GTEx.
Receptor genes of H266.3 (descriptive): MTNR1A, MTNR1B (melatonin), MC2R, MRAP (ACTH), NR3C1 (glucocorticoid), HSD11B2.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[2] / "A264_fetal_adrenal_covariates" / "scripts"))
sys.path.insert(0, str(HERE.parents[2] / "A244_fetal_adrenal_classes" / "scripts"))
import design_a264 as d264  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 266
WF = d264.WF
SIGNAL_INDUCED = ("PER1", "PER2", "NR1D1")
CORE = WF.drop(list(SIGNAL_INDUCED))
MIN_COUNT = 20.0
RECEPTORS = ("MTNR1A", "MTNR1B", "MC2R", "MRAP", "NR3C1", "HSD11B2")


def tests(obs) -> dict:
    """{name: (rows, covariates, strata, waveform)}: H266.1 and its descriptive companions."""
    rows = d264.d244.pooled(obs, obs.compartment == "cortex")
    return {"cortex:core": (rows, d264.cov(), "block", CORE),
            "cortex:core:steroid": (rows, d264.cov(score="steroid"), "block", CORE),
            "medulla:core": (d264.d244.pooled(obs, obs.compartment == "medulla"), d264.cov(), "block", CORE)}
