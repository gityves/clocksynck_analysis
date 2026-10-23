"""A30 gene classes, fixed before scoring (A30; `desync/validation/antiphase/method.md`). Mutually exclusive by priority (dict order).
Copied from the legacy `classes.py`; the Hallmark sets come from the MSigDB gmt under $SCCC_DATA/reference (v2024.1; the legacy
used v2026.1, the two cholesterol-homeostasis and mTORC1 sets used here are identical in both)."""

from __future__ import annotations

from common.paths import DATA

H = {}
for _f in sorted((DATA / "reference").glob("h.all.*.Hs.symbols.gmt"), reverse=True)[:1]:
    H = {l.split("\t")[0]: set(l.strip().split("\t")[2:]) for l in open(_f)}

FRAME = {"ARNTL", "BMAL1", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "RORA",
         "TEF", "CIART"}

# frozen D10 list: curated clock outputs (d10.CLOCK_OUT) + followers (>= 2 tissues), regenerated from d10.py
D10_32 = {"ADGRA2", "ADIRF-AS1", "ARHGAP42", "ARNTL2", "BHLHE40", "BHLHE41", "CDC25B", "CDK18", "CHST7", "CIPC",
          "DDIT4L", "EPB41L1", "FAM76A", "FOXO1", "GUCY1A1", "HLF", "JAM2", "KLF10", "LDB1", "LONRF3", "LPCAT3",
          "MAMDC2", "MPP2", "NAMPT", "NFIL3", "PPARA", "RAB3IL1", "RORC", "SUN2", "TSPAN4", "USP2", "WEE1"}

FASTING_UP = {"PDK4", "PDK2", "CPT1A", "CPT1B", "CPT2", "ACADVL", "ACADM", "SLC25A20", "ACOT1", "ACOT2", "UCP3",
              "HMGCS2", "ANGPTL4", "CD36", "PLIN2", "PLIN5", "LPIN1", "PNPLA2", "LIPE", "PPARGC1A", "IRS2", "TXNIP",
              "KLF15", "PCK1", "G6PC1", "IGFBP1", "FGF21", "FBXO32", "TRIM63", "FOXO1", "PPARA"}

# diag6 lists (../glucocorticoid/scripts/diag6_classes.py), copied verbatim
GC = {"FKBP5", "TSC22D3", "ZBTB16", "KLF9", "KLF15", "DDIT4", "SGK1", "DUSP1", "CDKN1C", "MT2A", "MT1X", "PDK4",
      "ERRFI1", "IRS2", "CEBPD", "SESN1", "FOXO3", "GLUL", "TXNIP", "PER1", "ALOX15B", "MAOA", "IL1R2", "CD163",
      "SLC19A2", "RHOB", "FAM107A", "SULT1A1"}
HEAT = {"HSPA1A", "HSPA1B", "HSPA8", "HSPH1", "HSP90AA1", "HSP90AB1", "DNAJB1", "DNAJA1", "HSPB1", "HSPD1", "HSPE1",
        "BAG3", "CIRBP", "RBM3", "SERPINH1", "CHORDC1", "HSPA4L", "HSPA6", "DNAJA4", "CRYAB", "FKBP4", "STIP1"}
ACTIVITY = {"NR4A1", "NR4A2", "NR4A3", "FOS", "FOSB", "JUN", "JUNB", "EGR1", "EGR2", "EGR3", "ATF3", "PPARGC1A", "MYC",
            "ANKRD1", "XIRP1", "HES1", "IER2", "IER3", "IER5", "MAFF", "KLF2", "KLF4", "CYR61", "CCN1", "CTGF", "CCN2",
            "GADD45B", "PPP1R15A", "BTG2", "ZFP36"}
FEEDING_CORE = {"SREBF1", "SREBF2", "FASN", "INSIG1", "LDLR", "HMGCR", "HMGCS1", "SQLE", "SCD", "ACACA", "ACLY",
                "ELOVL6", "ME1", "FADS1", "FADS2", "MLXIPL", "G0S2", "TKT", "ACSS2"}
FEEDING_UP = FEEDING_CORE | H["HALLMARK_CHOLESTEROL_HOMEOSTASIS"] | H["HALLMARK_MTORC1_SIGNALING"]

RAW = {"D10 clock outputs": D10_32, "fasting-up": FASTING_UP, "glucocorticoid": GC, "heat shock": HEAT,
       "activity / immediate-early": ACTIVITY, "feeding-up": FEEDING_UP}


def assign(raw=RAW, frame=FRAME):
    """gene -> class (first listing class wins); plane genes dropped."""
    out = {}
    for c, gs in raw.items():
        for g in sorted(gs - frame):
            out.setdefault(g, c)
    return out


CLASS_OF = assign()
# sensitivity sub-classes (labels on top of the exclusive assignment)
SUB = {"fasting-up minus GC overlap": {g for g, c in CLASS_OF.items() if c == "fasting-up"} - GC,
       "feeding-up core": {g for g, c in CLASS_OF.items() if c == "feeding-up"} & FEEDING_CORE}

