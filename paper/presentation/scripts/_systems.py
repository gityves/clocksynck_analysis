"""Organ systems shared by the tissue figures (GTEx organ map, baboon atlas): groups, colours (seaborn default palette, fixed per
system, see the paper-figures skill) and the tissue -> system maps.

Baboon codes: Mure et al. 2018 (GSE98965) abbreviations; ocular tissues go with the brain (neural), the pineal and pituitary with
the endocrine glands (as GTEx's pituitary).
"""
import seaborn as sns

_TAB = sns.color_palette().as_hex()
COLOR = {
    "brain, eye": _TAB[6],
    "cardiovascular": _TAB[3],
    "digestive": _TAB[2],
    "adipose, skin, breast": _TAB[1],
    "endocrine": _TAB[4],
    "other": _TAB[8],
}
GTEX_PREFIX = {
    "brain, eye": ("Brain",),
    "cardiovascular": ("Artery", "Heart"),
    "digestive": ("Colon", "Esophagus", "Stomach", "Liver", "Pancreas", "Small Intestine"),
    "adipose, skin, breast": ("Adipose", "Skin", "Breast"),
    "endocrine": ("Adrenal", "Pituitary", "Thyroid"),
}
BABOON = {
    "brain, eye": {"AMY", "ARC", "CER", "DMH", "HAB", "HIP", "LGP", "LH", "MGP", "MMB", "OLB", "PON", "PRA", "PRC", "PUT", "PVN",
                   "SCN", "SON", "SUN", "THA", "VIC", "VMH", "COR", "ONH", "RET", "RPE", "IRI"},
    "cardiovascular": {"AOR", "HEA"},
    "digestive": {"ANT", "ASC", "CEC", "DEC", "DUO", "ILE", "LIV", "OES", "STF", "PAN"},
    "adipose, skin, breast": {"OMF", "WAM", "WAP", "WAR", "WAS", "WAT", "SKI"},
    "endocrine": {"ADC", "ADM", "PIT", "PIN", "THR"},
}


def gtex_system(tissue: str) -> str:
    return next((s for s, pre in GTEX_PREFIX.items() if tissue.startswith(pre)), "other")


def baboon_system(code: str) -> str:
    return next((s for s, codes in BABOON.items() if code in codes), "other")


def legend_row(ax, fontsize=10, step=None):
    """One-row legend of the systems on an empty axis (``ax`` is turned off)."""
    step = step or 1 / len(COLOR)
    for j, (s, c) in enumerate(COLOR.items()):
        x = 0.02 + j * step
        ax.scatter([x], [0.5], s=34, color=c, transform=ax.transAxes, clip_on=False)
        ax.annotate(s, (x, 0.5), xycoords="axes fraction", xytext=(8, 0), textcoords="offset points", va="center",
                    fontsize=fontsize, color="#222222")
    ax.axis("off")


def legend_in(ax, fontsize=9, loc="upper left"):
    """The systems as a frameless one-column legend inside ``ax`` (for a panel with free space in a corner)."""
    from matplotlib.lines import Line2D

    h = [Line2D([], [], ls="", marker="o", ms=6, color=c, label=s) for s, c in COLOR.items()]
    ax.legend(handles=h, loc=loc, frameon=False, fontsize=fontsize, handletextpad=0.1, borderaxespad=0.2, labelspacing=0.3)
