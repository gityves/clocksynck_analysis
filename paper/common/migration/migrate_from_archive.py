"""One-off migration of the clean analyses of the archive (scCircadianCancer/paper/code) into this repository.

Copies the analysis folders cited by paper/docs/provisory_plan_paper.md to results/<figure>/scripts/<group>/, the
dataset builders to common/datasets/ and the shared helpers to common/, and rewrites every path expression so that
it resolves to the same file as in the archive:

- code paths (`Path(__file__)...`, `HERE...`, `CODE / ...`) -> their new place (`PAPER`, `DATASETS`, `group_dir()`);
- archive files outside the code (legacy tables, outputs) -> `ARCHIVE / ...`; the archive `data/` -> `PAPER_DATA`;
- outputs of another folder (`OUTPUTS / "figure_X"`) -> `group_out("<group>")`;
- folder names (`FOLDER = "figure_S07_..."`) -> group names.

Usage: python -I migrate_from_archive.py [--dry-run]. Writes a report of every rewrite to migration_report.md.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
from pathlib import Path

ARCHIVE = Path("/home/maxine/Documents/paychere/scCircadianCancer")
OLD_CODE = ARCHIVE / "paper" / "code"
PAPER = Path(__file__).resolve().parents[2]
DATASETS = PAPER / "common" / "datasets"
RESULTS = PAPER / "results"

FOLDERS = {  # archive folder -> (figure, group)
    "figure_1_method": ("figure_1", "method"),
    "figure_S02_windows": ("figure_1", "windows"),
    "figure_S12_species": ("figure_1", "species"),
    "figure_2_organs": ("figure_2", "organs"),
    "figure_3_backbone": ("figure_2", "backbone"),
    "figure_3_embryo": ("figure_2", "embryo"),
    "figure_S05_atlas_fits": ("figure_2", "atlas_fits"),
    "figure_S08_gene_rhythms": ("figure_2", "gene_rhythms"),
    "figure_4_cues_lost": ("figure_3", "cues_lost"),
    "figure_S07_gtex_covariates_cues": ("figure_3", "gtex_cues"),
    "figure_5_resilience": ("figure_4", "resilience"),
    "figure_S06_spread_multiplicity": ("figure_4", "spread_multiplicity"),
    "figure_S11_natural_experiments": ("figure_4", "natural_experiments"),
    "figure_6_cancer": ("figure_5", "cancer"),
    "figure_S10_cancer_details": ("figure_5", "cancer_pipeline"),
    "figure_6_lung_breaks": ("figure_5", "lung_breaks"),
    "figure_7_psoriasis": ("figure_5", "psoriasis"),
    "figure_S09_psoriasis_details": ("figure_5", "psoriasis_details"),
    "figure_S01_calibration": ("extended_data", "calibration"),
    "figure_S03_gate": ("extended_data", "gate"),
    "figure_S04_host_identity": ("extended_data", "host_identity"),
    "figure_S13_simulations": ("extended_data", "simulations"),
}
# Not cited by the plan (stay in the archive): A5, A6, A7, A12, A13, A29, A30, A32 and the 2026-10-04 cycling checks.
EXCLUDE = {
    "figure_6_lung_breaks": {"compute_a5_ild.py", "compute_a6_stress.py", "compute_a6_tstates.py",
                             "compute_a7_activation.py", "compute_a13.py", "_stress.py"},
    "figure_S04_host_identity": {"compute_axis_reading.py", "compute_common_axis.py"},
    "figure_6_cancer": {"compute_cycling_cophaser.py", "compute_cycling_scores.py"},
    "figure_S08_gene_rhythms": {"compute_bulk_clock.py", "compute_antiphase.py", "compute_antiphase_modes.py",
                                "compute_blood_gtex_arms.py"},
}
EXCLUDE_DATASETS = {"hlca_tstates_a7", "natri_tstates_a7"}  # A7 only
COMMON_FILES = ["__init__.py", "cli.py", "data.py", "gtex.py", "anchors.py", "anchors_a17.csv", "benchmark"]

BASE = r"(?:Path\(__file__\)\.resolve\(\)|C\.HERE|(?<![\w.])HERE|(?<![\w.])CODE|(?<![\w.])REPO)"
SUFFIX = r"(?:\.parent(?!s)|\.parents\[\d\])*"
LITS = r"(?:\s*/\s*\"[^\"{}]*\")*"
CHAIN = re.compile(f"({BASE})({SUFFIX})({LITS})")
LIT = re.compile(r"\"([^\"]*)\"")
NAMES = ("PAPER", "DATASETS", "PAPER_DATA", "ARCHIVE", "group_dir", "group_out")


def new_location(old: Path) -> Path:
    """New path of an archive code file or folder."""
    rel = old.relative_to(OLD_CODE)
    if not rel.parts:
        return PAPER
    head, rest = rel.parts[0], rel.parts[1:]
    if head == "data":
        return DATASETS.joinpath(*rest)
    if head == "common":
        return PAPER.joinpath("common", *rest)
    fig, group = FOLDERS[head]
    return RESULTS.joinpath(fig, "scripts", group, *rest)


def evaluate(base: str, suffix: str, lits: list[str], file_old: Path, here_old: Path | None,
             c_here_old: Path) -> Path:
    p = {"Path(__file__).resolve()": file_old, "HERE": here_old, "C.HERE": c_here_old,
         "CODE": OLD_CODE, "REPO": ARCHIVE}[base]
    if p is None:
        raise ValueError("HERE used but not defined")
    for m in re.finditer(r"\.parent(?!s)|\.parents\[(\d)\]", suffix):
        p = p.parent if m.group(1) is None else p.parents[int(m.group(1))]
    return p.joinpath(*lits)


def express(target_old: Path, file_new: Path, bootstrap: bool) -> tuple[str, set[str]]:
    """New expression for an archive path, and the common.paths names it needs."""
    def join(head, parts):
        return head + "".join(f' / "{x}"' for x in parts)

    try:
        rel = target_old.relative_to(OLD_CODE)
    except ValueError:
        rel = None
    if rel is not None:
        parts = rel.parts
        if not parts:
            if bootstrap:
                n = len(file_new.relative_to(PAPER).parts) - 1
                return f"Path(__file__).resolve().parents[{n}]", set()
            return "PAPER", {"PAPER"}
        if parts[0] == "data":
            return join("DATASETS", parts[1:]), {"DATASETS"}
        if parts[0] == "common":
            return join("PAPER", parts), {"PAPER"}
        if parts[0] in FOLDERS:
            return join(f'group_dir("{FOLDERS[parts[0]][1]}")', parts[1:]), {"group_dir"}
        raise ValueError(f"unmapped code path {target_old}")
    rel = target_old.relative_to(ARCHIVE)
    if rel.parts[:1] == ("data",):
        return join("PAPER_DATA", rel.parts[1:]), {"PAPER_DATA"}
    return join("ARCHIVE", rel.parts), {"ARCHIVE"}


def rewrite(text: str, file_old: Path, file_new: Path, log: list[str]) -> str:
    here_old = None
    m = re.search(r"^HERE = (Path\(__file__\)\.resolve\(\)(?:\.parent(?!s)|\.parents\[\d\])*)\s*$", text, re.M)
    if m:
        here_old = evaluate("Path(__file__).resolve()", m.group(1)[len("Path(__file__).resolve()"):], [],
                            file_old, None, OLD_CODE / "data" / "_atlas")
    c_here_old = OLD_CODE / "data" / "_atlas"  # `import census as C` in data/_atlas
    needed: set[str] = set()
    # A local `CODE = <path>` (the archive code root) becomes the new paper root, put on sys.path at once.
    n_up = len(file_new.relative_to(PAPER).parts) - 1
    text, n_local = re.subn(rf"^CODE = {CHAIN.pattern}\s*$",
                            f"PAPER = Path(__file__).resolve().parents[{n_up}]\nsys.path.insert(0, str(PAPER))",
                            text, flags=re.M)
    if n_local:
        log.append("local `CODE = ...` -> `PAPER = ...` (put on sys.path)")
    out, pos = [], 0
    for m in CHAIN.finditer(text):
        base, suffix, lits_s = m.group(1), m.group(2), m.group(3)
        lits = LIT.findall(lits_s)
        if base in ("HERE", "C.HERE") and here_old is None and base == "HERE":
            continue
        if base == "Path(__file__).resolve()" and not suffix:
            continue  # the file itself
        line_start = text.rfind("\n", 0, m.start()) + 1
        line = text[line_start:text.find("\n", m.end()) if text.find("\n", m.end()) >= 0 else len(text)]
        if line.lstrip().startswith(("#", "from common.paths import", "PAPER = Path(__file__)")):
            continue
        old_target = evaluate(base, suffix, lits, file_old, here_old, c_here_old)
        if base in ("Path(__file__).resolve()", "HERE", "C.HERE"):
            same = evaluate(base, suffix, lits, file_new,
                            None if here_old is None else new_location(here_old) if _in_code(here_old) else here_old,
                            DATASETS / "_atlas")
            if _in_code(old_target) and new_location(old_target) == same:
                continue
        bootstrap = "sys.path" in line and old_target == OLD_CODE
        expr, names = express(old_target, file_new, bootstrap)
        needed |= names
        out.append(text[pos:m.start()] + expr)
        pos = m.end()
        log.append(f"`{m.group(0).strip()}` -> `{expr}`")
    text = "".join(out) + text[pos:]

    def outputs(mo):
        arg = mo.group(1)
        if arg.startswith('"'):
            folder = arg.strip('"')
            arg = f'"{FOLDERS[folder][1]}"' if folder in FOLDERS else arg
        log.append(f"`OUTPUTS / {mo.group(1)}` -> `group_out({arg})`")
        return f"group_out({arg})"

    text, n = re.subn(r"(?<![\w.])OUTPUTS\s*/\s*(\"[^\"]+\"|[A-Za-z_][\w.]*)", outputs, text)
    if n:
        needed.add("group_out")
    return _fix_imports(text, needed, log).replace("__N__", str(n_up))


def _in_code(p: Path) -> bool:
    try:
        p.relative_to(OLD_CODE)
        return True
    except ValueError:
        return False


def _fix_imports(text: str, needed: set[str], log: list[str]) -> str:
    """Drop the removed names (CODE, OUTPUTS, REPO) from `from common.paths import`, and import the needed names
    right after the line that puts the paper root on sys.path (or after the docstring of an imported module)."""
    imp = re.compile(r"^from common\.paths import (\([^)]*\)|[^\n]*)", re.M)

    def clean(m):
        names = set(re.sub(r"#.*|[()]", " ", m.group(1)).replace(",", " ").split())
        names -= {"CODE", "OUTPUTS", "REPO", "FIGURES"}
        tail = "  # noqa: E402" if "noqa" in m.group(0) else ""
        return f"from common.paths import {', '.join(sorted(names, key=str.lower))}{tail}" if names else ""

    text = imp.sub(clean, text)
    if not needed:
        return text
    # insert before the first code line that uses a needed name, with a paper-root bootstrap when possible
    doc = re.match(r'\s*("""[\s\S]*?"""\s*\n)?', text)
    body_start = doc.end()
    use = re.compile(r"^(?!\s*#)[^\n]*\b(" + "|".join(sorted(needed)) + r")\b", re.M)
    m = use.search(text, body_start)
    at = m.start() if m else len(text)
    indent = re.match(r"\s*", text[at:]).group(0)
    if indent.strip("\n"):
        log.append("WARNING: first use is indented; import placed at module level before it")
        at = text.rfind("\n", 0, text.rfind("\n", 0, at)) + 1 if False else at - len(indent.lstrip("\n"))
        at = text.rfind("\n", 0, at) + 1
        prev = re.compile(r"^\S", re.M)
        k = [mm.start() for mm in prev.finditer(text, body_start, at)]
        at = k[-1] if k else at
    head = text[:at]
    lines = ""
    if re.search(r"^import sys\b", head, re.M) and re.search(r"^from pathlib import [^\n]*\bPath\b", head, re.M):
        n_up = None  # filled by the caller
        lines += "sys.path.insert(0, str(Path(__file__).resolve().parents[__N__]))  # paper/\n"
    lines += f"from common.paths import {', '.join(sorted(needed, key=str.lower))}  # noqa: E402\n"
    return text[:at] + lines + text[at:]


def rename_folders(text: str) -> str:
    for old, (fig, group) in sorted(FOLDERS.items(), key=lambda kv: -len(kv[0])):
        text = text.replace(f"outputs/paper/{old}", f"paper/results/{fig}/output/{group}")
        text = text.replace(f"paper/code/{old}", f"paper/results/{fig}/scripts/{group}")
        text = re.sub(rf"\bfrom {old}([ .])", rf"from results.{fig}.scripts.{group}\1", text)
        text = re.sub(rf"(?<![\w/]){old}\b", group, text)
    text = text.replace("paper/code/data/", "paper/common/datasets/").replace("paper/code/common", "paper/common")
    text = text.replace("docs/PLAN_internal_desync.md", "docs/PLAN_internal_desync.md (archive)")
    return text


def migrate_file(src: Path, dst: Path, log: list[str], dry: bool):
    dst.parent.mkdir(parents=True, exist_ok=True) if not dry else None
    if src.suffix in (".py", ".sh", ".md", ".jl", ".R", ".txt") and src.name != "anchors_a17.csv":
        text = src.read_text()
        entry: list[str] = []
        if src.suffix == ".py":
            text = rewrite(text, src, dst, entry)
        text = rename_folders(text)
        if src.name == "run.sh":
            excl = EXCLUDE.get(src.parent.name, set())
            text = "\n".join(l for l in text.split("\n") if not any(x in l for x in excl))
        if entry:
            log.append(f"### {dst.relative_to(PAPER)}\n" + "\n".join(f"- {e}" for e in entry))
        if not dry:
            dst.write_text(text)
            shutil.copymode(src, dst)
    elif not dry:
        shutil.copy2(src, dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    dry = ap.parse_args().dry_run
    log: list[str] = []
    jobs = []
    for name in COMMON_FILES:
        jobs.append(OLD_CODE / "common" / name)
    for d in sorted((OLD_CODE / "data").iterdir()):
        if d.is_dir() and d.name not in EXCLUDE_DATASETS and d.name != "__pycache__":
            jobs.append(d)
    for f in sorted((OLD_CODE / "data").glob("INDEX_*.md")):
        jobs.append(f)
    for folder in FOLDERS:
        for f in sorted((OLD_CODE / folder).iterdir()):
            if f.name not in EXCLUDE.get(folder, set()) and f.name != "__pycache__":
                jobs.append(f)
    for job in jobs:
        files = [job] if job.is_file() else [f for f in job.rglob("*") if f.is_file() and "__pycache__" not in f.parts]
        for f in files:
            migrate_file(f, new_location(f), log, dry)
    report = Path(__file__).with_name("migration_report.md")
    report.write_text("# Path rewrites of the migration (generated by migrate_from_archive.py)\n\n" + "\n\n".join(log)
                      + "\n")
    print(f"{len(jobs)} items; {sum(1 for l in log if l.startswith('###'))} files rewritten; report {report}")


if __name__ == "__main__":
    main()
