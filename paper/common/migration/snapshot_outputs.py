"""Copy the archive's full-run outputs of the migrated folders (outputs/paper/<folder>/, not subsample/, not runs/)
to results/<figure>/output/<group>/, so that scripts reading another group's tables find them before a rerun.

Every copied file is listed with its sha256 in results/<figure>/output/<group>/ARCHIVE_SNAPSHOT.tsv; a file whose hash
still matches is an archive number, not a rerun. Tables with a per-donor column (donor, subject, SUBJID) are not
copied (GTEx controlled access). Usage: python -I snapshot_outputs.py
"""

from __future__ import annotations

import gzip
import hashlib
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from migrate_from_archive import ARCHIVE, FOLDERS, RESULTS  # noqa: E402

SRC = ARCHIVE / "outputs" / "paper"
PER_DONOR = {"donor", "subject", "subjid", "subj.id", "sample_id"}


def header(f: Path) -> set[str]:
    if f.suffix not in (".csv", ".tsv", ".gz"):
        return set()
    opener = gzip.open if f.suffix == ".gz" else open
    try:
        with opener(f, "rt") as fh:
            line = fh.readline()
    except (OSError, UnicodeDecodeError):
        return set()
    return {c.strip().strip('"').lower() for c in line.replace("\t", ",").split(",")}


def main():
    for folder, (fig, group) in FOLDERS.items():
        src = SRC / folder
        if not src.exists():
            continue
        dst = RESULTS / fig / "output" / group
        rows, skipped = [], []
        for f in sorted(src.rglob("*")):
            rel = f.relative_to(src)
            if not f.is_file() or {"subsample", "runs", "cache"} & set(rel.parts):
                continue
            if header(f) & PER_DONOR:
                skipped.append(str(rel))
                continue
            (dst / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dst / rel)
            rows.append(f"{rel}\t{hashlib.sha256(f.read_bytes()).hexdigest()}")
        if rows:
            (dst / "ARCHIVE_SNAPSHOT.tsv").write_text(
                f"# copied from {src} (archive outputs, not rerun here)\nfile\tsha256\n" + "\n".join(rows) + "\n")
        print(f"{folder} -> {fig}/output/{group}: {len(rows)} files; per-donor tables skipped: {skipped}")


if __name__ == "__main__":
    main()
