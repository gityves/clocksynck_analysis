"""Check that every quote of every literature card occurs in that paper's extracted full text.

Usage: python -I check_quotes.py <text dir> <notes dir> [key ...]
A quote is any straight-double-quoted string of >= 20 characters. Matching ignores case, whitespace, line-break
hyphenation, ligatures and typographic quote/dash variants. Exit status 1 if any quote is missing.
"""
import re
import sys
import unicodedata
from pathlib import Path

SUB = {"ﬁ": "fi", "ﬂ": "fl", "ﬀ": "ff", "ﬃ": "ffi", "ﬄ": "ffl", "’": "'", "‘": "'",
       "“": '"', "”": '"', "–": "-", "—": "-", "−": "-", " ": " ", " ": " "}


def norm(s):
    s = unicodedata.normalize("NFKC", s)
    for a, b in SUB.items():
        s = s.replace(a, b)
    s = re.sub(r"-\s*\n\s*", "", s)  # hyphenation at line breaks
    return re.sub(r"[^a-z0-9]", "", s.lower())  # letters and digits only


def main(text_dir, notes_dir, keys):
    cards = sorted(Path(notes_dir).glob("*.md"))
    if keys:
        cards = [c for c in cards if c.stem in keys]
    bad = 0
    for card in cards:
        txt = Path(text_dir) / f"{card.stem}.txt"
        body = card.read_text()
        quotes = [q for q in re.findall(r'"([^"\n]{20,})"', body)]
        if not txt.exists():
            print(f"{card.stem}: no text ({len(quotes)} quotes unchecked)")
            continue
        full = norm(txt.read_text(errors="replace"))
        miss = [q for q in quotes if norm(q.replace("...", "")) not in full and
                not all(norm(p) in full for p in q.split("...") if p.strip())]
        bad += len(miss)
        print(f"{card.stem}: {len(quotes) - len(miss)}/{len(quotes)} quotes found")
        for q in miss:
            print(f"    MISSING: {q[:120]}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], set(sys.argv[3:]))
