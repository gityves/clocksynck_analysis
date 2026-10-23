"""WP0 item 3: pre-registration scorecard.

Parses every amendment of the archive PLAN (A1-A108c) and RESULTS, and of `registered_experiments.md` /
`results_preliminary.md` (A109 onwards), into one row per (amendment, hypothesis id), or one row per amendment when
the plan names no hypothesis id.  Rules (deterministic; a field that no rule fills stays empty, with a `parse_note`):

* amendment units: markdown headings (any level) whose text starts with an amendment id (`A12`, `A31b`, `Amendment A5`,
  `A36-A38`; ranges expand).  Plan text of a unit runs to the next heading of equal or higher level.
* hypothesis ids: `H<n>[a-z]?(.<k>)?`, `R<n>.<k>`, `V<n>.<k>` in the plan text before any "run note"; mapped to the amendment
  of their number (else to the unit's first amendment).
* label: first plan line holding the id decides: "post hoc" > "primary" > "secondary" > "exploratory|descriptive".  Else
  "post hoc" when the heading says so; else empty.
* verdict source lines: result-note lines (RESULTS sections naming the amendment, or the plan text after "run note")
  that hold the hypothesis id (no id: lines holding "verdict"). Tokens: pass = pass|passes|passed|supported|holds|met;
  fail = fail|fails|failed|not met|not supported|not replicated|not confirmed; uninformative; infeasible = infeasible|not
  feasible|cannot be run.  One class over all lines -> verdict; several -> the line(s) containing "verdict" decide;
  still several -> empty ("ambiguous").  Estimate, interval, p and MDE are read from the line that gave the verdict.
* commit/date: first commit whose diff adds the heading text (`git log -S --reverse`), best effort.
"""
from __future__ import annotations
import re, subprocess
from pathlib import Path
import pandas as pd

ARCH = Path("/home/maxine/Documents/paychere/scCircadianCancer")
REPO = Path("/home/maxine/Documents/paychere/scCircadianOverview")
OUT = Path(__file__).resolve().parents[1] / "output"

AID = r"A(\d+)([a-z]?)"
HID = re.compile(r"(?<![\w.])([HRV]\d+[a-z]?\.\d+|H\d+[a-z]?)(?![\w]|\.\d)")
VERD = [("fail", r"\bfail(?:s|ed|ure)?\b|not met|not supported|not replicated|not confirmed"),
        ("uninformative", r"uninformative"),
        ("infeasible", r"infeasible|not feasible|cannot be run"),
        ("pass", r"\bpass(?:es|ed)?\b|\bsupported\b|\bholds\b|\bmet\b")]
NUM = r"[+\-−]?\d+(?:\.\d+)?(?:e-?\d+)?"
INTERVAL = re.compile(rf"({NUM})\s*[\[(]\s*({NUM})\s*(?:,|to|;)\s*({NUM})\s*[\])]")
PVAL = re.compile(rf"\bp(?:_\w+)?\s*(?:=|<|≈|~|:)?\s*(\d*\.?\d+(?:e-?\d+)?)")
MDE = re.compile(rf"MDE[^0-9\n]{{0,40}}({NUM})")


def sections(path):
    """[(level, heading, [lines])] of a markdown file; a section ends at the next heading of any level."""
    out, cur = [], None
    for ln in Path(path).read_text().splitlines():
        m = re.match(r"(#+)\s+(.*)", ln)
        if m:
            cur = (len(m.group(1)), m.group(2), [])
            out.append(cur)
        elif cur:
            cur[2].append(ln)
    return out


def expand(text):
    """Amendment ids named in a heading (ranges expanded); None if the heading is not an amendment unit."""
    head = re.sub(r"^(Amendments?|Note)\s+", "", text.strip())
    m = re.match(rf"{AID}(?:\s*[-–]\s*{AID})?", head)
    if not m:
        return None
    n1, s1, n2, s2 = m.groups()
    if n2 and int(n2) > int(n1):
        return [f"A{n}" for n in range(int(n1), int(n2) + 1)]
    return [f"A{n1}{s1}"] + ([f"A{n2}{s2}"] if n2 else [])


def ids_in_heading_any(text):
    """All amendment ids named in a RESULTS heading (anywhere in it)."""
    ids = set()
    for m in re.finditer(rf"{AID}(?:\s*[-–]\s*{AID})?", text):
        n1, s1, n2, s2 = m.groups()
        if n2 and int(n2) > int(n1):
            ids |= {f"A{n}" for n in range(int(n1), int(n2) + 1)}
        else:
            ids.add(f"A{n1}{s1}")
            if n2:
                ids.add(f"A{n2}{s2}")
    return ids


def amend_of(hid):
    m = re.match(r"[HRV](\d+)([a-z]?)", hid)
    return f"A{m.group(1)}{m.group(2)}"


def classify(line):
    if re.search(r"gates? (?:fail|did not pass)|failed gate|gate failed", line, re.I):
        return ["uninformative"]   # CLAUDE.md: a failed gate is uninformative
    if re.search(r"not tested|not run\b", line, re.I):
        return ["infeasible"]
    if re.search(r"partly|narrowly|nominal|mixed|indeterminate", line, re.I):
        return []   # qualified verdicts are not mapped to pass/fail
    got = [k for k, rx in VERD if re.search(rx, line, re.I)]
    if "fail" in got and "pass" in got:   # "not supported" contains "supported", "not met" contains "met"
        stripped = re.sub(r"not met|not supported|not replicated|not confirmed", " ", line, flags=re.I)
        got = [k for k, rx in VERD if re.search(rx, stripped, re.I)]
        if re.search(r"not met|not supported|not replicated|not confirmed", line, re.I):
            got.append("fail")
    return sorted(set(got))


def verdict_from(lines):
    cand = [(ln, classify(ln)) for ln in lines]
    cand = [(ln, c) for ln, c in cand if c]
    if not cand:
        return None, "", "no verdict token on source lines"
    allc = {x for _, c in cand for x in c}
    if len(allc) == 1:
        return allc.pop(), cand[0][0], ""
    vl = [(ln, c) for ln, c in cand if re.search("verdict", ln, re.I)]
    vc = {x for _, c in vl for x in c}
    if len(vc) == 1:
        return vc.pop(), vl[0][0], "several lines; the 'verdict' line decides"
    return None, cand[0][0], "ambiguous: " + "/".join(sorted(allc))


def keyword_pass(res_secs, a, hid, plan_lines=()):
    """Second pass: last verdict-bearing sentence (single class) of the amendment's RESULTS sections.
    Sections naming several amendments contribute only sentences naming `a`; with a hypothesis id, sentences naming it
    are tried first (else the amendment-level sentence is used, and flagged)."""
    sents = []
    for h, lines in res_secs:
        txt = re.sub(r"<!--.*?-->|[*`]", "", " ".join(l.strip() for l in lines))
        ss = re.split(r"(?<=[.;!?])\s+(?=[A-Z(\-0-9])", txt)
        if len(ids_in_heading_any(h)) > 1:
            ss = [x for x in ss if re.search(rf"\b{a}\b", x)]
        sents += ss
    def last(pool):
        for x in reversed(pool):
            c = classify(x)
            if len(c) == 1:
                return c[0], x
        return None, ""
    if hid:
        pat = re.compile(re.escape(hid) + r"(?![\w]|\.\d)")
        v, x = last([y for y in sents if pat.search(y)])
        if v:
            return v, x, ""
        hit = id_window_verdict(list(plan_lines) + [l for _, ls in res_secs for l in ls], hid)
        if hit:
            return hit[0], hit[1], "keyword: id-anchored verdict word"
        v, x = last(sents)
        return v, x, "keyword: amendment-level sentence, not hypothesis-specific" if v else ""
    v, x = last(sents)
    return v, x, ""


def id_window_verdict(lines, hid):
    """Third rule: a line in which the hypothesis id is followed within 25 characters by a verdict word (`H54.1 fails`,
    `H105.2 not met (uninformative)`).  Rule statements (`->`, `if`, `Reading`) are skipped; the last such line decides.
    'uninformative' written next to a fail word wins (the text labels it)."""
    pat = re.compile(re.escape(hid) + r"(?![\w]|\.\d)(.{0,25})", re.I)
    hit = None
    for ln in lines:
        if re.search(r"->|\bif\b|Reading|Supported if|called", ln, re.I):
            continue
        for m in pat.finditer(re.sub(r"[*`]", "", ln)):
            w = m.group(1)
            w = re.split(r"[;]|(?<=\w)\.\s", w)[0]
            c = classify(w)
            if "uninformative" in c:
                c = ["uninformative"]
            if len(c) == 1:
                hit = (c[0], ln.strip())
    return hit


def numbers(line):
    est = lo = hi = p = mde = None
    m = INTERVAL.search(line)
    if m:
        est, lo, hi = (x.replace("−", "-") for x in m.groups())
    pm = PVAL.search(line)
    p = pm.group(1) if pm else None
    mm = MDE.search(line)
    mde = mm.group(1) if mm else None
    return est, (f"[{lo}, {hi}]" if m else None), p, mde


def label_for(hid, plan_lines, head):
    defs = [ln for ln in plan_lines if hid and re.search(r"\*\*\s*" + re.escape(hid) + r"(?![\w]|\.\d)", ln)]
    mention = [ln for ln in plan_lines if hid and re.search(re.escape(hid) + r"(?![\w]|\.\d)", ln)]
    first = (defs + mention)[:1]
    for ln in first:
        if True:
            i = plan_lines.index(ln)
            ln = " ".join(plan_lines[i:i + 3])
            for key, rx in [("post hoc", r"post[- ]hoc"), ("primary", r"\bprimary\b"), ("secondary", r"\bsecondary\b"),
                            ("exploratory", r"exploratory|descriptive")]:
                if re.search(rx, ln, re.I):
                    return key, ""
            break
    if re.search(r"post[- ]hoc", head, re.I):
        return "post hoc", ""
    return "", "label not stated"


def git_first(repo, path, needle):
    needle = needle[:60]
    r = subprocess.run(["git", "-C", str(repo), "log", "--reverse", "-S" + needle, "--format=%h %as", "--", path],
                       capture_output=True, text=True)
    first = r.stdout.strip().splitlines()
    return tuple(first[0].split()) if first else ("", "")


def src_files(lines):
    txt = "\n".join(lines)
    hits = re.findall(r"src: ([^>]*?)\s*-->", txt) + re.findall(r"Outputs?: (`[^\n]*)", txt)
    return "; ".join(hits)[:300]


def errata():
    """amendment id -> headings of correction/erratum sections naming it (archive RESULTS and current results_preliminary)."""
    out = {}
    for path in (ARCH / "docs/RESULTS_preliminary.md", REPO / "paper/docs/results_preliminary.md"):
        for lvl, h, lines in sections(path):
            if re.search(r"correction|erratum|errata|^20\d\d-\d\d-\d\d", h, re.I):
                for a in ids_in_heading_any(h):
                    out.setdefault(a, []).append((h[:90], lines))
    return out


ERR = errata()


def build(plan_path, res_path, repo, plan_rel, res_rel, only=None):
    plan = sections(plan_path)
    res = sections(res_path)
    # result lines per amendment id (sections whose heading names the id)
    res_by = {}
    for lvl, h, lines in res:
        for a in ids_in_heading_any(h):
            res_by.setdefault(a, []).append((h, lines))
    rows = []
    for lvl, head, lines in plan:
        ids = expand(re.sub(r"\(.*", "", head)) if re.match(r"(Amendments?\s+)?A\d", head) else None
        if not ids or (only and not any(only(a) for a in ids)):
            continue
        # plan text: this section plus deeper sub-sections is already separate (sections end at any heading)
        k = next((i for i, ln in enumerate(lines) if re.search(r"run note|after running", ln, re.I)), len(lines))
        pre, note = lines[:k], lines[k:]
        hyp = {}
        for ln in pre:
            for m in HID.finditer(ln):
                hyp.setdefault(m.group(1), None)
        hdate = (re.search(r"\((\d{4}-\d\d-\d\d)", head) or [None, ""])[1]
        commit, cdate = git_first(repo, plan_rel, head)
        for a in ids:
            bold = {h for h in hyp if any(re.match(r"\s*(?:[-*]\s*)?\*\*\s*" + re.escape(h) + r"(?![\w]|\.\d)", ln) for ln in pre)}
            mine = [h for h in hyp if amend_of(h) == a or (len(ids) == 1 and h in bold)]
            res_secs = res_by.get(a, [])
            res_lines = [ln for _, ls in res_secs for ln in ls]
            res_commit = ""
            outf = src_files(res_lines)
            targets = mine or [""]
            for h in targets:
                notes = []
                if h:
                    pat = re.compile(re.escape(h) + r"(?![\w]|\.\d)")
                    src = [ln for ln in res_lines if pat.search(ln)] + [ln for ln in note if pat.search(ln)]
                else:
                    pat = re.compile(rf"\b{a}\b")
                    if len(ids) > 1:
                        src = [ln for ln in res_lines if pat.search(ln) and re.search("verdict|pass|fail", ln, re.I)]
                    else:
                        src = [ln for ln in res_lines if re.search("verdict", ln, re.I)] + \
                              [ln for ln in note if re.search("verdict", ln, re.I)]
                    notes.append("no hypothesis id in plan")
                if not res_lines and not note:
                    notes.append("no result note found")
                v, line, vn = verdict_from(src)
                vsrc, vtext, sup = ("structured", line, "") if v else ("", "", "")
                if not v:
                    k, x, kn = keyword_pass(res_secs, a, h, lines)
                    if k:
                        v, vsrc, vtext, line = k, "keyword", x, x
                        notes = [n for n in notes if not n.startswith(("no verdict token", "ambiguous"))] + ([kn] if kn else [])
                elif vn:
                    notes.append(vn)
                for eh, elines in ERR.get(a, []):
                    etxt = re.sub(r"<!--.*?-->|[*`]", "", " ".join(l.strip() for l in elines))
                    es = [x for x in re.split(r"(?<=[.;!?])\s+(?=[A-Z(\-0-9])", etxt) if len(classify(x)) == 1]
                    es = [x for x in es if re.search(r"\b(PASS|FAIL)\b", x)] or [x for x in es if re.search(rf"\b{a}\b", x)]
                    mentions = (not h) or re.search(re.escape(h) + r"(?![\w]|\.\d)", etxt) or eh.startswith("Correction of " + a)
                    if es and mentions:
                        v, vsrc, vtext, sup, line = classify(es[-1])[0], "erratum", es[-1], eh, es[-1]
                est, iv, p, mde = numbers(line) if line else (None, None, None, None)
                if a in ERR:
                    notes.append("erratum exists: " + " | ".join(h_ for h_, _ in ERR[a]))
                lab, ln_ = label_for(h, pre, head)
                if ln_:
                    notes.append(ln_)
                rows.append(dict(amendment=a, date=cdate or hdate, commit=commit, hypothesis_id=h, label=lab,
                                 estimate=est or "", interval=iv or "", p=p or "", MDE=mde or "", verdict=v or "",
                                 verdict_source=vsrc, verdict_text=vtext[:200], superseded_by=sup,
                                 output_file=outf, heading_date=hdate, source=plan_rel, parse_note="; ".join(n for n in notes if n)))
    return rows


def main():
    rows = build(ARCH / "docs/PLAN_internal_desync.md", ARCH / "docs/RESULTS_preliminary.md", ARCH,
                 "docs/PLAN_internal_desync.md", "docs/RESULTS_preliminary.md")
    rows += build(REPO / "paper/docs/registered_experiments.md", REPO / "paper/docs/results_preliminary.md", REPO,
                  "paper/docs/registered_experiments.md", "paper/docs/results_preliminary.md")
    df = pd.DataFrame(rows)
    # curated verdicts (hand-read run notes), source "curated"; they override parsed values
    cur = pd.read_csv(Path(__file__).with_name("scorecard_curation.csv"), keep_default_na=False)
    for _, c in cur.iterrows():
        m = (df.amendment == c.amendment) & (df.hypothesis_id == c.hypothesis_id)
        if not m.any() and (c.hypothesis_id == "") and (df.amendment == c.amendment).any():
            m = (df.amendment == c.amendment)   # amendment-level curation covers all its rows
        if not m.any():
            df = pd.concat([df, pd.DataFrame([dict(amendment=c.amendment, hypothesis_id=c.hypothesis_id, source="curated")])], ignore_index=True)
            m = (df.amendment == c.amendment) & (df.hypothesis_id == c.hypothesis_id)
        df.loc[m, ["verdict", "verdict_source", "verdict_text"]] = [c.verdict, "curated", c.quote[:200]]
        extra = f"curated ({c.results_line})" + (f"; {c.note}" if c.note else "")
        df.loc[m, "parse_note"] = df.loc[m, "parse_note"].fillna("").astype(str).str.cat([extra] * int(m.sum()), sep="; ")
        if c.superseded_by if "superseded_by" in c else False:
            df.loc[m, "superseded_by"] = c.superseded_by
    df = df.fillna("")
    has_h = set(df[df.hypothesis_id != ""].amendment)
    df = df[(df.hypothesis_id != "") | ~df.amendment.isin(has_h)]
    df = df.drop_duplicates(subset=["amendment", "hypothesis_id"], keep="first")
    cols = ["amendment", "date", "commit", "hypothesis_id", "label", "estimate", "interval", "p", "MDE", "verdict",
            "verdict_source", "verdict_text", "superseded_by", "output_file", "heading_date", "source", "parse_note"]
    df[cols].to_csv(OUT / "scorecard.csv", index=False)
    n_am = df.amendment.nunique()
    has = df.groupby("amendment").verdict.apply(lambda s: (s != "").any())
    print(f"amendments {n_am}; rows {len(df)}; amendments with >=1 parsed verdict {has.sum()} ({has.mean():.0%});"
          f" rows with verdict {(df.verdict != '').mean():.0%}")
    print("verdict counts", df.verdict.value_counts().to_dict())
    print("by source", df.verdict_source.value_counts().to_dict(), "superseded", df[df.superseded_by != ""][["amendment", "hypothesis_id", "verdict"]].values.tolist())
    print("no verdict:", " ".join(has[~has].index))


if __name__ == "__main__":
    main()
