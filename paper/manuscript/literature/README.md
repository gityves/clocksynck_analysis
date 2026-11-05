# Literature notes

One card per cited paper, read from the full text of the Zotero export (`clocksync` collection). The PDFs and the
extracted text are outside the repo (copyright):
`/mnt/max8tb/Paychere/literature/zotero_export/clocksync/{files,text}/`.

| file | content |
|---|---|
| `index.tsv` | bib key -> Zotero item, file, extracted characters, status (`missing` / `no file` = no full text) |
| `lit_index.py` | builds `index.tsv` and `text/<key>.txt` (`python -I lit_index.py <export> <references.bib> <index.tsv>`) |
| `check_quotes.py` | checks that every quote of every card occurs in that paper's text |
| `notes/<key>.md` | one card per paper (template below) |
| `synthesis.md` | claims matrix, methods table, comparison candidates, verification ideas, missing literature |

## Card template

```markdown
# <key>
**Ref:** First author et al., Journal year. DOI. Type: method | data | biology | clinical | review. Lab: <PI>.
**Source:** full text | abstract only

## Question and design
Organism, tissue / cell types, n, assay (bulk / sc / snRNA). **Sampling:** timed? time coverage (full 24 h, day
only, hospital hours, unknown)?

## Method (if any)
What is inferred, from what; assumptions (full-cycle coverage, fixed gene acrophases, cohort vs within person,
supervised training data); outputs.

## Key results
- Result with numbers. "verbatim quote" (section / figure)

## Not addressed / limits
- ...

## Link to us
- **Intro claim(s):** which statement of the introduction it supports.
- **Relation to clocksync:** background | supports | competitor (compare?) | data we could use.

## Ideas
- Checks we could run; methods to compare; references it cites that we lack.
```

Quotes are verbatim, in straight double quotes, at least 20 characters; `check_quotes.py` must report none missing
before a card is used.
