---
name: paper-figures
description: House style and procedure for every figure of the sccc paper and presentation (matplotlib + paperfig). Use before writing or editing any plotting script in paper/presentation/scripts or paper/results/*/scripts, and when the user comments on a figure.
---

# Paper and slide figures

A living document: when the user corrects a figure, add the rule here (with the date) so it holds for every later
figure, then apply it to the figure at hand.

## Tools

- Layout and checks: `paperfig` (the user's package, `/home/maxine/Documents/paychere/paperfig`, editable install in the
  `sccc` env). Slides: `pf.use_style("slide")` + `pf.Layout.slide(mosaic, kind="half" | "wide" | "third")`; journal:
  `pf.use_style("journal")` + `pf.Layout(mosaic, width=pf.SINGLE_COLUMN | pf.DOUBLE_COLUMN)`.
  `fig = lay.finalize()`, `print(pf.collisions(fig))` (must be `[]`), `pf.save(fig, OUT / name)` (svg + png).
- When paperfig lacks something generally useful (a plot type, a check, a legend helper), propose adding it to
  paperfig rather than working around it in the script.
- Python: `$(conda info --base)/envs/sccc/bin/python`.

## Palette (user, 2026-11-16)

- Categorical: the seaborn default palette, `sns.color_palette()` (matplotlib tab10), in its order, blue first. Never
  cycle past 10 colours: fold extra groups into "other" or facet.
- One colour per entity, the same in every figure. Methods:

  | entity | colour (tab10 index) |
  |---|---|
  | clocksync | blue (0) |
  | CHIRAL | orange (1) |
  | CoFE, clock genes | green (2) |
  | CoFE, genome-wide | red (3) |
  | CYCLOPS 2.0 | purple (4) |
  | TimeTeller (supervised reference) | grey (7), dashed |
  | timed oracle / recorded-time truth | black, dotted |

  | morning surgery / sampling window (e.g. 06-10 h) | cyan (9) |
  | afternoon surgery / sampling window (e.g. 14-18 h) | brown (5) |
  | peripheral tissues (vs neural) | olive (8) |
  | neural / ocular tissues | pink (6) |

  Other recurring entities (organ systems, cohorts, cell classes) get a fixed mapping the first time they appear;
  write it in this table.
- Neutral elements (identity lines, references, connectors, grids): grey `pf.GREY` or lighter. Text stays in ink
  (`#222222`) or grey, except a direct label next to a mark, which may take the mark's colour.
- Sequential: one hue, light to dark (`viridis` if a full map is needed). Diverging: `RdBu_r` centred on 0.
- Figures made before 2026-11-16 use Okabe-Ito (`pf.OKABE_ITO`): migrate them when they are next edited.

## Text

- **Units in square brackets:** `offset error [h]`, `amplitude error [|log2|]`, `age [years]`. Never parentheses.
- **Titles are the condition only:** `full day`, `hospital hours`, `control`, `knockout`. No statistics, sample sizes
  or sentences in titles: the numbers go into the slide bullets or the caption. A panel whose content is clear from its
  axes gets no title.
- **No processing jargon on axes:** drop words like "aligned", "rotated", "normalised", "posterior"; say it once in the
  caption or the speaker notes.
- Short, lower-case axis labels; a line break rather than a smaller font when a label is long.
- Legends: one row under the panels, or direct labels next to the lines when there are at most four series.

## Marks

- Show the data and its uncertainty, nothing decorative: no background density squares, shading or boxes behind
  summary points unless they carry a quantity the reader needs.
- Summary point + interval (median and IQR over draws, mean and circular SD per bin, estimate and bootstrap 95 % CI);
  say which in the caption.
- Bands at low alpha (≈ 0.1); the method of interest drawn thicker and on top.

## Statistics shown in figures

- Prefer robust, interpretable error measures in hours (median absolute error, with a donor bootstrap 95 % CI) over
  correlation coefficients.
- Do not compare Jammalamadaka-SenGupta circular r between methods when the phases span the whole day: its mean
  direction is ill-defined and the paired bootstrap interval of the difference is very wide (GTEx DIP, 2026-11-16:
  [−0.46, +0.40]).
- Every number on a figure or in its text is read from an output table written by the compute script.

## Data rules

- GTEx (time of death, protected phenotypes) is controlled access: plot aggregates only (bins or cells with ≥ 5
  donors), never one point per donor. Split compute (`*_compute.py`, writes aggregates to `output/`) from plotting
  (reads them).

## Procedure

1. Compute script writes aggregate CSVs (long runs through `paper/common/run.sh`); plot script reads only those.
2. Docstring of the plot script: what each panel shows, the data source, the status (registered / post hoc), output path.
3. Render, then **look at the PNG** (Read the image) for overlaps, clipped labels, wasted space, and check
   `pf.collisions(fig) == []`.
4. Commit scripts and the md text only (`figures/` and `output/` are gitignored); update the numbers in
   `paper/presentation/results.md` from the outputs.
