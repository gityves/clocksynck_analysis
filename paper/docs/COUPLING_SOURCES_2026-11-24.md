# Where an organ's coupling to the body's time comes from: literature, our evidence, what we add (2026-11-24)

On the user's request (2026-11-24): extend the health part with all the evidence we have on the source of an organ's
coupling to the body's time (human and mouse), and state what we add on top of the literature. Builds on
`MODEL_federation_2026-11-19.md` (parts 3-6), `PAPER_PLAN_2026-11-19.md` sections 3 and 11.1, and slides 3.3-3.4b of
`presentation/results_outline_2026-11-23.md`. Literature from the cards in `manuscript/literature/notes/`; numbers
from `results_preliminary.md` entries cited. No new analysis in this document.

## 1. The question, split in two

An organ's state is x_{p,o} = G_o exp(i phi_o) theta_p + delta_{p,o} (person p, organ o). Coupling to the body's time,
c_o = S_o / (S_o + D_o) with S_o = G_o^2 V, has two sources:
- **G_o, the transfer** of the body's time to the organ (its gain): how strongly the organ expresses the shared time.
- **D_o, the departure** added at the organ, after the transfer: nearly the same size in every organ (A294b, slope 0.07),
  so coupling ranks organs by G_o (slides 3.4a-b).

So "where does coupling come from" is "what sets G_o" (and what sets phi_o, the organ's place in the body's time).
The literature answers phi_o well in mouse; it says little about G_o.

## 2. What the literature establishes (mouse unless stated)

| topic | consensus | cards |
|---|---|---|
| peripheral clocks run on their own | explants and cells oscillate for days to > 20 cycles; tissue-specific periods; phase not explained by period, so the inputs or the responsiveness differ by tissue | balsalobre1998serum, yoo2004periodluciferase, dibner2010mammalian |
| drives set peripheral phase | feeding inverts the liver (8-12 h) without moving the SCN; glucocorticoids shift peripheral clocks but are not needed for the steady-state phase (they brake feeding-driven shifts); body temperature cycles entrain the periphery (type 0 resetting), the SCN resists | damiola2000restricted, stokkan2001entrainment, balsalobre2000resetting, leminh2001glucocorticoid, brown2002rhythms, buhr2010temperature |
| organs weigh the drives differently | under daytime feeding the liver shifts about 12 h, kidney and heart 4-8 h, lung about 0; amplitude falls in every tissue ("a conflicting signal paradigm", effects "beyond the additive") | stokkan2001entrainment, manella2021liverclock, izumo2014differential |
| the SCN holds the organs together, not the cells within an organ | without the SCN, organs desynchronise from each other ("internal desynchronization") and the ensemble amplitude falls; hepatocytes stay coherent within the liver; light-dark keeps the organ map without an SCN clock; in constant dark organs drift and damp over days | yoo2004periodluciferase, sinturel2021circadian, finger2021peripheral, husse2014lightdark, izumo2014differential, tahara2012in |
| the SCN opposes conflicting drives | the liver re-phases to feeding in 1-2 days without the SCN, 6-7 h in 7 days with it | saini2013realtime, leminh2001glucocorticoid |
| jet lag | the SCN shifts at once, organs over days at organ- and gene-specific speeds; liver slowest after a delay; transient internal misalignment | yamazaki2000resetting, kiessling2010adrenal, husse2014lightdark |
| local clocks need the system for full rhythms | a clock in only the liver or epidermis keeps 10-15 % of rhythmic genes, phase advanced 4 h (liver), lost in constant dark; feeding restores part | koronowski2019defining, welz2019bmaldriven, greco2021integration, kornmann2007system |
| human | metabolic tissues carry the largest clock amplitudes in GTEx, brain and testis the lowest; organs within a few hours of each other (CHIRAL); a 12 h behavioural shift blunts the muscle clock without a measured shift in 3.5 days | talamanca2023sexdimorphic, wefers2018circadian |

What the literature does not give: (i) any measure of how strongly each organ is coupled to the body's time in people;
(ii) what sets an organ's clock amplitude (G_o) beyond "metabolic tissues are loud"; (iii) whether removing a drive and
conflicting drives act alike on amplitude, quantitatively; (iv) whether the human organ arrangement depends on the
brain, or on the external day.

## 3. What we add

| # | claim | our evidence | beyond the literature | status |
|---|---|---|---|---|
| 1 | Human organs share one body time with fixed places in it (the organ map) | A123 (CHIRAL concordance); reference map (slide 3.1); A293 Part 1: order kept in slow deaths (rho 0.784) and ventilated donors with a living brain (0.690) although the body clock no longer reads the hour of death (R 0.16-0.25 against 0.63) | the arrangement is internal, held without the external day; first human test | registered re-tests (not blind) |
| 2 | The central drive sets the arrangement, not the amplitude | brain death: order changed (D_rho 0.621; 0.584 adjusted, A317 / A318), shared clock state kept (A317 H317.2 FAIL), amplitude equal to equally ventilated living-brain donors (A315, clock-specific ratio 1.05; A335) | the human counterpart of SCN lesion, with the opposite amplitude result on the short term (mouse: damping over days to weeks) | aliased with procurement; not blind |
| 3 | Intensive care costs amplitude, not order | fast over ventilated: clock-specific log variance ratio +0.404 (A315 / A335); falls with ventilation time (descriptive) | amplitude depends on the systemic state (behaviour, feeding, sleep), not on the brain per se (with 2) | registered |
| 4 | Removing a drive keeps amplitude; conflicting drives cost it | mouse liver without a feeding rhythm keeps 0.724 of its amplitude, far above what additive drives predict (A297); conflict costs amplitude in every tissue, most where the weights are balanced (WAT 0.32, VAT 0.11), least where one drive dominates (lung, heart, kidney 0.50-0.65; liver 0.44-0.71); a pure mix of drives is rejected in two datasets (A322, A324) | quantifies Manella's "beyond additive": self-sustained cells re-synchronise at an intermediate phase (lung keeps about half where a mix leaves a quarter); one weight b per organ | A297 not blind; A322 not blind, A324 blind for heart, kidney, VAT |
| 5 | In people, a behavioural shift moves the clock | night shift, day 3: muscle core loop +4.06 h, excursion 0.49 (A257, A335) | Wefers 2018 saw blunting without a measured shift; we measure the shift and the amplitude cost | core loop blind (A335) |
| 6 | Organs differ 9.2-fold in gain (transfer), not by their cells and not by post-mortem attenuation | class gains conserved across organs (A293 H293.4; core loop in four of five classes, A350); composition predicts 1.29-fold against 2.89-fold observed (A323); gain does not move with ischaemic time or RIN (A356, 1.01 / 0.99) | Talamanca 2023 ranked amplitudes; we show the ranking is an organ property beyond cell content | A323, A356 registered; living reproduction not shown (A355) |
| 7 | Coupling to the body ranks organs by gain, because the departure has a fixed size | S spans 85-fold, D 4.6-fold, slope 0.07 (A294b); implied phase scatter 3.0 h (skin) to 1.1 h (oesophagus muscularis), matching A259 (Spearman 0.945); departure clock-specific (2.6 times matched), a coherent clock-state shift (A329 erratum), as much amplitude as phase (A311) | synchrony and coupling metrics used in the field read loudness (A291); louder organs stay closer to the body's time in hours | A294b registered; slide 3.4b displays post hoc |
| 8 | Each person has a volume that scales all organs | person amplitude a strong person-level property (half-host radii 0.75-0.86); every organ follows it nearly in proportion (A312, elasticity median 1.19) | a systemic amplitude factor in people; its source is open (could include donor-level technical scale) | registered |
| 9 | A person-level glucocorticoid-like state is read unequally by organs | it explains a small, clock-specific part of each organ's departure (A357 H357.2, R^2 0.027 against 0.008; muscle 0.195) and tracks the glucocorticoid-target state (r 0.77, descriptive) | the human counterpart of glucocorticoid resetting (Balsalobre 2000), as a source of person-specific internal desynchrony; small | registered, blind |
| 10 | Within an organ, cell classes share its time; composition moves gain, not phase | classes within about 1 h (A273, A279, A293); composition 0.13 h per SD (A295) | human cell-class resolution of Sinturel's "synchrony within an organ" | registered |

## 4. The missing piece: what sets an organ's gain G_o

Candidates, each with what would test it:

| candidate | prediction | test | data |
|---|---|---|---|
| (a) drive coherence: an organ is loud when the drives it weighs agree in phase in daily life | gain lower where the weights are balanced across drives with different phases; the mouse weight b (A322, A324) predicts the human gain ranking of matched organs | descriptive cross-species (about 6 organs: liver, white / visceral fat, lung, heart, kidney, muscle); too few organs for a strong test | have |
| (b) drive reception: organs with more receptors for the systemic signals are louder | gain rises with the organ's responsiveness to the systemic signals (glucocorticoid-target response, A357 per-organ R^2; NR3C1 expression), beyond matched genes | GTEx, per organ; risk of circularity (glucocorticoid-target genes include clock outputs) | have |
| (c) cell-level: louder or better-synchronised cells | per-cell amplitude or synchrony differs by organ for the same class | needs deep single-cell data; public 10x cannot decide (A334) | missing |
| (d) drive dynamics after jet lag: weights set transient amplitude and path | organs with balanced weights (fat) lose most amplitude transiently; feeding-dominant organs (liver) follow feeding and activity, SCN-dominant organs (lung, heart, kidney) follow the SCN; the arrangement is disrupted transiently and restored (unlike brain death) | out-of-sample: b from A322 / A324 predicts the jet-lag amplitude dip, b(1 - b) at a given conflict angle | in-house jet-lag RNA-seq (PAPER_PLAN 6.4; details needed from the user) |
| (e) person volume (part 8) | the person's amplitude follows systemic state: intensive care, sleep, activity, glucocorticoid state | GTEx: person amplitude against the A357 person factor, Hardy scale, ventilation (aggregates); age and sex already equal (A208, A231) | have |

Readings that matter for the paper: (a) or (d) would make the organ gain a property of the drive configuration (the
federation's "loudness follows how strongly and coherently the drives agree", MODEL part 5 and the model rules),
linking the healthy organ map to the conflict results; (b) would make it a property of reception; (c) stays open.

## 5. Proposed changes to the paper plan (for the user)

1. Section 2.1 gets a sub-block "Where coupling to the body's time comes from" with the table of section 3 (claims
   1-10) in a shorter form, and the decomposition of section 1 (coupling = transfer and fixed departure).
2. Section 4 (what is new) adds: the human organ arrangement is internal and set by the brain, but amplitude is not
   (claims 1-3); removal against conflict, quantified (claim 4); coupling ranks organs by transfer (claim 7); the
   person volume (claim 8).
3. Figure 4 or 5: a panel joining the conflict results (mouse weights and amplitude kept) with the human gain ranking,
   only if test (a) or (d) is run; otherwise the conflict panel stays in Figure 5 as now.
4. Registrations to consider: (d) jet lag (needs the data), (e) person volume (cheap), (b) reception (cheap, with a
   circularity check); (a) descriptive only.
5. Literature to card: Wehrens 2017 (human meal timing shifts adipose PER2 about 1 h), cited by the human part of claim
   5 if used.

## 6. Outcomes (2026-11-24)

- (a) A358, descriptive: not seen (rho -0.6 between mouse balance b(1 - b) and GTEx gain, 5 organs; criterion -0.7).
- (b) A359: H359.1 FAIL (GC reception slope against gain rho 0.20, p 0.16); the IEG reception slope tracks gain
  (rho 0.44, p 0.016; companion, a lead only).
- (e) A360: H360.1 PASS (person volume falls with the person GC-target state, coefficient -0.235, partial R^2 0.061,
  p 0.0006); H360.2 FAIL (no technical scale from RIN or ischaemic time).
- What sets the organ gain stays open; the person volume has a first systemic correlate.
