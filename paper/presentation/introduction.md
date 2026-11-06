# Introduction (talk section 1; source for `paper/manuscript/introduction.tex`)

Draft, 2026-11-06. One block per slide: **message** (one line), bullets with `\cite{}` keys from
`paper/manuscript/references.bib`, **visual**, and **notes** where useful. Every claim traces to a card in
`paper/manuscript/literature/notes/<key>.md`; numbers marked ✓ were re-read in the full text (2026-11-06).

Order of the argument (framing, 2026-11-05): the biology (organisation in health, desynchrony in disease), then the data
constraint (hospital hours), then clocksync as the tool that makes the question answerable.

---

## 1. The circadian cycle

### Slide 1.1: A clock in every cell

**Message:** Circadian rhythms are generated cell-autonomously, by a transcriptional feedback loop present in virtually all cells.

- Interlocked loops: CLOCK:BMAL1 drive PER/CRY, which repress them; REV-ERB/ROR and PAR bZip/NFIL3 loops give each
  target a characteristic phase \cite{takahashi2017transcriptional}.
- "Virtually all cells express these genes and have the capacity to generate circadian oscillations"
  \cite{takahashi2017transcriptional}; cultured fibroblasts oscillate after a serum shock \cite{balsalobre1998serum}.
- Single fibroblasts oscillate independently, with diverse periods and no coupling; population damping reflects loss of
  synchrony between cells, not damping of cells \cite{welsh2004bioluminescence, nagoshi2004circadian}.

**Visual:** schematic of the core loop (BMAL1/CLOCK → PER/CRY ⊣; REV-ERB ⊣ BMAL1), and a sketch of single-cell traces
drifting apart while the population average damps.

**Notes:** dephasing vs damping returns later: an averaged signal can lose amplitude without any oscillator losing it.
Same point at the cell scale in scRitmo \cite{salati2026inferring}.

### Slide 1.2: The body is a hierarchy of clocks

**Message:** The SCN keeps peripheral clocks in phase with each other; without it they keep ticking but drift apart.

- The SCN entrains peripheral clocks through autonomic, endocrine (glucocorticoid), temperature and feeding signals
  \cite{mohawk2012central, dibner2010mammalian}.
- Peripheral tissue explants oscillate for > 20 cycles in isolation, each with its own phase and period
  \cite{yoo2004periodluciferase}.
- After SCN lesion, rhythms persist but phase synchrony is lost "even from tissue to tissue within individual animals,
  which is a form of 'internal desynchronization'" ✓ \cite{yoo2004periodluciferase}; the SCN is required for synchrony
  between organs \cite{sinturel2021circadian}.

**Visual:** body schematic: SCN → arrows (nerves, hormones, temperature, feeding) → liver, heart, adipose, skin,
each a small oscillator dial; second panel with SCN removed and dials pointing in different directions.

**Notes:** this is the definition of synchrony we use: compartments of one individual agree in phase. `sinturel2021`
is Schibler/Dibner, not Naef.

### Slide 1.3: One clock, many outputs

**Message:** The core clock is shared and nearly in phase across organs, while its outputs are tissue-specific.

- In mouse, 43 % of protein-coding genes cycle in at least one of 12 organs ✓, but only 10 genes cycle in all of them,
  mostly core clock genes \cite{zhang2014circadian}. In baboon, up to 81.7 % of protein-coding genes cycle somewhere ✓
  \cite{mure2018diurnal}.
- Core clock genes peak within about 3 h of each other across mouse organs \cite{zhang2014circadian}, and within about a
  6 h window across baboon tissues, with conserved relative phases \cite{mure2018diurnal}.
- A shared clock module oscillates in synchrony across tissues with tissue-dependent amplitude (smallest in brain)
  \cite{yeung2018transcription}; the clock drives the high-amplitude rhythms, systemic and feeding signals the
  low-amplitude ones \cite{yeung2018transcription, weger2021systematic, atger2015circadian}.

**Visual:** schematic: one clock-gene phase wheel shared by all tissues (same gene order), with tissue-specific
amplitudes; output genes scattered differently per tissue.

**Notes:** this is the structural assumption clocksync uses (fixed gene acrophases, a compartment-specific amplitude
and offset), taken from CoPhaser \cite{paychere2025cophaser}.

### Slide 1.4: The human clock across tissues

**Message:** In humans, tissue clocks are synchronised within a few hours, but between-person estimates disagree, and
no one has measured it within one person.

- From untimed GTEx samples, "clocks across tissues were well synchronized, showing relative offsets of only a few
  hours" ✓; amplitude highest in metabolic tissues, lowest in brain and testis; damped with age; about twofold more
  rhythmic genes in women \cite{talamanca2023sexdimorphic}.
- Nearly half of protein-coding genes cycle in at least one of 13 human tissues ordered by CYCLOPS
  \cite{ruben2018database}.
- From timed data compared between people, brain clock gene expression is delayed by about 8.5 h relative to blood ✓
  \cite{hughey2016differential}.

**Visual:** two-panel cartoon: tissues as dials on one circle within a few hours (CHIRAL) vs brain and blood 8.5 h
apart (Hughey), labelled "between-person comparison".

**Notes:** the discrepancy motivates a within-person readout. CHIRAL assumes a common tissue shift across donors, so
it cannot tell whether a given person's organs agree. Our A-series answer comes in the results section.

---

## 2. Internal desynchrony

### Slide 2.1: What desynchrony means

**Message:** Desynchrony is a relation between clocks, distinct from misalignment with the environment or a broken clock.

- Mello et al. separate "environmentally induced desynchrony, misalignment between internal and external cues, or
  molecular disruption of the clock" \cite{mello2026rhythms}; Galinde et al. distinguish several levels, from cells
  within a tissue to tissues within the body \cite{galinde2023different}.
- "Although inter-tissue desynchrony is an established concept, experimental evidence is still sparse"
  \cite{galinde2023different}.
- Our usage: one compartment's clock relative to the clock state shared by the same person's other compartments.

**Visual:** three-row diagram: (a) external cycle vs body clock offset (misalignment); (b) organ dials out of phase with
each other (desynchrony); (c) one dial with no hand (disruption).

### Slide 2.2: Desynchrony can be induced, and it has costs

**Message:** In animals, conflicting cues pull organs apart at different speeds.

- Daytime feeding inverts peripheral clocks (8-12 h) while the SCN stays locked to light \cite{damiola2000restricted};
  liver shifts by 10 h within 2 days, lung more slowly \cite{stokkan2001entrainment}.
- Under conflicting light and feeding cycles, tissues lag behaviour by different amounts (after 4 days, eWAT 1.6 h vs
  7.1 h depending on the feeding schedule) and phase coherence across tissues is lost \cite{heyde2022induction}.
- Repeated light-schedule advances kill aged mice (survival 47 % vs 83 % in unshifted controls at day 56); the authors
  speculate on internal desynchrony \cite{davidson2006chronic}.

**Visual:** time course schematic: SCN dial fixed, liver dial rotating fast, lung slowly, after a feeding inversion.

**Notes:** the mechanism behind "internal desynchrony": different tissues weight different cues and reset at
different rates.

### Slide 2.3: In humans, desynchrony is measured one compartment at a time

**Message:** Human misalignment studies show tissue clocks that fail to follow behaviour, but only one tissue at a time,
under laboratory protocols.

- Mistimed sleep cuts rhythmic blood transcripts from 6.4 % to 1.0 % \cite{archer2014mistimed}; a simulated night shift
  dampens most blood rhythms rather than shifting them \cite{kervezee2018simulated}.
- After an inverted behavioural cycle, the muscle clock keeps its original pattern ("had not adapted to the behavioral
  rhythm"), with lower insulin sensitivity \cite{wefers2018circadian}; misalignment raises glucose and lowers leptin
  \cite{scheer2009adverse}.
- Long-term rotating night work is associated with breast cancer (≥ 30 y: RR 1.36) and heart disease (≥ 10 y: HR 1.18)
  \cite{schernhammer2001rotating, vetter2016association}, though large studies find no increased cancer risk
  \cite{sancar2021clocks}.

**Visual:** small table: study, tissue sampled, protocol, outcome; one column "two compartments of the same person?"
with "no" in every row.

**Notes:** `kervezee2019individual` and `skene2018separation` are plasma or between-subject: not two compartments in
the same people (synthesis §3).

### Slide 2.4: Disease without an external cause

**Message:** Disease may desynchronise the body's clocks, but the human evidence is cross-region, post-mortem and
indirect.

- Major depression: consistent rhythms across six brain regions in controls, but only 2 of the top 16 genes significant
  in more than one region in patients \cite{li2013circadian}.
- Ageing changes rhythms in human prefrontal cortex (loss, gain and phase shifts) \cite{chen2016effects}.
- In Alzheimer's disease, core clock rhythms are nearly unchanged across cell types, while outputs are disrupted
  \cite{hollis2025reconstructed}.

**Visual:** none, or a one-row summary of the open question: "Are a patient's organs out of phase with each other?"

**Notes:** `li2013` and `chen2016` are timed post-mortem cohorts (time of death); `hollis2025` assigns one donor phase
to all cell types, so cell-type offsets are zero by construction.

---

## 3. Interaction with cancer

### Slide 3.1: Host disruption promotes tumours in mice

**Message:** Disrupting the host's clock or its synchrony accelerates tumour growth in mice.

- SCN lesion: tumours grow two to three times faster and survival is shorter; the authors propose loss of temporal
  organisation across tissues as a cause \cite{filipski2002host}.
- Chronic jet lag increases lung tumour burden and shortens survival in a Kras model
  \cite{papagiannakopoulos2016circadian}, and causes liver cancer in wild-type mice (8.75 % vs 0 %)
  \cite{kettner2016circadian}.
- Clock-gene loss can act inside the tumour too (Per2, Bmal1 in lung) \cite{papagiannakopoulos2016circadian,
  fu2002circadian}, but clock genes are not general tumour suppressors \cite{sancar2021clocks}.

**Visual:** schematic mouse with SCN lesion or shifting light-dark bars, arrow to a larger tumour; citation under each.

### Slide 3.2: The tumour talks back

**Message:** A tumour can rewire the host's peripheral rhythms without touching the host's core clock.

- Lung adenocarcinoma in mice rewires hepatic rhythms: 46 % of shared rhythmic liver genes change phase by ≥ 1 h,
  while core clock genes in the liver are unchanged \cite{masri2016lung}.
- The effect is organ-specific (WAT and muscle differ) and is linked to systemic inflammation (IL-6)
  \cite{masri2016lung}.

**Visual:** schematic: tumour → serum cytokines → liver; liver core-clock dial unchanged, output genes rotated.

**Notes:** this is why we read the core clock only: a host-referenced core-clock readout separates clock coupling from
output rewiring.

### Slide 3.3: Tumour clocks range from absent to intact

**Message:** Tumour clocks have been reported as arrhythmic, shifted or in phase; the cell cycle is one expected source of offset.

- The literature reports tumour rhythms that are arrhythmic, rhythmic but out of phase, or in phase with normal tissue
  \cite{sancar2021clocks}.
- Clock-gene co-expression is weaker in tumours than in normal tissue across cancers \cite{ye2018genomic,
  shilts2018evidence}, but strong-clock luminal A breast tumours keep a phase matched to their normal tissue
  \cite{li2024tumor}; glioblastoma stem cells keep an intact clock and depend on it \cite{dong2019targeting}.
- The cell cycle couples to the clock, predominantly cell cycle → clock, and accelerates it \cite{bieler2014robust,
  droin2019low}; across mouse tissues, phase offsets (spanning 1.5 h ✓) correlate with Mki67 and Myc
  \cite{droin2019low}.

**Visual:** three dials (tumour vs host): no hand / rotated hand / aligned hand, with the clocksync readout under each
(|z| ≈ 0; low c with an offset; high c). Mapping proposed in synthesis §4: needs your OK (see open questions).

**Notes:** `feillet2014phase` (Delaunay, Rand, van der Horst; not Naef) shows multiple locking modes (1:1, ~3:2, ~5:4).

### Slide 3.4: Timing treatment: motivation, not proof

**Message:** Treatment timing changes toxicity and sometimes response, but survival benefit is inconsistent and depends on sex.

- Chronomodulated 5-FU, folinic acid and oxaliplatin: response 51 % vs 29 %, severe mucositis 14 % vs 76 %; 3-year
  survival 22 % vs 21 % ✓ \cite{levi1997randomised}.
- Of 18 randomised chronochemotherapy trials (2,547 patients), 3 improved efficacy and 11 reduced toxicity ✓; in
  colorectal cancer, survival was prolonged in men and reduced in women \cite{levi2024circadian, levi2010circadian}.
- Earlier-in-the-day immune-checkpoint infusions associate with longer survival (13 studies, 1,663 patients, OS HR 0.50)
  \cite{landre2024effect, qian2021effect}.
- Between people, rest-activity and skin-temperature acrophases differ by up to 12 h \cite{levi2024circadian}; "a major challenge … is the development of
  noninvasive methods to assess the circadian parameters of cancers" \cite{sancar2021clocks}.

**Visual:** forest-style summary of the three published aggregates (no new analysis), or text only.

**Notes:** all ICI data are retrospective. Do not imply chronotherapy works; the point is that a patient's and a tumour's
clock state are unknown at treatment time.

---

## 4. Technical challenges

### Slide 4.1: How human clock phase is measured

**Message:** Every existing approach assigns a time to each sample, from timed training data or from a population that covers the cycle.

- **Supervised predictors**, trained on timed series (mostly blood): ZeitZeiger, TimeSignature, BodyTime, PLSR, TimeTeller
  \cite{hughey2016zeitzeiger, braun2018universal, wittenbrink2018high, laing2017blood, vlachou2024timeteller}.
- **Orderers**, which place untimed samples on a cycle: CYCLOPS and CYCLOPS 2.0, CHIRAL, CoFE, DCPR
  \cite{anafi2017cyclops, li2024tumor, talamanca2023sexdimorphic, ananthasubramaniam2025rhythm, han2025dcpr}.
- **Single-cell methods**, which order cells within a sample: Tempo, CoPhaser, scRitmo \cite{auerbach2022tempo,
  paychere2025cophaser, salati2026inferring}.

**Visual:** methods table condensed to three columns (unit, supervision, needs 24 h coverage), from synthesis §2.

**Notes:** CHIRAL, CoPhaser and scRitmo are the lab's; clocksync takes from CoPhaser the fixed gene acrophases with
fitted amplitude and mean shifts.

### Slide 4.2: Orderers need the whole day

**Message:** By their authors' own account, orderers need samples spread over the 24 h cycle.

- "It requires data from the entire periodic cycle to form an ellipse. Biopsies are almost exclusively obtained during
  the day." ✓; "> 250 samples are required to order biopsy samples" ✓ \cite{anafi2017cyclops}.
- "Most computational methods for quantifying circadian rhythmicity and inferring the status of the clock require that
  samples be acquired over the entire 24-h cycle." ✓ \cite{shilts2018evidence}.
- "clinical samples are collected during surgical working hours that are not evenly spread across the day" ✓
  \cite{grabe2024quantification}; orderers degrade on incomplete cycles \cite{han2025dcpr}.

**Visual:** the three quotes, large type, one per line, with author and year.

### Slide 4.3: Clinical tissue is collected in hospital hours

**Message:** Surgical and biopsy cohorts cover a working-day window, while GTEx deaths cover the whole day.

- Biopsies are taken "during surgical working hours (∼6:00 AM to 6:00 PM)" ✓ \cite{anafi2017cyclops}.
- Recorded surgery times of 302 glioblastoma patients analysed from 10:00 to 18:00 ✓ \cite{gonzalezaponte2025circadian};
  phases inferred for TCGA biopsies "were concentrated during the presumed working day" ✓ \cite{li2024tumor}.
- GTEx deaths span day and night: 351 donors between 08:00 and 17:00, 315 between 21:00 and 05:00, 266 in between ✓
  \cite{wucher2023day}. GTEx is a non-cancer post-mortem reference, not a healthy cohort \cite{lonsdale2013genotype}.

**Visual:** 24 h clock-face (or linear 0-24 h) bar: shaded hospital window (10-18 h, glioblastoma; 6-18 h stated)
vs GTEx day / night / twilight counts as three bands with donors per hour (published aggregates only).

**Notes:** GTEx is therefore where we can simulate hospital hours with known truth (A123). Only published aggregates
appear on this slide; no per-donor times.

### Slide 4.4: Hospital-hours data make "uncoupling" ambiguous

**Message:** Published tumour "uncoupling" comes from methods that need full coverage, applied to working-hours cohorts.

- CoFE places TCGA tumour times throughout the cycle in all cancers, although surgery is daytime, and finds
  tumour-normal differences "distributed throughout the circadian cycle" ✓ → "uncoupled"
  \cite{ananthasubramaniam2025rhythm}.
- TimeTeller places 93.2 % of breast tumours in the working day ✓, partly after reassigning off-window samples, and
  reads r ≈ 0 with sampling time (n = 108, narrow window) as uncoupling \cite{vasilyev2026robust}.
- CYCLOPS 2.0: tumour and paired-normal phases agree in strong-clock luminal A tumours ✓, poorly overall
  \cite{li2024tumor}; noise in weak clocks predicts the same.

**Visual:** schematic: a 10-18 h band of samples forced onto a full circle → random tumour-normal differences.

**Notes:** this is a reason for the method, not a result. Our own tests (A145, A204) belong in the results section.

### Slide 4.5: What is needed

**Message:** A within-person readout that compares each compartment with the person's own clock, without assigning a time to any sample.

- Closest precedents: CCD is time-free and robust to 8 h windows but group-level and "cannot detect a phase difference"
  \cite{shilts2018evidence}; Grabe et al. order each tumour by its paired normal but assume full coverage
  \cite{grabe2024quantification}.
- CHIRAL "assumes that tissues obtained from the same individual at the same time are at the same molecular phase"
  ✓ \cite{li2024tumor, talamanca2023sexdimorphic}; clocksync estimates that relation instead.
- Requirements: within-person compartments; no sample time; per-compartment coupling c, relative amplitude |z| and
  offset arg z; a count noise model; a check that the shared signal is clock-shaped.

**Visual:** transition to the model animation (`figures/model_animation/model_01.mp4`): donor clock state shared by
compartments.

**Notes:** constant offsets (sampling hour, chronotype) are common to all compartments of a donor and cancel; this is
why hospital-hours samples are usable.

---

## Proposed figures (to build after the text is agreed)

| slide | figure | type | data |
|---|---|---|---|
| 1.1 | core loop and dephasing traces | schematic | none |
| 1.2 | SCN → organs, intact vs lesion | schematic | none |
| 1.3 | shared clock wheel, tissue amplitudes | schematic | none |
| 1.4 | within-person vs between-person dials | schematic | numbers from cards |
| 2.1 | misalignment / desynchrony / disruption | schematic | none |
| 2.3 | human studies table | table | cards |
| 3.3 | three tumour scenarios → c, \|z\|, arg z | schematic | none |
| 3.4 | three published aggregates | text or forest | levi1997, levi2024, landre2024 |
| 4.3 | 24 h axis: hospital window vs GTEx day/night/twilight | data plot (published aggregates) | wucher2023day, gonzalezaponte2025, anafi2017 |
| 4.4 | partial window forced onto a circle | simulated schematic (`clocksync.simulate`) | script output |

## Open questions for you

1. **Slide 3.3 mapping** of the three tumour scenarios onto |z| ≈ 0 / low c with offset / high c: it puts a claim about
   our readouts in the intro. Keep, or move it to the results?
2. **Slide 2.3 epidemiology:** keep the shift-work associations (with Sancar's null large studies), or drop them? IARC
   cards (`iarc2020night`, `straif2007carcinogenicity`) are in the bib but not carded.
3. **Slide 4.4:** name CoFE and TimeTeller on the slide, or describe them as "orderer-based" and "supervised" only?
4. **Length:** 15 slides. Merge 1.1 + 1.2 and 2.1 + 2.2 for a short talk?
