# kettner2016circadian
**Ref:** Kettner NM, Voicu H, Finegold MJ, Coarfa C, Sreekumar A, Putluri N, Katchy CA, Lee C, Moore DD, Fu L, Cancer Cell 2016;30:909-924. DOI 10.1016/j.ccell.2016.10.007. Type: biology (mouse, metabolism and carcinogenesis). Lab: Fu / Moore. Data: GEO GSE75475.
**Source:** full text (main text; STAR-methods jet-lag schedule not in the extract)

## Question and design
C57BL/6J mice, lifelong (4-90 weeks) chronic jet lag (CJ) versus steady LD (Ctrl), in wild type and in Per1/Per2 double knockout, Cry1/Cry2 double knockout, liver-specific Bmal1 (Alb-cre;Bmal1fl/fl), Car-/- and Fxr-/- lines. Host (whole-body, with liver clock) perturbed; spontaneous liver tumours, so the tumour clock is not measured and not compared with host. **Sampling:** timed but sparse: liver microarray at ZT2, 10, 18 in 12- and 30-week mice; metabolomics rhythmic at 12 and 30 weeks; western and RT-PCR at ZT2 reference; HCC incidence by 90 weeks.

## Method (if any)
Kaplan-Meier survival, histology, microarray (hierarchical clustering, overlap with human HCC signatures, hypergeometric test), targeted metabolomics, ChIP, reporter assays. No circadian inference method.

## Key results
- Cancer outcome: Table 1 shows HCC 0 % in WT Ctrl (n = 110) vs 8.75 % in WT CJ (n = 80), first found at 78 weeks; HCC is male-dominant. "unirradiated jet-lagged WT and mutant mice displayed a similar risk of lymphoma and HCC" (Results).
- Clock disruption in the host liver: "these results demonstrated that chronic jet lag is sufficient to disrupt the liver clock independent of circadian gene mutations"; Bmal1, Clock, Per1, Per2, Cry1 and Nr1d1 were deregulated by RT-PCR (Fig 4C).
- In liver-specific Bmal1 knockout, "Per2 mRNA still displayed a shifted and dampened circadian expression profile" while Cry1 was constitutively high.
- Chronic jet lag "abolished BMAL1 binding to both Car and Per1 promoters" (ChIP) and made Car arrhythmic and elevated in WT liver.
- Genome-wide: "we carried out a large-scale microarray analysis using total liver RNA prepared from 12 and 30 week old control and jet-lagged WT mice at ZT2, 10, and 18"; deregulated signature overlaps human HCC; nuclear-receptor, cholesterol, bile-acid and xenobiotic pathways are top.
- Sympathetic route: "chronic jet lag significantly increased SNS tone in the sleep phase" (AP1/CREB activation of Car).
- Genetics: Fxr-/- under CJ has more than 2-fold higher HCC incidence; Car-/- is completely resistant to spontaneous HCC (Table 1, Fig 5), showing the oncogenic effect runs through CAR.
- Authors' framing: "In the hierarchical circadian system, peripheral organs rely on signaling from circadian output pathways to maintain synchrony with the central clock."

## Not addressed / limits
- Only the liver clock (and a few other tissues) is measured; no other organ and no tumour clock, hence no host-tumour comparison or inter-organ synchrony metric.
- Three time points per condition cannot give acrophase or amplitude per gene; "clock disruption" is read from expression differences at ZT2/10/18.
- Incidence low in WT CJ (8.75 %, n = 80); disruption is confounded with metabolic syndrome (obesity, leptin resistance, cholestasis).
- Jet-lag protocol and per-timepoint sample sizes are not in the extracted text.

## Link to us
- **Intro claim(s):** the refs_*.md row is supported (chronic jet lag alone in WT mice produces HCC via NAFLD and FXR/CAR). It supports systemic misalignment being carcinogenic, not any claim about tumour clocks.
- **Relation to clocksync:** data we could use (GSE75475, liver ZT2/10/18, Ctrl vs CJ; planned as scout A125 in paper/docs/provisory_plan_paper.md, not registered or scored). Related done analysis: A201 (aorta, 6 h advance, Auerbach 2026).

## Ideas
- Register A125 on GSE75475: in the same liver, does the clock-gene panel show lower amplitude |z| and offset under CJ at ZT2/10/18 (three known times give a truth check), and does the gate remain valid?
- Missing reference: Lee 2010 (Cancer Cell, Per mutants, jet lag, SNS, p53 and c-Myc), Kettner 2015 (leptin resistance under circadian disruption).
