# druzd2017lymphocyte
**Ref:** Druzd, Matveeva, Ince, Harrison, He, Schmal, Herzel, Tsang, Kawakami, Leliavski, Uhl, Yao, Sander, Chen, Kraus, de Juan, Hergenhan, Ehlers, Koletzko, Haas, Solbach, Oster, Scheiermann, Immunity 2017 (46:120). DOI 10.1016/j.immuni.2016.12.011. Type: biology (mouse immunology). Lab: Scheiermann (LMU Munich).
**Source:** full text

## Question and design
Do lymphocytes traffic through lymph nodes (LN) in a circadian way, and do lymphocyte-intrinsic clocks control it? Mouse (C57BL/6), 12:12 light-dark (ZT), constant darkness (DD, CT), reversed light (DL). Inguinal and other LN, blood, efferent lymph (cannulated), adoptive transfers at ZT5 or ZT13, homing block (anti-integrin) and egress block (FTY720), sorted CD4+ T cells for clock gene qPCR and Ccr7 / S1pr1, lineage-specific knockouts: Bmal1flox x Cd4-cre (T cells), Bmal1flox x Cd19-cre (B cells), S1pr1+/flox x Cd4-cre, Ccr7-/-; EAE immunization at ZT8 vs ZT20, H. pylori and influenza. n typically 3 to 34 mice per point. **Sampling:** timed, full 24 h in ZT or CT.

## Method (if any)
Counts and flow cytometry per time point, qPCR (whole LN and sorted cells), promoter-luciferase assay (S1pr1 in HEK293 with Clock and Bmal1), simple mathematical homing/egress model.

## Key results
- LN cellularity oscillates and is circadian: "Lymphocyte homing to lymph nodes peaked at night onset, with cells leaving the tissue during the day." Blood lymphocytes peak near ZT5, LN counts of CD4, CD8 and B cells about 8 h later, at ZT13; rhythm persists in DD and inverts under reversed light.
- Both homing and egress are rhythmic; efferent lymph peaks at ZT9, low at ZT21; LN half-life of CD4 cells 12 h (injected ZT13) vs 9 h (ZT5); B cells 16 vs 13.5 h. The model fits best with both rhythmic.
- Lymphocyte-intrinsic clock: "Together, these data provide evidence for a functional role of cell-autonomous clocks in lymphocyte migration." T-cell or B-cell Bmal1 deletion ablates LN count rhythms, rhythmic homing, rhythmic CCR7 (surface and mRNA) in CD4+ T cells, and rhythmic lymph counts; BMAL1-deficient CD4+ T cells show no time-of-day difference in LN half-life.
- Clock-controlled targets: Ccr7 (high at night onset) and S1pr1 (high by day, ablated in Bmal1-deficient CD4+ T cells). BMAL1 and CLOCK reduced S1pr1 promoter luciferase in HEK293 cells, so "expression of S1pr1 is regulated by BMAL1 and CLOCK".
- Not purely cell-autonomous: adoptive transfers of day and night cells into day and night recipients show "a mixed contribution of both lymphocyte and microenvironment timing"; CCL21 on high endothelial venules oscillates (non cell-autonomous).
- Clock genes: "LNs exhibit oscillations of clock genes" (Fig. 3A, whole LN) and sorted CD4+ T cell clock gene profiles in controls versus Bmal1 knockout (Fig. 3B); the figure panels are labelled Bmal1, Clock, Per1, Per2, Cry1, Nr1d1, Dbp (which genes oscillate, and amplitudes, are not stated in the text).
- Function: EAE severity depends on immunization time (ZT8 worse than ZT20) and this is lost in T cell Bmal1 knockout; also H. pylori (LN counts highest at ZT7) and influenza (more CD8+IFN-g+ T cells at ZT8 than ZT20).
- Other organs: thymus (data not shown) and bone marrow show no overt oscillation of absolute numbers.

## Not addressed / limits
- Mouse only. Clock-gene rhythms are shown in whole LN and sorted CD4+ T cells; no single-cell data, no tissue-resident lymphocytes in non-lymphoid organs, no human data.
- Cell-autonomy is inferred from knockout phenotypes and trafficking readouts; free-running clock gene oscillations in isolated lymphocytes (e.g. ex vivo luciferase reporter) are not shown. Bmal1 deletion removes the clock but also non-clock BMAL1 functions.
- The microenvironmental (CCL21), autonomic or humoral drivers are outside the scope of the paper. A previous lymphocyte Bmal1 knockout (Hemmers and Rudensky) "yielded no obvious phenotype", explained by testing at one time of day.

## Link to us
- **Intro claim(s):** immune cells carry clocks that gate trafficking and adaptive responses (counterpoint to the idea that lymphocytes have no readable clock of their own).
- **Relation to clocksync:** counterpoint / background. Evidence is for clock-gene-dependent expression of Ccr7 and S1pr1 in sorted lymphoid T cells (blood and LN of mice, fully timed). It does not show a self-sustained cycle in tissue-resident lymphocytes; mixed cell-autonomous and microenvironmental timing is also consistent with a driven (below Hopf) state in which BMAL1/CLOCK targets follow a systemic signal, but a Bmal1 phenotype shows the clock machinery is functionally used.

## Ideas
- In human sc atlases, test whether Ccr7 / S1PR1 / CD69 and other Bmal1/Clock targets track the tissue's clock state in T cells (paired with PER1, NR1D1 as hormone-induced genes), and compare LN or blood T cells with tissue-resident ones.
- Check which of Druzd's panel genes (Bmal1, Clock, Per1, Per2, Cry1, Nr1d1, Dbp) are detectable in T cells in our atlases; the weak core-loop signal in some atlases could be this.
- Cited but not in export: Hemmers and Rudensky 2015 (lymphocyte Bmal1 null, no phenotype), Fortier 2011 (Clock mutant T cells).
