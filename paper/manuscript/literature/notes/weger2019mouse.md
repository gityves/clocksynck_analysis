# weger2019mouse
**Ref:** Weger, Gobet, Yeung, Martin, Jimenez, Betrisey, Foata, Berger, Balvay, Foussier, Charpagne, Boizet-Bonhoure, Chou, Naef, Gachon, Cell Metab 29:362 (2019). DOI 10.1016/j.cmet.2018.09.023. Type: biology. Lab: Gachon (Nestle, EPFL, Queensland) with Naef (EPFL).
**Source:** full text (main text with STAR Methods)

## Question and design
Does the microbiome shape daily rhythms of host gene expression and metabolism, and sex differences in them? C57BL/6J germ-free (GF) vs conventionally raised (ConvR) male and
female mice, 12 h light : 12 h dark, 9 to 16 weeks old; RNA-seq of liver, duodenum, ileum, perigonadal white adipose tissue (WAT); liver and serum metabolomics; GH secretion kinetics; ghrelin and GH injections;
fecal microbiota profiling and co-housing. "All boxplots represent at least 12 independent biological replicates per condition" (physiology panels).
**Sampling:** timed, time series across the day (zeitgeber time; the exact interval and number of time points are in Figure S1A and are not in the extracted text).

## Method (if any)
Rhythmicity by harmonic regression with BIC-based model selection across 2 conditions (GF, ConvR) or 4 (GF/ConvR x sex), as in Atger 2015: models 1 to 5 (e.g. model 4 rhythmic in both, model 2 males only, model 3 females only, model 5 different phase or amplitude);
genes with BIC weight below 0.4 labelled ambiguous. Phase change below 4 h or amplitude fold below 2 called mild.

## Key results
- The clock is nearly unaffected, rhythmic physiology is: "genes showed only subtle changes in phase and amplitude of mRNA abundance across" tissues; "key metabolic processes exhibited an altered rhythmic expression in GF mice."
- Liver of ConvR mice: 71 % of rhythmic genes rhythmic in both sexes, 9 % only in males, 16 % only in females, 4 % with different parameters (model 5), which were enriched for core clock genes; "genes, namely Cry1, Npas2, Nr1d2, and Nfil3, showed a higher" amplitude in females.
- In GF mice sex-specific rhythms collapse: "a majority lost their sex-biased rhythmic expression and" became similar in males and females; male GF livers are feminised (genes up in GF male enriched for female-biased genes), female GF masculinised.
- Metabolome: "Around 34% of metabolites were found to oscillate in abundance in the liver and 20% in the serum"; sexual dimorphism of metabolite rhythms is damped in GF.
- Mechanism: altered sexual maturation and GH secretion; ghrelin injection rescues feminisation of GF males (Fig. 7E-H).

## Not addressed / limits
- Mouse, young adults; sex effect only within one strain and housing; the main rhythm data are of liver, intestine and WAT, not brain or SCN.
- Because the clock genes hardly change, this does not show that organ clocks desynchronise; it shows that output programs and metabolism rhythms change.
- Sex-specific rhythmic genes are a minority (about 29 % of rhythmic liver genes differ in rhythmicity by sex); effect on the core clock phase is small.
- Time points and replicates per time point are not retrievable from the extracted text.

## Link to us
- **Intro claim(s):** refs_circadian_desynchrony.md claim: microbiome required for sex-specific diurnal liver rhythms; sex effects in rodent rhythms. SUPPORTED (liver and WAT, GF vs ConvR); the refs file labels this Naef-lab coauthor, which is true (Naef is a co-author; correspondence is Gachon).
- **Relation to clocksync:** background for sex as a covariate of rhythm amplitude and for the finding that the molecular clock stays similar while output rhythms differ (supports reading the clock-gene panel as the shared state).

## Ideas
- Stratify clocksync coupling c and |z| by sex in GTEx (aggregates only) and compare with the female-higher amplitude of Cry1, Npas2, Nr1d2, Nfil3 here; specificity companion on the same permutations.
- Microbiome-rich compartments (colon, ileum, liver) vs others in GTEx: check whether offsets arg z differ, using the GF intestine and liver here as the mouse expectation (subtle clock phase change).
