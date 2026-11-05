# koronowski2019defining
**Ref:** Koronowski et al., Cell 2019 (177:1448-1462). DOI 10.1016/j.cell.2019.04.025. Type: biology (tissue-specific reconstitution). Lab: Sassone-Corsi / Aznar Benitah.
**Source:** full text

## Question and design
How autonomous is the liver clock? Bmal1-stopFL knockout (KO) mice crossed with Alfp-Cre give hepatocyte-only Bmal1 re-expression (Liver-RE); WT, KO, Liver-RE compared. Liver RNA-seq (3 per time point), metabolomics (757 metabolites, n = 4), BMAL1 ChIP, then DD. **Sampling:** timed, six time points across a 12:12 LD cycle (exact ZTs not stated in the text read), ad libitum feeding; DD cohort sampled after 7 days in DD. Single tissue (liver, plus epidermis imported from Welz 2019); no paired multi-tissue samples.

## Method (if any)
JTK-type rhythmicity per genotype; autonomous = rhythmic in WT and Liver-RE; comparison of peak phase and amplitude. No inference of phase from untimed data.

## Key results
- Core clock genes oscillate in Liver-RE with advanced timing: "Thus, under physiological conditions, oscillations within the liver are independent from other body clocks, revealing tissue-level autonomy." BMAL1 phosphorylation and ChIP peak at ZT8 in Liver-RE vs ZT12 in WT, i.e. ~4 h earlier (Fig 1J,K); acetylation ~4 h advanced (Fig 5E).
- Transcripts: 2,010 Bmal1-dependent rhythmic in WT; "Of these, 218 (10.8%) also oscillated in Liver-RE mice" (Fig 3A). "In Liver-RE, 92.2% of transcripts were phase advanced" and 78.4% reduced in amplitude.
- Metabolites: 200 of 757 rhythmic in WT; "WT oscillatory metabolites retained oscillation in Liver-RE mice" for 38 (19%); glycogen and NAD+ salvage are autonomous.
- Epidermis-RE (from Welz): 392 of 2,606 rhythmic genes (15.0%) autonomous, half in phase with WT; only 22 autonomous genes (3.7%) shared with liver, so autonomy is tissue-specific.
- DD: "After 7 days in DD, the WT molecular clock oscillated robustly, in stark contrast to Liver-RE mice" (Fig 7); cosinor 20-26 h found no clock gene rhythm. Authors argue arrhythmia, not phase dispersion: "Thus, reconstituted livers are still coherent between animals in DD but do not exhibit circadian oscillations."
- Conclusion: "Hence, full circadian function in the liver depends on signals emanating from other clocks, and light contributes to tissue-autonomous clock function." (Summary)

## Not addressed / limits
- Liver-RE mice are behaviourally arrhythmic with flat feeding; the 4 h advance is relative to WT liver at the same ZT in LD, not to another tissue in the same animal.
- Group-mean sampling across different animals per time point cannot distinguish a non-oscillating clock from inter-animal phase dispersion; Sinturel 2021 (real-time recording) reaches the opposite conclusion for hepatocyte-only mice in DD.
- No other tissue profiled in the same mouse; no human data.

## Link to us
- **Intro claim(s):** inter-tissue relations and tissue-specific rescue as precedent for uncoupling: supported. Liver clock partially depends on extra-hepatic clocks and light.
- **Relation to clocksync:** background and a biological truth for a tissue whose clock is cut from the others: liver loses ~89% of Bmal1-dependent rhythms and runs ~4 h early. Related test done on a different line: A111/A129 (hepatocyte Bmal1 KO, Manella 2021: liver offset to WAT + lung host -1.31 h untimed, p 0.001).

## Ideas
- Done: knockout-tissue readout (A129, A111). The A129 scout judged Koronowski/Welz single-tissue and unused; GEO accessions not given in the text read, so the matrices are not confirmed open.
- Use the ~4 h liver-RE advance (and the preserved phase in Epidermis-RE, Welz 2019) as the expected sign/size of an offset if a paired multi-tissue Bmal1-stopFL set appears.
- Missing reference: Husse 2014 (SCN-specific Bmal1 loss, liver in DD), Izumo 2014 (SCN-less peripheral oscillation).
