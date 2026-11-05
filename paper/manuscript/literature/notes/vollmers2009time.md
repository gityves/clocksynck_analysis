# vollmers2009time
**Ref:** Vollmers, Gill, DiTacchio, Pulivarthy, Le, Panda, PNAS 2009. DOI 10.1073/pnas.0909591106. Type: biology (data: GEO GSE13093). Lab: Panda.
**Source:** full text

## Question and design
How much of the hepatic circadian transcriptome is driven by feeding rather than the liver clock? Mouse liver (Affymetrix arrays) from WT C57BL/6 and
oscillator-deficient Cry1-/-;Cry2-/- mice under ad libitum feeding, 8 h daytime-restricted feeding (tRF), and fasting (24 h from CT4 or CT14) with
refeeding; respiratory exchange ratio (n = 8 WT, 4 or 3 knockouts). **Sampling:** timed, 1 h (WT) or 2 h (knockout) resolution over one 24 h cycle
for ad lib and tRF; fasting arms 2 h resolution for up to 10 h after CT2 or CT12 (22 h of coverage in all). One tissue (liver); no SCN or other organ.

## Key results
- Feeding drives the daily RER and phosphorylation rhythms even without a clock: Cry1-/-Cry2-/- under tRF show food anticipation and an RER rhythm "indistinguishable from those of the WT" mice (Fig. 1).
- Most WT hepatic rhythms need food: "only 368 of the 2,997 rhythmic transcripts including the core oscillator components continue self-sustained rhythms in the absence of feeding" (Fig. 6C,D); more than 80% stop oscillating when food is withheld.
- Feeding time moves WT phases: "Feeding consolidation under tRF condition improved rhythms in a set of 3,217 probe sets, while dampening the rhythms of a set of 1,254" (Results); 84% of common cyclers (1,743) show a shift of more than 6 h in phase from night-fed ad libitum to daytime tRF.
- In a clockless liver, feeding alone restores rhythms: "feeding consolidation restored rhythmicity for 617 probe sets" in Cry1-/-Cry2-/- under tRF, none in ad libitum; most WT cyclers are only low-amplitude in the knockout, so the clock sharpens phase and amplitude.
- Rhythms of liver-Bmal1-knockdown mice (Kornmann 2007 set) are rhythmic in the clockless Cry-deficient liver under tRF: "nearly all transcripts previously shown to be rhythmic in the liver of liver-specific Bmal1 knockdown mice" regain rhythm (Fig. S3).
- Core clock Per genes respond acutely to refeeding: Per1 mRNA declines and Per2 rises, also in knockouts.

## Not addressed / limits
- Liver only: nothing here on other tissues, so no inter-tissue offset; the effect of feeding on tissue phase is measured in absolute ZT within the one organ.
- Mouse, young males, no light shift; Cry1/2 double knockout is whole body, so SCN output is lost too.
- Phase of clock genes under tRF is not tabulated against SCN or another tissue, so how many hours the liver shifts relative to what has to be read from Damiola 2000 / Stokkan 2001, not from here (the paper reports a shift of more than 6 h for common cyclers vs ad libitum).

## Link to us
- **Intro claim(s):** CLAIM WRONG IN DETAIL. refs_circadian_desynchrony.md says "Bmal1 KO vs feeding". The oscillator-deficient animals are Cry1-/-;Cry2-/-; liver-specific Bmal1 knockdown data come from earlier work and are only re-read here. The qualitative claim (feeding and clock are separate drivers of hepatic rhythms) is supported.
- **Relation to clocksync:** background and caution. A non-clock driver (feeding) shifts hundreds of liver rhythms by more than 6 h; relevant to hospital patients (fasting, nutrition changes) and to the cue/locked gene sets (A81b, A111, A112).

## Ideas
- Done: feeding-driven versus clock-locked liver response with mouse ortholog gene sets (A111 Manella, A112 Greenwell; A111 passed H111.1/H111.2 in part, A112 H112.1 failed with wrong sign, results_preliminary.md). Not yet done: GSE13093 itself (Cry-DKO under tRF) as a negative-control for the clock gate in a clockless but feeding-rhythmic liver (grep GSE13093 shows no hit in the docs).
- Reference we lack: Kornmann 2007 (PLoS Biol, in references.bib) is cited here for liver clock vs systemic cues; check it is carded.
