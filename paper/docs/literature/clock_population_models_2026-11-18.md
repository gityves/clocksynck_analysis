# Population models of tissue clocks: what is established (literature check, 2026-11-18)

Stage 0 of the clock-model task (`docs/prompts/clock_model_2026-11-18.md`). Search by a Sonnet subagent (about 40 web
calls, PubMed / PMC / publisher pages / arXiv; no Google Scholar, no systematic forward-citation search), verified by the
lead: every DOI below was resolved on Crossref (journal, year, first author, title match). "Passage read" means the full
text or the relevant passage was opened; "metadata" means the bibliographic record and abstract only. Absence of a
source is not proof of absence.

Candidate model: bulk clock state of cell class c in organ o of person p, x = a_c R_o exp(i(theta_p + phi_o + delta_op)),
with a_c the single-cell amplitude, R_o the organ's synchrony, phi_o a locked lag, delta_op a fixed-size additive
departure (angle ~ 1 / (a R)).

## (i) Bulk amplitude = single-cell amplitude x synchrony

| claim | source | evidence | new for us |
|---|---|---|---|
| Ensemble damping of cultured fibroblasts is loss of synchrony, not single-cell damping | Welsh, Yoo, Liu, Takahashi, Kay, *Curr Biol* 2004, 14:2289, doi:10.1016/j.cub.2004.11.057 | passage read: "rapid ensemble damping is much better explained by loss of synchrony among cells than by damping of individual cells" | the founding case of the factorisation (cell culture, after a reset; not tissue in vivo) |
| Fibroblasts carry self-sustained, cell-autonomous clocks; mitosis shifts phase | Nagoshi, Saini, Bauer, Laroche, Naef, Schibler, *Cell* 2004, 119:693, doi:10.1016/j.cell.2004.11.015 | metadata | supports a_c as a cell property; not itself the desynchrony claim (correction to the prompt) |
| Single fibroblasts keep their period and amplitude over six weeks | Leise, Wang, Gitis, Welsh, *PLoS ONE* 2012, 7:e33334, doi:10.1371/journal.pone.0033334 | metadata | stability of a_c in culture |
| Damped noise-driven and self-sustained single-cell models both fit single-cell traces | Westermark, Welsh, Okamura, Herzel, *PLoS Comput Biol* 2009, 5:e1000580, doi:10.1371/journal.pcbi.1000580 | metadata | caveat: "self-sustained" single cells is not settled from traces alone |
| Peripheral cells are weakly coupled; a stochastic phase-population model fits bulk bioluminescence (dispersion, noise, coupling) | Rougemont, Naef, *Mol Syst Biol* 2007, 3:93, doi:10.1038/msb4100130 | passage read: "coupling is found to be weak and insufficient to synchronize cells" | the closest population model of peripheral cells (one cell type, culture) |
| Population amplitude depends on cellular amplitude and the population's synchronisation | Masuda, Kon, Iizuka, Fukada, Sakurai, Hirano, *Nat Commun* 2023, 14:2819, doi:10.1038/s41467-023-38392-x | passage read (qualitative; no product formula) | the statement in words, cultured cells |
| Population-level and single-cell amplitude changes are distinguishable in reporter data | St John, Taylor, Abel, Doyle, *Biophys J* 2014, 107:2712, doi:10.1016/j.bpj.2014.10.026 | passage read | the separation is measurable in time series; our setting has none |
| In explants, single-cell amplitude and coherence are reported separately; coherence fades | Lande-Diner, Stewart-Ornstein, Weitz, Lahav, *Mol Biol Cell* 2015, 26:3940, doi:10.1091/mbc.E15-06-0403 | passage read | tissue (explant) version of Welsh 2004 |
| Cells gain amplitude as they synchronise (SCN) | Webb, Taylor, Thoroughman, Doyle, Herzog, *PLoS Comput Biol* 2012, 8:e1002787, doi:10.1371/journal.pcbi.1002787 | passage read: "As the system synchronizes, the cells gain amplitude" | **caveat for the model:** with coupling, a_c and R_o are not independent |

Not for (i): Bieler 2014 and Feillet 2014 concern clock / cell-cycle locking (see iv), not damping against desynchrony.

## (ii) Amplitude-dependent entrainment and phase of entrainment

| claim | source | evidence | new for us |
|---|---|---|---|
| The entrainment range is set by the ratio of zeitgeber strength to oscillator amplitude and by rigidity; coupling raises amplitude and rigidity; weak tissue clocks (lung) entrain to extreme T-cycles, SCN does not | Abraham, Granada, Westermark, Heine, Kramer, Herzel, *Mol Syst Biol* 2010, 6:438, doi:10.1038/msb.2010.92 | passage read: "the ratio between the zeitgeber strength and oscillator amplitude indeed determines the entrainment range" | the theory behind "a fixed-size input turns a weak clock more" (our angle ~ 1 / (a R)) |
| The phase of entrainment is monotonic in the period mismatch and depends on the zeitgeber / amplitude ratio; small-amplitude rhythms have shorter transients | Bordyugov, Abraham, Granada, Rose, Imkeller, Kramer, Herzel, *J R Soc Interface* 2015, 12:20150282, doi:10.1098/rsif.2015.0282 | passage read; the term is "Arnold tongue", not "Arnold onion" (correction) | phi_o follows detuning and gain: an organ's lag can be read as its detuning, if periods were known |
| Time to synchronisation | Granada, Herzel, *PLoS ONE* 2009, 4:e7057, doi:10.1371/journal.pone.0007057 ("How to achieve fast entrainment? The timescale to synchronization") | metadata | correction: "Granada 2009" in *Methods Enzymol* 454 (doi:10.1016/S0076-6879(08)03801-9) is a chapter on phase-response curves |
| Resonance (amplitude peaking at a detuning) in entrained circadian oscillators | mentioned only in passing in Bordyugov 2015 | not verified as a circadian result | the "no resonance" wording rule stands |

## (iii) Population, Kuramoto and mean-field models of peripheral clocks

| source | content | makes our factorisation? |
|---|---|---|
| Huang, Zhang, Braun, arXiv:2212.13579 (2022, preprint), "A minimal model of peripheral clocks reveals differential circadian re-entrainment in aging" | two Kuramoto phase populations (SCN, one periphery) with Daido mean-field reduction; "collective amplitude" = phase coherence; ageing as lower coupling / light sensitivity | no (full text checked by the lead): phase oscillators, so no single-cell amplitude factor; one generic peripheral population; no human tissue data |
| Hannay, Booth, Forger, *J Biol Rhythms* 2019, doi:10.1177/0748730419878298; Hannay, Forger, Booth, *Sci Adv* 2018, doi:10.1126/sciadv.1701047 | Ott-Antonsen macroscopic models: the SCN's amplitude is its order parameter; fitted to human light-response data | no: central clock only, amplitude = synchrony with no cell-amplitude factor (metadata) |
| Burckard, Teboul, Delaunay, Chaves, *Interface Focus* 2022, 12:20210087, doi:10.1098/rsfs.2021.0087 | coupled network of liver clocks (piecewise-linear) | no |
| Delaunay, Fortuné, Krawczyk, Anastasia et al., arXiv:2607.21214 (2026, preprint) | Kuramoto model of hepatocyte spheroids | no (abstract) |

## (iv) Turnover and cell-cycle coupling

| claim | source | new for us |
|---|---|---|
| Clock and cell cycle lock 1:1 in single cells; division at a fixed clock phase | Bieler, Cannavo, Gustafson, Gobet, Gatfield, Naef, *Mol Syst Biol* 2014, 10:739, doi:10.15252/msb.20145218 (metadata); Feillet et al., *PNAS* 2014, 111:9828, doi:10.1073/pnas.1320474111 (metadata) | daughters inherit phase: proliferation need not dephase a population; it may damp or rewire instead |
| Low-dimensional torus dynamics of the two coupled oscillators | Droin, Paquet, Naef, *Nat Phys* 2019, 15:1086, doi:10.1038/s41567-019-0598-1 (metadata) | same |
| Epidermal stem cells sit at opposite clock phases (heterogeneity created by the clock) | Janich et al., *Nature* 2011, 480:209, doi:10.1038/nature10649 (metadata) | a known within-tissue phase spread, relevant to R in epidermis and psoriasis |
| A tissue renewal model predicting bulk clock amplitude from turnover | **not found** | open; "Gaucher 2018" not found (probably a review) |

## (v) Human multi-tissue work

| source | what it claims |
|---|---|
| Ruben et al., *Sci Transl Med* 2018, 10:eaat8806, doi:10.1126/scitranslmed.aat8806 (correction: not *Science*) | CYCLOPS on GTEx; tissue-specific rhythmic genes; amplitude used as a threshold; no population model |
| Talamanca, Gobet, Naef, *Science* 2023, 379:478, doi:10.1126/science.add0846 (CHIRAL) | GTEx phases across tissues; core clock synchronised across the body; low brain amplitudes; sex and age effects; no coupling or order-parameter model |
| Wucher, Sodaei, Amador, Irimia, Guigó, *PLoS Biol* 2023, 21:e3001986, doi:10.1371/journal.pbio.3001986 | day-night and seasonal variation from time of death (correction: not CHIRAL) |
| Saint-Antoine, Anafi, bioRxiv 2025, doi:10.1101/2025.11.14.688559 (abstract read) | in human skin time series, rhythm amplitudes cluster by subject (a person-level global amplitude), and a person's epidermis amplitude predicts their dermis amplitude | a person-level gain factor shared between tissues: our candidate model has none (theta_p carries phase only); it must be added or excluded |

## Verdict

No published model found writes the organ clock as single-cell amplitude x organ synchrony at a locked organ phase,
spans many organs, or is fitted to human tissue. The parts are established separately:

- the factorisation in culture and explants (Welsh 2004; Lande-Diner 2015; Masuda 2023);
- the 1 / amplitude sensitivity to a fixed input (Abraham 2010; Bordyugov 2015);
- phase-only Kuramoto models for one peripheral population (Huang 2022) or the SCN (Hannay 2018, 2019).

The stop condition (a published model already making the factorisation in human tissue) is not met.

Points the model must handle:
1. **Coupling makes a_c and R_o dependent** (Webb 2012; Abraham 2010). A product a_c R_o with independent factors is a
   weak-coupling approximation. In strongly coupled tissue, a "same single-cell amplitude across organs" prediction
   can fail without the model being wrong.
2. **A person-level gain** shared between tissues is reported in skin (Saint-Antoine and Anafi 2025, preprint).
3. **Proliferation locks rather than dephases** (Bieler 2014; Droin 2019). The tumour gain loss (A304) is therefore
   not dephasing by default.
4. **Self-sustained single cells** are not settled from time series (Westermark 2009). A297's tissue-level result
   is consistent with them but is not proof.
