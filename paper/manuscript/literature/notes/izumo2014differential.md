# izumo2014differential
**Ref:** Izumo et al., eLife 2014 (3:e04617). DOI 10.7554/eLife.04617. Type: biology (ex vivo real-time bioluminescence). Lab: Takahashi (with Bradfield).
**Source:** full text (eLife article JSON) and source data (Figures 5, 8, 9)

## Question and design
What happens to peripheral clocks when the brain clock is removed but the animal is otherwise normal? CamKIIa-Cre;
Bmal1 fx/fx (BKO: forebrain and SCN) vs Fx/Fx, crossed to PER2::LUC. Explants of pituitary, liver, kidney, heart, lung,
spleen from the same mouse (n = 11-15 per arm), recorded in LumiCycle; arms LD, DD (> 30 days), DD + restricted feeding
(FR, CT/ZT 6-10, 2 weeks). Heart imaging on 40 um grids (3 pairs). **Sampling:** continuous explant recording;
phase = peak of the second cycle; all tissues of a mouse cut together at ZT/CT 10.5-13.

## Method (if any)
FFT-NLLS on detrended traces (first 20 h trimmed); circular statistics across mice per organ (Watson-Williams,
bootstrap variance) and, per mouse, circular variance over its six organs (Mann-Whitney between genotypes).

## Key results
- Behaviour arrhythmic in DD and LL; SCN explants blunted.
- LD keeps the organ map: "BKOs displayed well-phased rhythms from tissue to tissue when maintained in LD" (Fig 5F).
- DD: "a loss of phase coordination occurred in peripheral clocks within individual BKO animals when placed in DD" (Fig 5G, p = 0.0008).
- Amplitude: "the relative amplitude was significantly decreased when a light-driven signal was removed from BKOs in DD, but the severe reduction was not observed in BKOs from LD" (Fig 6A,B). Heart imaging: "the amplitude decrease in the BKO tissue is caused by a combined effect of amplitude reduction and phase desynchrony at the cellular level" (Fig 6H,I).
- FR: liver and kidney re-synchronise at a new phase; heart, lung, spleen shift but cluster less; within-mouse variance stays high (Fig 8H, p = 0.0002). Amplitude stays low except liver: "with the exception of the liver, oscillatory amplitude of peripheral rhythms of BKOs remained lower than that of Fx/Fx control mice even after 2 weeks of FR" (Fig 8I).

## Not addressed / limits
- Explant phases include the dissection reset; BKO phases in DD cluster across arrhythmic mice, so some common cue (dissection, housing) remains.
- Within-mouse circular variance mixes a consistent but spread organ order with independent organs (our A319 separates them).
- Amplitude only as bar graphs (no source data); FFT-NLLS amplitude on explants, not in vivo.
- Forebrain deletion is broad (DMH, pituitary partly), not SCN-only.

## Link to us
- **Intro claim(s):** the central clock keeps the organs' order; a substitute drive (light without SCN, feeding) sets an order of its own.
- **Relation to clocksync:** supports (mouse counterpart of A317 / A315); data we used (A319). Tension: tissue amplitude falls after > 30 days without any shared drive, against "amplitude kept" in brain-dead donors (hours to days, ICU cues).

## Ideas
- A319: reanalysis of the source data (order variance vs an independence null).
- Per-mouse per-organ phases and amplitudes are not public: a request to the authors would allow a paired test (outward-facing: ask the user first).
