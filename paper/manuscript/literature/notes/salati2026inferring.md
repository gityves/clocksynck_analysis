# salati2026inferring
**Ref:** Salati, Paychère, Hahaut, Gobet, Naef, bioRxiv (posted April 1, 2026). DOI 10.64898/2026.03.30.715278.
Type: method (scRitmo). Lab: Naef (EPFL), our lab.
**Source:** full text

## Question and design
Infer single-cell circadian phase with uncertainty from scRNA-seq and separate biological desynchrony from technical
noise. Data: 695 unsynchronised 3T3 fibroblasts (FLASH-seq + Rev-Erbα-YFP reporter), SABER-FISH 10-55 h after Dex,
mouse liver, aorta, skin time series (LD, ZT-sampled), Drosophila clock neurons LD vs DD.
**Sampling:** time series (ZT known) for validation; inference itself is unsupervised.

## Method
Negative binomial counts, harmonic mean profiles; gene peak phases and amplitudes fixed from bulk multi-organ time
series (Clock-Reference, 15 genes) or extended per cell type (70-100 genes); marginal likelihood over latent phase;
posterior mode + circular SD per cell. σbio from "simulation-calibrated variance decomposition" (σdata vs σtechnical
simulated for a perfectly synchronised population).

## Key results
- Fixing gene phases: "fixing these parameters based on established biological knowledge often provides a useful constraint that stabilizes inference".
- Low depth bias: at 3,000 UMI/cell a "phase attractor" near the Bmal1 trough ("circadian counts are approximately 4-fold lower").
- Reporter validation: Nr1d1 mRNA vs YFP R2 = 0.73; reconstructed protein cycle "(R2 = 0.55)".
- Desynchrony grows after Dex: "the circular standard deviation grows by approximately 0.1 hours (5.4 minutes) per elapsed hour"; period spread "(σT = 1.78h)".
- In vivo: "For fibroblasts, smooth muscle cells, and hepatocytes, σbio ranged between approximately 1.1-1.6 hours".
- Depth: "at sufficient depth (>20,000 UMI/cell), the Clock-Reference alone provides accurate and unbiased phase inference".
- Limits stated: "the estimated σtechnical can exceed σdata , rendering the subtraction undefined"; relative comparisons are the robust use.

## Not addressed / limits
- Within-sample, between-cell desynchrony; not between organs or compartments of one person.
- Validation in vivo needs ZT time series (extended sets selected against known collection time).
- Drosophila pooled animals: intra- vs inter-individual dispersion not separable (stated).

## Link to us
- **Intro claims:** desynchrony is multi-scale (cells within tissue, organs within body); damping at population level reflects dephasing of intact oscillators; separating technical noise is essential (same logic as our Poisson V_ik term vs compartment noise σ_k²).
- **Relation:** our lab's companion method at cell scale; clocksync is the donor/compartment scale. Not a competitor.

## Ideas
- Within-compartment dephasing (σbio ≈ 1.1-1.6 h) lowers pseudobulk amplitude by a factor exp(-σ²/2) in radians: 0.96 (1.1 h) to 0.92 (1.6 h), small; state this in the methods to pre-empt reviewers.
- A joint figure: cell-level σbio (scRitmo) vs compartment-level coupling (clocksync) in the same human scRNA atlas, where depth allows.
