# hughey2016zeitzeiger
**Ref:** Hughey JJ, Hastie T, Butte AJ, Nucleic Acids Res 2016;44(8):e80. DOI 10.1093/nar/gkw030. Type: method (supervised sparse learning for a periodic variable). Lab: Butte (UCSF) / Hastie (Stanford).
**Source:** full text

## Question and design
Method paper. Training and validation data are all mouse: GSE54650 (12 organs, 288 samples, every 2 h for 2 days, constant darkness after 12:12 LD entrainment, Affymetrix), plus 20 independent mouse datasets (nearly 800 samples) including phase-shifted (diet) and clock-mutant mice. No human data analysed in this paper (human applications are in later work, e.g. Wittenbrink 2018, Laing 2017).
**Sampling:** timed, full 24 h in training (48 h at 2 h spacing); validation sets cover different organs and circadian times, LD or DD.

## Method (if any)
ZeitZeiger: (1) time-dependent mean and variance of each feature by periodic spline, (2) sparse principal components (penalised matrix decomposition) of the time-dependent means, (3) Gaussian likelihood of each SPC given time, (4) maximum-likelihood prediction of time of a test observation. Output: predicted time (phase) and log-likelihood (fit to the expected trajectory). Training needs observations with known time spanning the cycle; time-points need not be evenly spaced; test data need no missing values for the selected features. The multi-organ predictor was built after ComBat treating each organ as a batch, and uses 13 genes (8 core clock genes) on average; 2 genes already give MAE about 1.1 h. Mouse tissues only; transfer to human tissue is not tested here.

## Key results
- Leave-one-organ-out CV (12 mouse organs): MAE 0.6-1.1 h for nSPC >= 2 (random predictor 6 h); abstract: "ZeitZeiger predicted circadian time (internal time of day) in each of 12 mouse organs to within ∼1 h".
- Versus molecular timetable: "Predictions by MT had an overall MAE of 1.2 h, 34% higher than those by ZeitZeiger".
- Per organ: cerebellum mean error -1.2 h (MAE 1.4 h), skeletal muscle +1.3 h (MAE 1.5 h); all other organs MAE below 1.2 h.
- Independent datasets: "the median absolute error was <1.5 h in 9 of 10 datasets".
- Phase shifts: high-fat-diet mice show errors about 2 h higher with similar log-likelihood ("they did not impair its operation"); genetic mutants show larger errors in 5 of 7 datasets and lower log-likelihood in 4 of 7; SPC1-SPC2 trajectory shrinks in liver Nr1d1/Nr1d2 knockout and shifts in muscle Arntl knockout.
- Design note: "If the oscillator is phase-shifted but otherwise functioning normally, observations will tend to have constant, non-zero prediction error and high log-likelihood."

## Not addressed / limits
- No human data and no disease; the claim that the predictor transfers to humans rests on later papers.
- Needs timed training data in the same tissue (or a multi-organ training set with batch correction); the multi-organ predictor still needs the test data merged and ComBat-adjusted with the training data, so it is not a single-sample tool without that step.
- Error is defined against external time; the authors stress that this assumes external time equals true circadian time (mice in DD drift).
- A predicted time cannot flag extrapolation: "the predicted value alone cannot indicate when the method is extrapolating."
- Per-organ errors are reported but no inter-organ phase offsets within an animal are estimated (organs are treated as batches to remove).

## Link to us
- **Intro claim(s):** supervised phase inference from timed data. Our refs line (supervised, timed training spanning the cycle, single sample, bulk) is mostly right but should say mouse multi-organ training, no human application in this paper.
- **Relation to clocksync:** background / competitor (supervised family). The phase-shift versus dysfunction logic (constant error with high likelihood versus variable error with low likelihood) parallels our offset arg z versus gate/coupling failure.

## Ideas
- Compare: ZeitZeiger log-likelihood as a per-sample dysfunction score against our coupling c / gate on the same mouse multi-organ data (GSE54650, per-organ phase offsets in the same animal).
- Mouse multi-organ GSE54650 gives organ-wise phase offsets from ZeitZeiger residual mean error (cerebellum -1.2 h, muscle +1.3 h): a reference for the sign and size of offsets we should recover.
- Missing reference: Ueda 2004 molecular timetable (the baseline here).
