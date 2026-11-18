# grabe2022synergies
**Ref:** Grabe, Mahammadov, Del Olmo, Herzel, Front Netw Physiol 2022 (1:803011). DOI 10.3389/fnetp.2021.803011. Type: model (theory). Lab: Herzel (Charite / HU Berlin).
**Source:** full text

## Question and design
How do two coexisting 24-h zeitgebers acting on different clock nodes combine? Light on Per (square wave, amplitude r_light) and a Rev-Erb agonist drug on Rev-Erb (r_+REV), phase difference DeltaPhi = difference between light onset and drug onset. Two mouse-based gene-regulatory models: a modified Korenčič model (15 ODEs from the 5-DDE liver "LD" model, free running period 24.95 h in Methods and 24.8 h in Fig. 3) and the Almeida model (7 proteins). **Sampling:** simulations; 80 zeitgeber cycles, last 10 analysed by Lomb-Scargle.

## Method (if any)
Systematic scans of zeitgeber period, strength, DeltaPhi; entrainment range, phase of entrainment (Per crossing its mean relative to light onset), entrained amplitude, Arnold tongues, phase-response-curve (PRC) overlay, transient time.

## Key results
- Single zeitgeber: at r = 0.2 the range of entrainment is 1.7 h for light and 3.6 h for the Rev-Erb agonist; a 12 h rule holds (phase spans about 180 degrees across the range); "large periods lead to earlier phases"; stronger zeitgeber gives a weaker slope of phase vs period; amplitude peaks at 25 h (light) and around 26 h (drug).
- Synergy depends on DeltaPhi: there are "synergistic phase differences between the two zeitgebers, which result in larger ranges of entrainment than for solitary zeitgebers, and antagonistic phase differences with small ranges of entrainment". Synergy follows from overlap of the PRCs of the two solitary zeitgebers: optimal when the +REV PRC is shifted about 15.5 h earlier (DeltaPhi = 15.5 h), antagonism near 4 h.
- Numbers: DeltaPhi = 12 h vs 3 h: "1.5-fold higher amplitude" and range 22.2-28.2 h vs 23.8-26.8 h; transient to entrainment median 9 cycles at 15.5 h vs 20 at 4 h; at DeltaPhi = 4 h strong zeitgebers (r > 0.2) are needed for entrainment, at 15.5 h r = 0.1 suffices.
- The entrained phase mostly follows the drug cycle in the Korenčič model; in the Almeida model phase is mixed, with a sign change near DeltaPhi = 9 h and period doubling near 14 h.
- Different periods (24 h light, 28.8 h drug): 5:6 locking with a folded limit cycle and amplitudes that "vary dramatically from day to day".
- Discussion proposes the phase relation between Rev-Erb modulation and light as a way to strengthen low-amplitude rhythms (septic shock, cancer).

## Not addressed / limits
- Both zeitgebers are applied to different nodes (Per and Rev-Erb) of a single cell-autonomous oscillator; cues on the same node, or an effective single input, are not studied. Synergy is a statement about relative phase in the PRC sense (15.5 h optimum), not about "12 h apart" per se.
- Zeitgeber amplitudes are equal (r_light = r_+REV = 0.2 in the main figures); one drug (Rev-Erb agonist) and square waves; no cell-level or organ-network coupling; no data. The authors note there is "not a gold standard amplitude mesor", so amplitude is hard to define.
- Entrainment of oscillators above the Hopf point (self-sustained limit cycles), not damped or driven states.

## Link to us
- **Intro claim(s):** conflicting zeitgebers shift or damp peripheral clocks (cites Manella 2021, Mukherji 2015); amplitude after conflicting cues depends on how cues combine.
- **Relation to clocksync:** background and a caveat to our drive-cancellation assumption. In their model synergy vs antagonism is set by relative phase and by which nodes receive the cues: 12 h apart on different nodes (light on Per, agonist on Rev-Erb) is near synergy (1.5-fold amplitude vs 3 h), while the antagonistic DeltaPhi is about 3 to 4 h. Cancellation of two cues 12 h apart, as in our one-effective-input model, requires that both cues act on the same node with the same PRC (a single input); with two nodes at 12 h the result depends on the PRC overlap.

## Ideas
- State in the model that the cancellation assumption holds only if central and behavioural cues feed one effective input; a night shift (muscle amplitude to 0.23, shift 2.7 h) tests this since muscle data show lower amplitude and not synergy.
- Check with their Almeida and Korenčič models whether feeding on one node and light on another gives amplitude reduction at 12 h for some node pairs (PRC overlap), to bound where cancellation is expected.
- Related references we may lack: Nikhil 2014 (flies, light plus temperature, 4 h optimum), Daminelli 2014 (transitions between dominating zeitgebers), Manella 2021 (already in literature).
