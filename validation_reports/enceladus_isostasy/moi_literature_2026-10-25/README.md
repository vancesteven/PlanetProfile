# Enceladus MoI and core-density literature check — 2026-10-25

Sources (in papers/): iess2014gravity(+SOM), mckinnon2015effect, beuthe2016enceladus(+SOM),
beuthe2018enceladus, hemingway2019enceladus. Text extracted with `pdftotext -layout` (Sonnet agent);
page, section and equation locations are as reported by the agent. Values are quoted, and interpretation
is marked as such.

## Moment of inertia

**Iess 2014**
- Source (SOL1, E9/E12/E19): J2 = 5435.2 ± 34.9e-6 and C22 = 1549.8 ± 15.6e-6 (1σ).
- How α ≈ 0.335 was obtained: from an Airy-compensation split of the degree-2 gravity and topography, with f30 fixing the compensation depth and assuming f20 = f22, plus the Radau-Darwin linearisation (SOM Eq. S7, which imposes 10/3).
- Uncertainty quoted by the paper (admittance test, SOM S3.5 / Fig. S6):
  - 1σ: α = 0.333–0.338
  - 2σ: α = 0.330–0.340
- **There is no ±0.001 in the paper.** The template's `Cuncertainty = 0.001` has no source.

**McKinnon 2015**
- Under rapid rotation the hydrostatic J2/C22 is about 3.25, not 10/3.
- Refit result: NMOI ≈ 0.330 (range 0.328–0.333).
- Z22 changes sign at NMOI ≥ 0.336, which he treats as a quasi upper limit.
- With a 925 kg/m³ shell: ρc ≈ 2450 (range 2375–2500) kg/m³.
- With a 1000 kg/m³ shell: 0.3305 ± 0.0025 and ρc = 2530 ± 70 kg/m³.

**Beuthe 2016**
- Gives no MoI.
- From the MCMC fit to SOL1 gravity plus TOPA shape: ρc = 2422 ± 46 kg/m³ (1σ). Case range in SOM Table S6: 2370–2479.

**Hemingway & Mittal 2019**
- "we cannot obtain the moment of inertia directly" (Sect. 3.1).
- Airy + Tajeddine shape: ρc 2340–2410 kg/m³ (68%), shell 19–24 km, ocean 35–39 km.

**Park 2024 (erratum)**
- 0.337 (0.335–0.338).

**Every one of these values is derived from the same Cassini J2/C22, under hydrostatic-core plus isostatic-shell assumptions.** None is independent of the data this campaign conditions on.

## Consequences for the campaign

1. Fixing moi_target = 0.335 double-counted the gravity data. This is now replaced by sampling ρ_s (`mass_rho_rock`, commit bfb227e1).
2. The ρ_rock ~ U[2200, 2600] box:
   - It contains every gravity-derived ρc range above (2320–2590). Its derived C/MR² of 0.323–0.347 also contains every MoI value quoted above.
   - **None of these papers supports the 2200 floor as an independent physical bound.**
   - Physical support has to come from outside these papers. Candidates: hydrated-silicate grain densities (serpentine/chlorite, about 2.5–2.7 g/cc) combined with core porosity (e.g. Choblet 2017, Neveu, Le Gall). These are not yet in papers/.
   - The only independent statements found are qualitative: H&M 2019 call their densities "suggestive of hydrated silicates and/or some degree of porosity", and the Iess SOM gives 2.5 g/cc as its low-density core example.
3. The PlanetProfile template `PPEnceladus.py` still sets Cuncertainty = 0.001. The Enceladus inference config overrides this with 0.08 for the cache build. Changing the template default is a PlanetProfile-wide choice, so it is flagged here and not changed.
