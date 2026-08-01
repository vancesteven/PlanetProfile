# Machine B handoff

Updated: 2026-08-01 at genai `54106fbd` (Machine A refresh after v5/v6/v7
delivery). Authoritative executable queue for compute-intensive work. Machine B
should pull the exact `origin/genai` commit named by Machine A before each
campaign, use the `PPcl` mamba environment, record package versions and seeds,
and return artifacts plus machine-readable reports. Never tune thresholds
after seeing a failure: stop and surface the evidence to Machine A.

Machine B model roster: Claude Opus 4.8 / Sonnet 4.6 / Haiku 4.5.

## 0. HOLD — v5/v6/v7 gate adjudication (Machine A)

v5, v6, and v7 are DELIVERED and received (artifacts, caches, reference
results, gate reports all on origin/genai — thank you). All three baseline
arms share the same failure signature (sbc pass; limits containment fail;
crosscheck fail on D_iceIh_km + log10_wOcean_ppt). This is under Machine A
scientific adjudication. Do NOT retrain, regenerate datasets, or adjust any
gate for these campaigns until Machine A sends a scoped verdict. Machine A may
send follow-up diagnostics requests (e.g. extended reference-MCMC chains or
containment sweeps) — treat those as the next v5/v6/v7 action, nothing else.

## 1. Titan NH3-ocean campaign — Phase 0 RE-RUN (top executable item)

The NH3 water-activity model was corrected 2026-07-28 (commit `6c5ee2af`,
Melinder-anchored; the CoolProp mixture excess had the wrong sign and the
"L-K 2002" shift was deleted). Any earlier Phase-0 feasibility scan is void —
it would bake a 24–32 km D_iceIh bias into the training set.

- Full executable spec: `plans/active/titan-nh3-ocean-campaign-spec.md`
  (read the "IMPORTANT — activity model corrected" section first).
- Prereqs: pull genai >= `6c5ee2af`; `pip install CoolProp` in PPcl;
  sanity `python -m pytest tests/coolprop_nh3_test.py` (15 pass, ~3 min).
- Provisional Phase-1 rectangle: Tb in [248, 257] K x w in [30, 100] ppt —
  confirm with the Phase-0 scan before building the 2D cache.
- Then Phases 1–3 per spec (2D cache -> config + reference MCMC -> SBI +
  gates). Send the Phase-1 rectangle and the production config to Machine A
  for freeze before generating the production dataset.

## 2. Cassini–Enceladus production

Status: queued after Titan Phase 0. Prior physics blockers closed
(hydrostatic ratio 3.25, correlated C20/C22, negligible elastic libration).

- base config: `PlanetProfile/Inference/configs/enceladus_cassini_smoke_6D.json`;
- constraints: `plans/mission-body-constraints.md`;
- dense Seawater 10 ppt Tb cache over ~271.8–272.6 K;
- cache-build overrides: `Cuncertainty = CuncertaintyUpper =
  CuncertaintyLower = 0.08`;
- production adds `dC20_nh ~ U[-1.5e-3, 1.5e-3]`, `dC22_nh ~ U[-1e-4, 1e-4]`;
- retain `observable_correlations = {"C20,C22": 0.47}`,
  `gravity_j2_over_c22 = 3.25`; libration is the primary decoupling channel.
- Phases: dense cache + inspection -> config freeze via Machine A ->
  reference MCMC + SBI + SBC/crosscheck/limits -> Ae sidecar -> ship all
  with seeds/versions/verdict. Do not deploy.

## 3. Callisto campaign (roadmap — no frozen config yet)

- Read `papers/cochrane2025stronger.pdf` before finalizing induction
  observables/uncertainties.
- Existing Callisto C2 structure cache has a P_MPa/r_m length mismatch
  (118 vs 122) — regenerate during setup.
- Do NOT inherit Europa's |Ae| support cuts: the Callisto prior must keep
  no-ocean/low-|Ae| AND saturated high-|Ae| corners inside support
  (Cochrane2025 posture; see 2026-07-25 audit note in git history of this
  file for the numbers).
- User directive: full coverage in ice thickness and ocean COMPOSITION;
  SeaFreeze 1.1.3 ships NaCl(aq) splines (NaClaq_LP/HP) — consider NaCl
  cache variants.

## Environment notes

- SeaFreeze pinned 1.1.3 (`pip install --upgrade SeaFreeze==1.1.3`); code at
  HEAD does not run against 1.0.0.
- CoolProp >= 8.0.0 required for NH3 work.
- Frozen reference env: `environment.yml` (PPcl, committed `11b514dd`).

## Queue boundary

No heavy campaign beyond the above is authorized. Anything else is roadmap
material until Machine A provides a frozen config and an explicit handoff.
Any `PlanetProfile/Test/` additions or `.gitignore` changes require explicit
user permission before commit.
