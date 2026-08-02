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

## 1. Titan NH3 — JOINT no-ocean+ocean cache + SBI (Task #68 governs)

RECONCILED by Machine A 2026-08-02 per the user's firm decision (2026-07-30,
re-affirmed 2026-08-02) and Machine B's addendum
`plans/STATUS-2026-08-01-machineB-joint-nh3.md`. The earlier "Phase-0
rectangle re-run" framing here is WITHDRAWN: the provisional rectangle
Tb [248,257] K x w [30,100] ppt would re-condition the posterior on "an ocean
exists" and exceeds the 70 ppt NH3 ceiling. The corrected-activity-model
directive (`6c5ee2af`) COMPOSES with, not supersedes, the joint design — its
requirement is satisfied because the joint validations and cache are built at
a HEAD that includes the correction. Any pre-correction feasibility numbers
remain void.

- JOINT posterior over the FULL range Tb in [249, 263] K x w in [1, 70] ppt
  NH3. Frozen (Tb, w) nodes build as REAL no-ocean interiors via
  `build_tbw_grid_cache(..., retry_frozen_as_no_ocean=True)`, tagged
  `has_ocean=False`, and STAY in support — do not truncate Tb, do not drop
  frozen nodes.
- Config `configs/test54_titan_nh3_freegrav.json`: no
  `phase_stability.enforce` guard in either direction; density-inversion
  guard + k2 support bounds still apply.
- Continue Task #68 as planned: 3x4 joint validation -> scientific-reviewer
  interpretation -> author config -> USER go-ahead -> full-range production
  cache + 1M NH3 build + gates.
- Commit and push the builder edits (typed `NoIceLiquidTransitionError`,
  `do_overrides` channel, `retry_frozen_as_no_ocean` flag) once the 3x4
  validation confirms them — Machine A needs them on origin before wiring
  anything that reads the new cache metadata (`has_ocean`,
  `n_no_ocean_nodes`).
- Spec: `plans/active/titan-nh3-ocean-campaign-spec.md` (joint-posterior
  supersession note at top; per-phase mechanics otherwise still apply).

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
