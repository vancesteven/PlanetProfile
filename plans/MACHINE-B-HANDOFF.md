# Machine B handoff

Updated: 2026-08-01 at genai `54106fbd` (Machine A refresh after v5/v6/v7
delivery). Authoritative executable queue for compute-intensive work. Machine B
should pull the exact `origin/genai` commit named by Machine A before each
campaign, use the `PPcl` mamba environment, record package versions and seeds,
and return artifacts plus machine-readable reports. Never tune thresholds
after seeing a failure: stop and surface the evidence to Machine A.

Machine B model roster: Claude Opus 4.8 / Sonnet 4.6 / Haiku 4.5.

## 0. v5/v6/v7 gate adjudication COMPLETE (2026-08-02) — scoped follow-ups

Full record: `plans/active/europa-v5v6v7-gate-adjudication.md` (manager +
independent Opus-5 scientific review). Verdicts: v5 NOT ratifiable
(D_iceIh shape-excess survives HEAD, 1.12x; zeta_Ih mean 1.074x; SBC
underpowered n=108); v6 conditionally ratifiable (clean at HEAD; SBC n=102
below the >=200 spec minimum); v7 blocked — its crosscheck FAIL is
adjudicated a REFERENCE-side artifact (true v7 posterior provably identical
to v5's at the fiducial; mass in the newly opened support region measured
0.000000 in both references; the 1.06 km reference disagreement has the
wrong sign and rigid-translation shape of nested-sampling log-volume error;
which reference wandered is OPEN — v5's log_Z_err 0.281 makes it the more
suspect run). NO retraining and NO dataset regeneration — execute exactly
these follow-ups, in order, and stop at any surprise:

- **B1** Regenerate ALL gate reports (v5/v6/v7 baselines + arms) at current
  HEAD: use the DEFAULT constructed sweep grid (do not pass
  --sweep-values), record the validate_sbi.py commit SHA in every report,
  apply BH-FDR uniformly, reconcile `v{5,6,7}_gate_summary.json`.
- **B2** SBC re-runs at n_sbc_pairs >= 500 for v5/v6 baselines; report raw +
  BH p for every param, explicitly D_iceIh_km / log10_wOcean_ppt / Tb_K.
- **B3** Reference wander: re-run BOTH v5 and v7 reference MCMC, >=3 fresh
  seeds each, SAME pinned environment (note scipy 1.16.3→1.17.1 straddled
  the originals), n_live >= 2000. Report per-seed D_iceIh/log10_w means and
  log_Z ± err; success = between-seed scatter brackets the observed
  1.06 km. Do NOT assume v7 is the outlier.
- **B4** Evidence-ratio check: ΔlogZ vs −ln(V7/V5), V7/V5 by Monte Carlo
  over the prior under the two induction-bound sets (expected ~2.5).
- **B5** v7 flow-side diagnostics (preregistered; required before any v7
  ratification): (i) flow-vs-flow v5/v7 at the fiducial; (ii) v7-flow vs
  v5-reference formal crosscheck incl. shape clause; (iii) synthetic-truth
  bias test; (iv) local expected-coverage (TARP) near the fiducial.
- **B6** Limits anchor mode (W1 <= 0.25 sigma vs ground-truth MCMC) at
  reachable Im_k2 values, all three campaigns.
- **B7** Record outcomes as FAIL-ADJUDICATED-ACCEPTABLE where applicable;
  never relabel a FAIL to PASS.

Priority vs Titan: Titan NH3 Task #68 (section 1) remains your active task;
run B1–B6 after the 3x4 validation completes, before the full NH3
production build.

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
