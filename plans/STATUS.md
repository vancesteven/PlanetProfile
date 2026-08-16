# PlanetProfile genai — current status

Big-picture program map: `plans/STRATEGY.md` (production line states,
what current work buys). This file is the tactical log.

Updated: 2026-08-16b (Machine A: **r5 = NOT RATIFIED — the gate battery
caught wiring-vs-API gaps before they became a built campaign.** All
physics verified sound (every constant/bracket/invariant reproduces);
the defects are PLUMBING: the ruled B2' treatment + correction +
sys_frac + rho_ice override exist as APIs but are not wired at the
production call sites (D3-D6, all on the sole thickness channel), the
ocean MoI window is an unrecorded builder default (D1), and the SBI
dataset generator cannot emit the observable vector (D2, no C30 arm).
Asymmetry RULED bound-and-caveat with B-A1/A2/A3 measurements required
in B's run design; Tajeddine double-count correction binding. Full
adjudication: validation_reports/enceladus_isostasy/r5_ADJUDICATION.md.
Path: D1 -> build can start; D3-D6 repaired together + verified through
the production dispatch; D2 before datasets; r6. A8 landed (f8ed04cd,
86+1xfail). C12 sweep FAIL routed as preregistered scope-bound outcome
(§0.25). Next session: D1-D8 repair batch (agents), then r6.)
Prior: 2026-08-16 (Machine B: **C12 amortized sweep for the deployed NH3 artifact: FAIL (29/149, 19.5%, threshold <=5%).** No driver for this existed; built one (plans/scripts/nh3_c12_amortized_sweep.py) after Opus flagged and I fixed a design trap first -- the likelihood closure is captured from config.observables at MCMCRunner construction and does NOT track a swept x_obs argument, so a naive sweep would have silently computed p(x_fiducial|theta)/q(theta|x_swept), a meaningless ratio with plausible-looking Pareto-k. Fixed by rebuilding a fresh config+runner per point. Preregistered every design choice (validation_reports/titan_freegrav_nh3_1m/is_correction/c12_sweep_PREREGISTRATION.md) before running; verified the driver reproduces the committed fiducial Pareto-k/ESS to full precision first. Result: 141/200 requested prior-predictive draws survived the support guard; 26/141 (18.4%) failed Pareto-k<=0.7, plus 3/8 (37.5%) axis endpoints including a hard C4 reject at the C22-hi corner. Failures concentrate at high z-distance from the training center (failed mean z=3.48 vs passed 1.59), overwhelmingly driven by extreme Re_k2/Im_k2 -- the same channel already under separate quarantine discussion, not a new C20/C22 finding. Does NOT affect the already-passed fiducial validations (C3/C5.3/C10/C11/C13/C16). Full findings + reviewer questions: validation_reports/titan_freegrav_nh3_1m/is_correction/c12_sweep_FINDINGS.md. Not self-adjudicated -- routed to the scientific-reviewer; does this bound scope (informative-only outside the already-narrower GUI x_obs_limits) or extend the k2 quarantine. Enceladus: confirmed A2-A7 landed + manager-verified while I was catching up (frozen-branch redesign); my next Machine B step (production cache -> 1M -> training, MACHINE-B-HANDOFF.md sec 0.24 item 4) stays gated on Machine A's r5 ratification pass, not yet run.)
closed: the deployed Titan free-gravity no-ocean campaign now has a generated,
source-backed methods snippet covering its 12 priors, 4 observables, 9-node
cache and SHA, 877,227-simulation artifact metadata, raw gate FAILs, and
committed adjudication. `all --verify` passes all nine campaigns; the new key
also verifies outside the repository working directory. Implementation
`2328f12b`. No queued Codex task.)
Prior: 2026-08-15c (Machine A Codex: **C24 VERIFIED.** Doc-doctor pass #3
records 11 PASS / 1 FINDING across the full 12-item checklist. The finding is
E11 only: the eight configured methods snippets verify, but the separately
deployed Titan free-gravity no-ocean campaign has no generated snippet. Slot
provenance 14/14 and results-panel tests 11/11 pass; all eight methods outputs,
README links, and the 8-cache deploy build-only snapshot verify. Report:
`validation_reports/doc_doctor/2026-08-15_third_pass.md`; implementation
`8f6cf33c`.)
Prior: 2026-08-15b (Machine A Codex: **C23 VERIFIED.** The
Cassini-Enceladus Seawater/isostasy analysis now has a stable, nondefault
awaiting-artifact slot that renders its frozen-design scope and disables
posterior generation without touching missing files. Focused AppTest 1/1 and
neighboring slot AppTests 2/2 pass. Implementation `856f50e6`; C24 next.)
Prior: 2026-08-15a (Machine A Codex: **C22 VERIFIED.** The permanent
Delta_rho libration invariant now covers 27 mass-conserved Enceladus stacks
plus a Titan-class stack at rel=1e-12 (exploratory maximum residual 1.64e-13).
Focused plus neighboring libration tests: 13 passed, 1 intentional
strict-xfail. No production physics changed. Implementation `093dc971`; C23
next.)
Prior entries follow.

Updated: 2026-08-15b (Machine A, three agents: **B2' REPAIR IMPLEMENTED
AND NUMERICALLY VERIFIED** — the Delta_rho-consistent Eq.-12 treatment
replaces the rescinded option-A path; every adjudicated number
reproduces (C2 bracket 25.99-27.34 km at Park values; surface-only ==
C2=0 bit-for-bit; default path byte-identical to HEAD at 6 structures);
pinning test moved in-commit; config records the retaken decision;
REVIEWER SIGN-OFF PENDING (dispatched) before B2' closes. COMMIT
HYGIENE NOTE: the repair's code+tests were swept into 5e5b9d78 (the
design-ruling commit) via concurrently-staged files — content verified
correct; message does not describe it; recorded here. **Sigma_model
wiring (B15) + config-schema reconciliation DONE** (22 new tests;
no-op bit-identical for existing campaigns; from_json now loads the
candidate). **FROZEN-BRANCH DESIGN RULING** (reviewer): config's frozen
branch was UNREACHABLE as written (F1: PP's argmin-MoI seafloor
selection = undeclared conditioning collapsing the rho_rock prior);
redesign APPROVED conditional on A1 — constant-rho route with exact
mass closure, separate frozen arrays (v3.2 schema), analytic
zb<->rho_rock bijection, TWO-RUN evidence-based branch comparison (no
theta indicator), Cuncertainty ruling (prefer NO MoI gate; 0.015 if
unavoidable; MODERATE_4's ±0.001 instruction COUNTERMANDED; ocean
branch must DECLARE its argmin selection). A1 smoke PARTIAL: constant-
rho solve works (rho_sil found=actual=2385, M=1.0000 M_E, no core) —
open: no-liquid variant + exact zb binding + window bypass, folded into
A2-A5 implementation. Suite 79 passed + 1 xfail.)
Prior: 2026-08-15 (Machine A: **B2' USER RULING EXECUTED + CODE GAPS
CLOSED + CODEX QUEUE DRAINED TWICE.** User adopted the reviewer
recommendation (§0.23): option A rescinded, Delta_rho-consistent
treatment ruled in, the Delta_rho identity now a STANDING INVARIANT
(C22: committed test, max residual 1.6e-13 across 28 stacks). Codex
C20-C24 all accepted (methods snippets for 8 campaigns; constraints-doc
Park/k2/C18 corrections + GSW extrapolation test; invariant test;
Enceladus GUI placeholder; doc-doctor pass 3 = 11 PASS / 1 FINDING vs
2/10 at first audit — C25 queued for the one finding). Sonnet
implementation pass closed the Enceladus code_gaps (manager-verified,
14 tests + baseline 40+1xfail green): registry params, C30 dispatch,
isostatic_hm2019 forward-model dispatch, B4 zb cache builder with real
smoke build. Still open: B2' repair (B, reviewer-signed), frozen-branch
builder + MAJOR-1 (reviewer-blocked), Sigma_model wiring, config-schema
reconciliation at ratification. Titan HF snapshot rebuilt+verified (all
3 compositions + no-ocean) — user ships. RESUME_NOTE ledger updated.)
Prior: 2026-08-14 (Machine B: **B2' ADJUDICATED — reviewer verdict BLOCK. Remains BLOCKING, but the block lands on the SHIPPED treatment, not on the measurement.** Ran the published-answer test the resume note called for (H&M section 3.3 publishes a libration-only shell thickness 16-22 km, against Thomas 2016's 0.120+/-0.014 deg, NOT the Park 0.092 the config conditions on). **CRITICAL:** the shipped `H22_obs_m` path — the config's adopted "option A / the H&M convention" — is NOT A VALID TREATMENT. PP's layer form telescopes into H&M Eq. 19's Delta_rho sum, `Ks/(3 omega^2) = rho_ice*f_top + (rho_ocean-rho_ice)*f_base` (verified ~1e-14), so the shell-base interface's physical weight is Delta_rho ~80-95 kg/m3 but appears in code as a difference of two ~10x larger cancelling terms. Scaling `(Bs-As)` scales only one half, applying an implicit, structure-dependent, SIGN-FLIPPING effective scale of +0.33 to -0.58 to the base figure while the docstring asserts interior interfaces stay hydrostatic — a code/docstring divergence on the SOLE thickness channel. Defensible treatments: hydrostatic and surface-only, Delta_rho-consistent Eq.-12 bracketed by C2 in 25.99-27.34 km; span 24.42-27.68 km over the mass-admissible D_ocean. **TWO prior headline numbers are ARTIFACTS of the same mis-weighting:** Eq.-12's "30.5 sigma / no zb matches" -> **+1.66 sigma at zb 25.99 km**, and the C2 sweep "0.1005->0.1824 deg" -> **1.6 sigma span, ~20x smaller**. So MAJOR-4's "gravity = compensation state, libration = sole thickness channel" split LARGELY SURVIVES; my earlier claim it needed revisiting is withdrawn. **HEADLINE, RATIFIED:** the H&M-vs-Park separation is essentially ENTIRELY the libration measurement revision (Thomas 0.120 -> Park 0.092) propagating through an unchanged forward model — **Delta = +5.6 km** (19.5-19.7 -> 25.1-25.4 km at H&M's own published core parameters, inside their JOINT Table 2 band and Park's band respectively) against a ~5.5 km band-centre separation. The two published bands are NOT competing answers; the campaign deliverable needs restating (proposed framing in the adjudication file). Frozen-branch ruling UNAFFECTED. Also: PP's 2x2 solve is ALGEBRAICALLY IDENTICAL to H&M Eq. 20 (rel <=2.2e-16) so the "not term-for-term identical" caveat was over-cautious about the form, though the COEFFICIENTS stay uncertified (K_int 8pi/15; Bsp_Asp linearized 1.7-2.0% low). The published-answer test did NOT discriminate (all treatments fit inside a band whose 1-sigma half-width ~2 km exceeds the whole treatment spread). Pole screen: right answer, wrong mechanism — the rejected zb~3.4-5.1 km are det(A_mat)=0 resonances and min|gamma| pre-pole is 0.40-0.48 deg, 3-5x above both targets, so no real branch was discarded; `det(A_mat)>0` is the correct invariant. **THIRD self-correction:** I claimed whole-difference and surface-only were the same code path — WRONG, `(Bs-As)` is a difference of consecutive ellipsoids so the base piece is -128% of the total; the prior three-distinct-treatments record was right. **USER DECISION NEEDED:** the 2026-08-13 "adopt option A" decision rested on two now-falsified premises; not reversed unilaterally, and the repair moves a test pinning the suspect +1.40-1.43 sigma shift. Campaign-level risk flagged: the recurring failure mode is perturbing a physical quantity through a DIFFERENCE rather than the interface weight that carries it — one root cause produced all three defects plus my own retracted claim; reviewer recommends the Delta_rho identity as a STANDING INVARIANT. Test baseline 11 passed + 1 xfail (K_int tripwire still xfail). No production code changed. Resume note: validation_reports/enceladus_isostasy/RESUME_NOTE.md; adjudication: validation_reports/enceladus_isostasy/b2prime_ADJUDICATED_drho_weighting.json. Prior entry below.)
Prior: 2026-08-14a (Machine A Codex: **C21 VERIFIED; DELEGATE QUEUE
COMPLETE.** Enceladus constraints now use the Park 2024 erratum gravity and
libration values, reject hypothetical k2, and record the C18 rigid-branch
adjudication (+0.8% correction, ±0.4% residual; 0.03% claim struck). C17's
standard-run C20/C22 output closes STRATEGY core-parity exception (b). The new
GSW extrapolation regression confirms finite, salinity-monotone rho/Cp/sigma
at 42/70/100 ppt, 3 MPa, 271–274 K; focused tests pass 6/6. Titan methods HTML
regenerates byte-identically from Markdown with its source stamp. Markdown,
reference, `py_compile`, and format checks pass. Implementation `a6a55631`.)
Prior: 2026-08-13a (Machine A Codex: **C20 VERIFIED.** A source-backed
methods-snippet generator now covers all eight requested campaigns (NH3,
MgSO4, NaCl, Europa v5/v6/v1.1/v4, Titan Test50). Each committed Markdown
output carries per-sentence source pointers and reproduces config priors and
observables, cache axes/flags/hashes, artifact training metadata, gate records,
and deployment state. The verifier regenerates byte-for-byte, checks
cross-source hashes/names, and renders all eight outputs through
Python-Markdown. Implementation `0fbbe6fe`; C21 next.)
Prior: 2026-08-12l (Machine A Codex: **C19 VERIFIED UNDER THE ADJUDICATED
RESPONSE-BAND DESIGN.** Enceladus now carries the 15.5592 hr PPO period for
response-only evaluation; no fixed PPO `Be1xyz` vector was invented. The
existing true-anomaly row resolves unchanged from arbitrary working
directories and cannot be silently reused as the synodic row. The inference
display reports `|Ae| x [1,2] nT` with the drifting PPO/no-stable-phase
caption. The Saur conductivity bracket, two-period finite-Ae smoke, loader
guard, band unit test, and AppTest all pass; focused induction/inference
regressions 52/52 pass. Implementation `eab1d788`; C20 next.)
Prior: 2026-08-12k (Machine A: **REBALANCING DIRECTIVE §0.22 — user
raised bog-down concern; manager agrees on the margin.** Salt C16
FAIL-track OPENED (reviewer: flow-dominated, 4.4/5.6 sigma) and PARKED
post-1.0 — reference ocean fractions are the citable numbers (MgSO4
~0.49, NaCl ~0.53: salt Titan models genuinely ocean-ambiguous, vs NH3
93%). NH3 lift gates -> background. ENCELADUS = foreground: Machine A
finishing freeze blockers (zb builder + libration scans on the
+0.8%-corrected model), then B runs cache->dataset->n_eff=8000
reference->training. Track 2b/2c post-1.0. B's repool-repair completion
accepted. Codex C18 verified under the revised spec (5 passed + 1
strict-xfail tripwire; Love-number-sum guard permanent).)
Prior: 2026-08-12j (Machine A Codex: **C18 VERIFIED UNDER REVISED SPEC.**
First libration regression suite now covers the VH16 published band
(524.38 m inside 466–590 m), a strictly decreasing 5–45 km shell sweep, the
adjudicated rigid/elastic discrepancy pin (0.110905164/0.111814157 deg), a
strict-xfail physical-direction tripwire for the future K_int repair, and a
permanent y1 partial-Love-number sum guard (0.0148% vs ALMA Re(k2), inside
0.1%). Focused + nearest Enceladus gate: 5 passed, 1 xfailed. No production
physics changed. Implementation `54bf8b4a`; C19 resumed next. Prior entry
below.)
Prior: 2026-08-12i (Machine A: **BOTH CODEX ESCALATIONS ADJUDICATED.**
(1) C18/libration: reviewer decomposed the 0.82% rigid-vs-elastic gap
exactly — the bug is the RIGID branch's K_int missing 8π/15 (+1.9%) plus
linearized figure moments (+1.6%) cancelling the REAL elastic reduction
(−2.6%; elastic branch validated end-to-end, Love-number sum matches
ALMA k2 at 0.015%, VH16 direction confirmed). Shipped rigid channel is
+0.75% LOW (0.21σ_obs at the Park solution) BY CANCELLATION — not
robust to rheology. Adopted: +0.8% multiplicative correction in the
frozen config + ±0.4% residual (freeze-doc B2); falsified 0.03% record
corrected in the 2026-07-09 handoff; K_int fix = post-campaign reviewed
task, must not land alone; C18 UNHELD with revised spec (published-value
+ monotonicity land; discrepancy pin + strict-xfail + permanent
Love-number-sum guard). (2) C19/induction: the PPO driver is
non-stationary — no stable Be vector exists; synodic channel represented
by response function + amplitude BAND A×[1,2] nT, no Be1xyz row
invented; orbital row stays the primary. C19 resumed under the amended
spec. Campaign remains rigid=True.)
Prior: 2026-08-12h (Machine A Codex: **C19 NOT IMPLEMENTED — SOURCE-DATA
ESCALATION.** Saur et al. 2024 §2.2 supplies the nominal PPO period
(0.6483 d = 15.5592 hr) but only a time-varying amplitude of up to 1–2 nT;
it does not supply the exact J2000 complex Cartesian vector required for a
`Be1xyz` row. Choosing its direction, phase, or a value within that range
would invent a scientific assumption, so Codex stopped at B8 without
production/test changes; B9/B10 remain not implemented because the required
two-period smoke cannot be completed. Claim `ee4b1293`; manager scientific-
data adjudication is required. Prior entry below.)
Prior: 2026-08-12g (Machine A Codex: **C18 NOT IMPLEMENTED — SCIENCE
ESCALATION.** The published Van Hoolst 2016 rigid model reproduces (524.38 m
inside 466–590 m) and the 5–45 km shell sweep is strictly monotone, but the
same-node DP-ALMA y1 comparison is rigid 0.110905 deg vs elastic 0.111814 deg
(0.8196%, failing the <0.1% criterion); a fresh full run independently gives
0.9452%. Both have the opposite sign from the paper's stated elastic
reduction. No physics or test files changed; manager/scientific adjudication
is required before C18 can resume. Prior entry below.)
Prior: 2026-08-12f (Machine A: **CODEX C17 VERIFIED.** Standard-run
degree-2 gravity is now available behind default-off
`Planet.Do.CALC_C20_C22`, returning direct-Clairaut `Gravity.kf/C20/C22`
with Europa/Titan inference-parity conventions and normal log printout.
Focused gravity verification is 20/20 PASS; full pytest is 162 PASS + 2
pre-existing Reaktoro-speciation failures (`SupcrtAqueousLookupByFormula`).
Implementation `da01ce26`; C18 next. Prior entry below.)
Prior: 2026-08-12e (Machine A: **C16 STOP RELEASED — re-ratified after
independent gate recomputation** (residual 1.20 sigma inside the
preregistered 2x bound; R4 recorded: the n_eff=2000 reference was the
biased estimator, LOW on branch mass, FAIL-ADJUDICATED with narrow
scope). §0.21: n_eff=8000 pooled set is the authoritative NH3 reference;
NH3 C12 sweep + corrected SBC UNBLOCKED; salt C16 stays evidence-only
until R3-class reference recomputes per composition. The corrected NH3
ocean fraction ~0.932 is now ratified. Codex C15 accepted earlier;
C16t in progress.)
Prior: 2026-08-12d (Machine B: **§0.20 R1→R3 COMPLETE — C16 RE-RATIFIED
BY REVIEWER (PASS-WITH-CONCERNS); STOP RELEASE IS MACHINE A's CALL.** The
C16 tension is CONSTRUCTIVELY RESOLVED: the n_eff=2000 reference was biased
LOW by sampling resolution, not the corrected side. R1 (PASS-w/-concerns):
neither side weighting-fragile; reference-side seed scatter (span 0.0121 at
n_eff=2000) named as suspect. R2 (PASS-w/-concerns, N=100k): corrected side is
a precisely-resolved fixed point (~0.933) that moved AWAY from the reference at
higher N — falsifies finite-N corrected bias. R3 (DECISIVE, n_eff=8000, 3 seeds
72/172/272, all R-hat=1.000, ESS~12.9k): reference rose +0.0114
(0.91725→0.92865), residual collapsed +0.0149(3.6σ FAIL)→+0.00352(1.20σ),
inside the committed 2×combined-SE bound (0.00588). |Im_k2| between-seed std
collapsed to 0.00015. **Scientific-reviewer (a8d0ea37) INDEPENDENTLY reproduced
the gate to the last digit, confirmed the reference move is like-for-like
(identical config_hash + cache sha256, regime-preserving n_active ratio 2),
verified seed-272 inclusion is CONSERVATIVE (dropping it shrinks residual to
+0.00149), and ruled PASS-WITH-CONCERNS re-ratifying C16 + supporting STOP
release.** All 4 reviewer required-validations (non-blocking) folded into
R3_decisive_reference_recompute.md; provenance seeds_note added to
matched_reference_report.json. Artifacts:
validation_reports/nh3_diagnosis/{R1,R2,R3}*.md +
matched_reference_neff8000/. **PENDING MACHINE A: release the MANAGER-GATE STOP
on NH3 Track-1 corrected compute + clear the deployed corrected NH3 ocean
fraction (~0.932) as ratified — flip `pass`/clear STOP in
is_validation_nh3.json.** R4 RECORDED (biased estimator = n_eff=2000 reference,
biased LOW; FAIL-ADJUDICATED) as `r4_resolution_record` in the C16 block WITHOUT
flipping `pass` (that flip is Machine A's). §0.20 R1→R4 sequence complete; only
the manager status-flip + STOP release remain. Not pushed — Machine A owns
pushes. Prior entry below.)
Prior: 2026-08-12c (Machine A: **C16 RESOLUTION AUTHORIZED (§0.20)**
answering B's ⛔ escalation — Phase 1 record accepted (correctly-fired
preregistered trigger; STOP was right). Order: R1 zero-compute
reference-weighting + PSIS/ablation audits from saved artifacts; R2 one
N=100k corrected run (persistence prediction preregistered); R3 decisive
n_eff>=8000 3-seed NH3 reference recompute, re-ratify only on
2x-combined-SE agreement; R4 biased estimator recorded
FAIL-ADJUDICATED. Salt fiducial validations authorized in PARALLEL with
C16 as evidence-not-gate + the same-sign-residual reading preregistered;
Track 2c authorized when cores free. v5 near-no-op control + CMI 0.605
nat accepted — the flow's Im_k2 under-update is confirmed RECOVERABLE
information, not a data limit. Non-blocking: methods .html should be
generated from the .md per the docs ruling. Codex C13 accepted
separately. Codex C14 verified: salt RATIFICATION consolidations +
future-run gate-manifest provenance schema, implementation `f282b669`.
Codex C15 verified: SBI methodology + GUI capability docs, implementation
`87f4ca31`. Codex C16t verified: 13-slot artifact/config/cache provenance +
fixed-seed Test50 reproduction regression, 15 tests passing in 2.20s,
implementation `9754646e`; C17 next.)
Prior: 2026-08-12b (Machine A: **ISOSTASY MODULE IMPLEMENTED — B13
REPRODUCTION GATE PASSES ON FIRST RUN.** PlanetProfile/Gravity/isostasy.py
(core-parity per user directive: available to CLI + App): H&M equal-
pressures Airy + finite-amplitude + interface gravity, unit tests 12/12
against reviewer-computed brackets; the B13 gate reproduces H&M's
published inversion from their inputs — misfit minimum shell 20 km,
ocean 39 km, core 193 km, rho_core 2373, ALL inside their box (19-24 /
35-39 / 192-195 / 2340-2410), L=1.91. Shape-channel fit path VALIDATED.
Remaining before config freeze: B1'/B2'/B5 scans on the assembled model,
B6 rheology survival, B7 deltaT check, reviewer consolidation (frozen-
branch support model + final sketch), zb cache builder (B4). §0.19
punch list posted for Machine B. Prior entry below.)
Prior: 2026-08-12 (Machine A: **ENCELADUS FIT PATH CONFIRMED — freeze
design complete pending consolidation.** User rulings: k2 omitted;
no-ocean cases STAY (overrides reviewer); zb axis + Seawater ratified;
"include libration as a parameter" interpreted as conditioned observable
+ per-sample derived display. Shape directive resolved: H&M 2019 fit
shape as forward-model INPUT, not observable — Airy equal-pressures
isostasy coupling adopted (dC20/dC22_nh free boxes → one sampled
compensation fraction; C30 ADDED as observable; Tajeddine 2017 shape
primary, Nimmo ablation; gravity becomes ~0.5σ/km shell channel).
New blockers B11-B16 incl. the H&M REPRODUCTION GATE (module must
reproduce shell 19-24/ocean 35-39/core 192-195 km from H&M inputs;
fallback = display-only + restored boxes). Implementation spec written
(plans/active/enceladus-isostasy-module-spec.md); B16 verified
(compute_eccentricities is full iterative Tricarico, per-interface).
Codex C13 verified (artifact/deploy/plan-index currency; implementation
`411167fb`); C14 next, then C15-C19 in order. C18 is the libration regression
test (Librations.py has zero coverage); C19 is induction plumbing B8-B10.
Frozen-branch-vs-Airy tension
resolution + final consolidated config sketch pending the reviewer
consolidation pass. Machine B §0.18 Phase 1 results still pending.)
Prior: 2026-08-11d (Machine B: **⛔ MANAGER-GATE STOP — NH3 C16 REOPEN.
ESCALATION TO MACHINE A.** [v5 control + NH3 Track-2a CMI since folded in below.] §0.18 P1.3 NH3 fiducial IS validation ran clean
on C3 (byte-exact, max_rel_diff 0.0), C5.3 (0.00094<0.01), C10 pushforward
(all gaps ≪ threshold; Im_k2 median 0.0449→0.109 recovers MCMC ceiling
~0.103) and the C13 pushforward-median stability gate (0.057/0.073 σ_obs,
tidal deliverable seed-stable — PASS). BUT C13 3-seed evidence
(offsets 0/1000/2000, N=20k) shows the corrected ocean fraction
0.9293/0.9310/0.9363 vs reference 0.9173 = residual +0.012/+0.014/+0.019:
**stable, all-positive, NOT shrinking.** scientific-reviewer
(agent a4ebcd7bd368dc8d6) OVERTURNED the earlier fiducial-only C16 PASS →
**BLOCK / C16 REOPEN**: the residual +0.0149 is ~19× the 1/ESS finite-N
bias ceiling (7.9e-4 at ESS~1270, Pareto-k 0.19 clean) so it is provably
structural, not finite-N — the validation-3 "shrinks with N" premise is
FALSIFIED. Between-seed std 0.00367 matches the SNIS delta SE and falsifies
the crude ESS SE; under the validated SE the residual is 2.2–3.4σ/run,
~2.3–3.5σ pooled = the preregistered reopen trigger. NOT the plan-149
catastrophic mode (no-ocean branch preserved 0.064–0.071, C5.3 passes) —
a ~1.5% branch-mass inconsistency between two same-target (C3 byte-identical)
estimators; escalation points at the **reference** side (n_eff=2000 modest;
pooled-ref weights repaired only in 0acff866). **MANAGER GATE: STOP all NH3
Track-1 corrected compute; Machine A to authorize resolution** (larger-N
corrected run N~100k–200k + higher-n_eff reference recompute n_eff≥8000 +
identify the biased estimator; re-ratify only on agreement within combined
SE). Reviewer CLEARED as independent (may proceed, must NOT consume the NH3
corrected ocean-fraction): Europa v5 positive control + NH3 Track-2a CMI —
**both now DONE.** (1) **Europa v5 IS positive control: clean near-no-op**
(verdict clean, ESS 13487 / ESS/N 0.674, Pareto-k 0.256, w_max 0.0013;
Im_k2 0.00264→0.00263 Δ2e-6; C3 byte-exact) — the IS machinery does NOT
distort a well-calibrated flow, corroborating that the NH3 C16 tension
(ESS/N 0.06) is reference-side, not a correction-machinery defect
[e78e21f9]. (2) **§0.18 P1.4 Track-2a NH3 CMI (GATE) = 0.605 nat**
(band [0.599,0.615], ~6× the 0.1-nat threshold; NH3 npz regenerated seed 72,
n_kept 689845 == train manifest) → **Track 2b/2c NOT cancelled**; NH3 lands
with the salts (NaCl 0.63, MgSO4 0.39) [bb791648]. Commits 5cc6650a
(fiducial + C16 PASS, later overturned), c2178927 (C13 evidence),
1c289854 (C16 REOPEN/STOP + escalation), e78e21f9 (v5 control),
bb791648 (NH3 CMI). **§0.18 Phase 1 now BLOCKED solely on Machine A's C16
resolution authorization** (all other P1 items — P1.1/P1.2/P1.3-cleared-parts/
P1.4 — discharged).)
Prior: 2026-08-11b (Machine A: **PROGRAM REVIEW + DOC AUDIT.** Opus
reviewer adjudicated plans/next steps: Enceladus starts NOW in parallel
(hold point moved to FLOW TRAINING with a preregistered release rule;
config freeze must evaluate 2D Tb×w — Enceladus has the program's
strongest composition data); Track 2 MORE urgent post-C1 (success metric
is ESS/N — the correction is a fiducial bridge, NOT viable for arbitrary
user-x at 1-3 CPU-h/click on HF); v5 IS control PROMOTED to Phase 1
(only preregistered positive control); 0.9326 ocean-fraction plumbing
check BLOCKING before C16 is read; corrected-SBC budget corrected to
~39 CPU-h/comp at deploy-N → NH3-only first (~80 CPU-h saved on
conditional extension); v7 recommended RETIRED-as-ablation (run B4,
skip B5 — user to ratify); **1.0 MILESTONE declared** = mission coverage
(Enceladus the only missing pair; canonical-model naming recommended v6
for Clipper–Europa, user decision; quarantine lift NOT a 1.0
prerequisite). Doc audit: two numeric slips fixed (MgSO4 containment
0.004σ not 0.007σ in INDEX + slot); salt RATIFICATION.md consolidation
identified as highest-value gap; C13–C17 queued for Codex (INDEX
deployment-gate rewrite + table promotion, DEPLOYING 8-cache update,
plans-index refresh/archiving, salt RATIFICATION compilation +
manifest schema, methodology doc IS-correction section + GUI README,
slot-reproduction regression test, C20/C22 standard-run output).
§0.18 B reorder written; STRATEGY production table + 1.0 section
updated.)
Prior: 2026-08-11 (Machine A: **PATH TO LIFTING THE SPLIT RATIFICATION
FOUND AND CONFIRMED.** User directive: amortized flows must be trustable as
representative of the MCMC. Remedy plan authored
(plans/active/tidal-sector-remedy-plan.md): Track 1 = importance-sampling
correction (flow as proposal, exact MCMC likelihood for weights — both
ingredients already computed by the deploy path). Reviewer adversarial
adjudication: APPROVE WITH CONDITIONS C1–C16 (Pareto-k primary gate,
byte-identity asserts, ocean-fraction gate, SBC of the corrected pipeline,
cost-architecture choice; several factual corrections to the draft folded
in). **C1 stop-gate EXECUTED and PASSED on NH3**: corrected Im k2
pushforward median 0.1064 vs MCMC matched ceiling 0.1037 (flow alone
0.044) — the flow defect CLOSES under correction; Pareto-k −0.18 (clean),
w_max 0.009, ESS 680@10k (N≈15–20k for the 1000 floor; 2.3 min/10k on A).
NEXT: Machine A implements Track 1 under C2–C16 (validation driver +
deploy path), Machine B validates per composition + runs Track 2
(information-gain comparison first, {Re,Im k2}-only conditioner ablation,
transform pilot at 3 seeds). Quarantine lifts per composition only after
the preregistered gates pass; the model-data k2 tension caveat SURVIVES
any lift.)
Prior: 2026-08-10e (Machine A: **MgSO4 + NaCl COUNTERSIGNED (split-status)
and GUI-WIRED — Titan production line now three compositions deep.** Manager
independently verified the committed gate numbers (SBC 13/13 both; MgSO4
crosscheck 13/13; NaCl eta_V FAIL reproduced from the report, 0.352 vs 0.30
dex, recorded FAIL-ADJUDICATED-ACCEPTABLE; pushforward medians + 0.24σ/0.53σ
median-to-median gaps confirmed to the digit). Both placeholder slots
replaced with full slots: split-status sector warnings carry the two
reviewer-mandated caveat-copy corrections (median-to-median gap quoted;
deployed SBI Re_k2 marginal stated CONSERVATIVE vs reference MCMC), MgSO4
high-w extrapolation caveat + NaCl eta(ice V) do-not-cite included;
ocean_comp threaded for wedge labeling. INDEX rows added, countersigned.
AppTest 10/10 (new salt-slot test + NaCl-is-newest-default fix). Deploy
snapshot rebuilt: 8 registry caches exact, 317M — USER ships HF to take
v5/v6 + both Titan salts + globe fallback + provenance exports live in one
deploy. Task #52 Phase C complete.)
Prior: 2026-08-10d (Machine B: **reviewer req-val #1 + #2 DISCHARGED.** Gate
interpretation verdict [agent a719438103ced11d0]: **PASS WITH CONCERNS, no STOP**
— gravity C20/C22 deploy TRUSTED, tidal k2 quarantine STANDS, NaCl eta_V FAIL
non-blocking (poorly-identified HP-ice-V nuisance, threshold-adjacent to MgSO4's
passing 0.209 dex), repool repair statistically sound (Kish ESS + n reproduce
crosscheck reports to the digit). Two required validations before Re_k2 deploys:
(#1) **MCMC→Re_k2 pushforward RUN** (titanG_mcmc_rek2_pushforward.py; pooled
repaired MCMC posterior systematically resampled unweighted → identical θ→x
forward loop, obs_noise=False, byte-identical to the SBI ppc). RESULT: both comps
**UNDER-PREDICT** Re_k2 — MgSO4 pp_med 0.570 (−0.79σ), NaCl 0.575 (−0.68σ) vs
obs 0.608; NEITHER centers on 0.608, so the flow-offset/quarantine trigger is NOT
met. Gravity pushforward centers exactly on obs (C20/C22 ≤0.05σ); Im_k2 also
under-predicts (−1.60/−1.34σ, consistent with the quarantined NH3 mechanism).
Median-to-median gap (the apples-to-apples statistic the decision rule uses):
MgSO4 0.24σ, NaCl 0.53σ. (#2) **limits monotone doc corrected** — recorded as
gated FAIL + adjudicator override, not N/A (below). **Pushforward reviewer verdict
[agent ad7ae4436e53cdc1e]: PASS — BOTH MgSO4 and NaCl may deploy under full
split-status (C20/C22 TRUSTED, Re_k2 informative-with-caveat, tidal k2
QUARANTINED). NO remaining STOP.** NaCl "INTERMEDIATE" confirmed a mechanical
threshold artifact (SBI more pessimistic than MCMC, both under-predict in the
same direction) — SAME model-data-tension verdict as MgSO4. Reviewer caught a
reporting defect: the reports mix deviation-of-median (MCMC) with
median-of-abs-deviation (SBI ppc, spread-inflated); apples-to-apples
median-to-median gap is 0.24σ/0.53σ, and the DECISION-RULE gap statistic already
used the correct median-to-median comparison so the verdict is unaffected. Two
GUI-caveat-copy prerequisites for Machine A (NOT blockers to physics): (1) quote
the median-to-median gap (0.24σ/0.53σ), not the mixed statistic; (2) state the
deployed SBI Re_k2 marginal is CONSERVATIVE relative to the reference MCMC
(tension-leaning, not a tightened bound — safe-side for informative-with-caveat).
NEXT: **campaign compute COMPLETE.** Phase C GUI slot wiring (task #52) + the two
caveat-copy corrections are a Machine A duty. commits 99ed8c42 (+ this).)
Prior: 2026-08-10c (Machine B: **§0.16 gate seq COMPLETE — SBC/limits/crosscheck
run for both MgSO4+NaCl; reviewer routing for interpretation.** Pushforward-
acceptance verdict [agent ab702bb410edb3da6]: **PASS WITH CONCERNS, PROCEED no
STOP** — tidal k2 quarantine correctly STANDS both comps; Im_k2 under-update
reproduces the established NH3 flow-defect+data-limit mechanism (no island/EOS
defect; Margules region carries no anomalous Im_k2); gravity+Re_k2 deployable
(Re_k2 = informative-with-caveat, NaCl 1.52σ). Gate outcomes (interpreted, never
tuned): **SBC (n_sbc=1500) PASS both** (BH-FDR corrected; gravity/mass-cons ranks
uniform ⇒ reviewer req-val #2 satisfied). **limits FAIL both but NON-BLOCKING** —
containment_pass=True (binding read, max shift ≤0.007σ ≪ 0.25σ). Correction
(reviewer 2026-08-10): the FAIL is the monotone clause, and it is a *gated FAIL*
(monotone_pass=False, monotone_gated=True, monotone_verdict=FAIL,
monotone_na_reason=None), NOT the code returning N/A — the sweep window extends
to 0.2/0.3 (above the LIMITS_MONOTONE_FALSIFIED_BELOW=0.15 auto-N/A trigger), so
the gate DID evaluate monotonicity and it failed (medians INCREASE with Im_k2:
MgSO4 12.61→12.86, NaCl 12.81→13.18 — opposite the decreasing premise). It is
rendered NON-BINDING by adjudicator OVERRIDE on the established falsified-premise
basis (Titan-joint monotone-decreasing Im_k2→η_Ih premise invalid; 2026-07-09
MCMC ground truth) + the fact that the deployment x_obs Im_k2=0.135 sits BELOW
the 0.15 falsified boundary (real obs is in the auto-N/A region; the sweep only
probes above 0.15 for limiting behavior). Recorded as adjudicated override, not
as monotone_pass=None.
**crosscheck: MgSO4 PASS (all 13 params); NaCl FAIL on log10_eta_V only** (shape+
median, median_diff 0.352 vs 0.30 dex tol) — poorly-constrained HP-ice-V nuisance;
all observable-relevant params (Tb,w,dC20,dC22) + primary eta_Ih PASS. Fixed a
pooled-reference weights bug (pooling never updated res.weights → crosscheck
length mismatch; repaired both samples+weights renormalized 1/n_seeds, no MCMC
recompute; titanG_repool_reference.py). commit 0acff866. NEXT: reviewer interprets
the full gate set (esp. req-val #1 Re_k2 model-data tension + NaCl eta_V FAIL) →
split-status deploy decision (Machine A GUI wiring). Reviewer routing dispatched.)
Prior: 2026-08-10b (Machine B: **§0.16 EXECUTED — all four jobs launched +
condition (a) CLEARED.** Acted on the FLOW TRAINING GO "everything in parallel":
NaCl+MgSO4 flow training (DEPLOYED nsf, seeds 74/73) and NaCl+MgSO4 B3 reference
MCMC (n_eff=2000/n_active=1024, 3 seeds each) all running — cache SHAs verified
against gen manifests, thread-pinned separate processes. Scientific-reviewer
[agent a909e4b1d07ba13ac] adjudicated the #4 corrected multi-seed report:
**PASS — does NOT overturn the eliminated-capacity verdict; MgSO4/NaCl training
CLEARED to continue on the deployed architecture** (max arm×seed 0.0471 ≪ 0.0862
bar; all concentration_ratio > 1; A/C converged FEWER epochs than D0, hit_ceiling
False — opposite of premature early-stopping). §0.16(a) binding gate SATISFIED.
Reviewer MODERATE provenance fix + MINOR notes applied to the manifest (commit
1a993b56). Reviewer archival follow-ups non-blocking. NEXT: as each training
completes, run its post-training pushforward + tidal-quarantine gate (§0.16 gate
seq step 1); ref MCMCs are the long pole (~hours). commits 1fb0af97/3c5e1cc9/1a993b56.)
Prior: 2026-08-10 (Machine A: **FLOW TRAINING conditional GO** for
MgSO4/NaCl, §0.16 — the per-composition quarantine re-verification is
restructured from a pre-training hold into the mandatory first
post-training gate (it requires a trained flow; NH3 pushforward method;
split-status default), and the #4 corrected-re-run report becomes a
commit-with-first-deliverable condition (stop+escalate if it overturns
the eliminated-capacity verdict). B also starts fresh per-comp reference
MCMC (B3 protocol) immediately in parallel. Codex C9–C12 all accepted +
pushed earlier; B has not pushed since 08-07 — likely awaiting this GO.)
Prior: 2026-08-07b (Machine A: §0.15 closeout `verified` — rebuilt
MgSO4/NaCl cache bytes RE-ACCEPTED (independent reload: counts exact, SHAs
match gen manifests, 0 sandwiched frozen nodes, 252 K repairs confirmed);
`ResetNearestExtrap` int-truncation fix adjudicated + landed with
regression tests (4/4; consumer smokes 7/7). GUI: WebGL-free static globe
fallback + SVG sample picker shipped after user report ('WebGL not
supported'; Safari worked — browser-side WebGL, app now browser-proof);
AppTest 9/9, static PNG visually verified. Deploy snapshot rebuilt —
USER ships HF.)
Prior: 2026-08-07 (Machine B: **GATE-3 GO from Machine A → re-run + reviewer
re-adjudication complete.** Gate-3 re-run fresh on the installed REBUILT
(float-coerced) MgSO4/NaCl caches; scientific-reviewer [agent a02a9038b2de1d382]
independently reconstructed both grids/onsets/mass-balance/252-K-repair from the
bytes → **BOTH PASS, no blocking issue for flow training.** NaCl half_cell 0.500 K
AT the 0.5 K bar (uniform 1.0 K Tb; non-strict ≤ = PASS; ocean-frac diag span
0.0218); MgSO4 half_cell 0.250 K PASS. MgSO4 w≥150 corner diag FAIL (span 0.05)
CONFIRMED benign (coarse-denominator, demoted metric). MgSO4 zero-None confirmed
correct (buildable box + typed-exception retry, NOT masked failure); NaCl 66 None
= genuine high-w/high-Tb corner. §0.15 updated. Flow TRAINING still HELD on #4
architecture-pilot verdict. Machine A: re-verify new cache bytes + latent
ResetNearestExtrap dtype bug (non-blocking for these caches per reviewer).
Prior: 2026-08-07 (Machine A manager session: **v5 RATIFIED (scoped) + v6
RATIFIED** — final adjudication appended to
plans/active/europa-v5v6v7-gate-adjudication.md; dC22_nh recorded
FAIL-ADJUDICATED-ACCEPTABLE (v5-local nuisance undertraining; v6 owns that
deliverable, clean); 2nd-seed SBC queued non-blocking. GUI: not-ratified
gating LIFTED for the six v5/v6 slots, ratified scope notes + always-visible
gate-status captions; AppTest verified, app suite 8/8. Deploy snapshot
rebuilt at 536ddfe2 (299M, registry cache set exact) — USER ships HF
(carries v5/v6 + Ae=0.9 synodic defaults + e-notation fix). MgSO4/NaCl
production caches ACCEPTED after independent pkl verification (§0.14:
node counts exact, deepest MgSO4 liquid 1371.1<1400 MPa, extrap column
monotone) — per-comp gate-3 + 1M dataset gens GO; flow training still HELD
on per-composition quarantine re-verification. C8 doc-doctor triaged:
C9–C12 queued for Codex (ledger/routing, assumption-text contradictions,
scope-note priors, figure provenance); v4 direct-Clairaut recorded as
STRATEGY core-parity exception (b); Geotherm-tab deviation recorded as
accepted. INDEX.md v5/v6 rows updated to ratified.)
Prior: 2026-08-06 (Machine B: §0.10 priority-(a) v5/v6 B1/B2/B6 gates EXECUTED
vs the FRESH pooled n_eff~2000 v5 reference, reviewer-adjudicated.
**v5 = PASS WITH CONCERNS** — ships for its DEPLOYED D_iceIh/ocean/salinity
deliverable with a mandatory scope-note (Machine A finalizes + countersigns):
SBC FAILs on dC22_nh ONLY (BH-adj p=0.0185; 621 kept pairs — the B2 n→1500 fix
EXPOSED a miscalibration the old n=108 masked), a real but LOCALIZED
v5-flow-specific undertraining on the ~100:1-compressed nuisance (rank-CDF mild
monotone; NOT free-gravity-family — v6's identical nuisance PASSES p=0.61);
deployed marginals doubly validated (SBC-clean D_iceIh p=0.996 + crosscheck PASS
even without the 0.36 km floor). **v6 = ALL GATES PASS.** B6 DEFERRED (empty
anchor set, below the 0.15 falsified boundary). **#4 architecture-pilot EXECUTED — capacity/embedding
ELIMINATED** (reviewer-adjudicated): uniform negative across D0/A/B/C, more
capacity made concentration WORSE. The Im_k2 gap PARTITIONS — 0.042→~0.093 is a
real flow-defect (training-signal/identifiability, NOT representational) and
~0.093→0.135 is a data-limit (obs sits in the upper tail of the EXACT posterior;
reachable frac≥obs ~0.19, not 0.5). Reviewer: MgSO4/NaCl proceed on the DEPLOYED
architecture under split-status (re-verify quarantine per-composition); next NH3
tidal diagnostic is upstream identifiability, NOT a bigger flow. FOR MACHINE A:
finalize v5 ratification + countersign; deploy v6. See MACHINE-B-HANDOFF §0.11
(v5/v6), §0.12 (#4), §0.13 (MgSO4/NaCl build prereqs). MgSO4/NaCl cache builds:
reviewer SIGN OFF on NaCl extrap_ocean (no-op≡clamp) + SIGN OFF WITH CONDITIONS
on the per-composition Tb grid; empirical onset probes run — NaCl is a ~39 K
DISJOINT ocean/frozen diagonal (→30-40 fine Tb nodes, not NH3's 12), MgSO4
~15 K [240,255] K with a REAL (fine-Pfreeze-confirmed) non-monotone frozen
island at the w=194 (2-molal) cap. HELD on reviewer adjudication of whether that
low-Tb ocean island is physical (keep) or pathological (exclude) before sizing
the MgSO4 grid + running gate-3 (half-cell Tb-shift). #4 corrected 4-arm×3-seed
re-run in flight (recovers best_val/train_val_gap; A@72/172 ~0.043, no pass —
corroborates the eliminated-capacity verdict).
(§0.10 priority-b, architecture-independent). MgSO4/NaCl build gates now CLEARED
(reviewer a3c9ed8ae24664527, 2026-08-06): MgSO4 w=194 island = PATHOLOGICAL
(Margules melting-monotonicity violation → excluded by construction + post-build
has_ocean invariant); NaCl w=290 monotone + retry-corner discriminator PASS;
MgSO4 extrap_ocean=True adopted (physics verified monotone/stiffening to 1500
MPa), ceiling raised 1200→1400 MPa to preserve the w=194 salinity cap, with
stratified eos_extrapolated flags + dρ/dP band invariant + T-spot-check wired
into the NEW committed build driver plans/scripts/titanG_build_ocean_cache.py.
Both configs authored + BOTH PRODUCTION CACHES BUILT, VALIDATED, and PUSHED for
Machine A review: NaCl (40 Tb × 15 w = 600 nodes; 314 ocean/219 no-ocean/67 None;
tilted-diagonal has_ocean map matches onset probe) at
Test54_nacl_ocean/titan_nacl_joint_structure_grid_2d.pkl; MgSO4 (17 Tb × 16 w =
272 nodes; 191 ocean/81 no-ocean/0 None) at
Test54_mgso4_ocean/titan_mgso4_joint_structure_grid_2d.pkl. All MgSO4 reviewer
guards report CLEAN in-manifest: island invariant OK (0 has_ocean=True at
w≥180 & Tb<248), extrap-ceiling OK (deepest w=194 node 1371<1400 MPa; 0 rejects),
dρ/dP∈[0.14,0.17] in all extrap columns, ocean-base T≤300 K everywhere; 45 deep +
7 mild eos_extrapolated flags on the high-w tail (carry the extrapolation caveat
in any high-salinity claim). NEXT (Machine B, no user input): per-comp gate-3
(half-cell Tb-shift → ocean-fraction, MgSO4 in the w≥150 corner), then 1M dataset
gens (seeds NaCl 74/7474/74, MgSO4 73/7373/73). Flow TRAINING still HELD on
per-composition quarantine re-verification (§0.10). Prior NH3 #1 → PERSIST
(salinity ELIMINATED) below.)
Prior: 2026-08-06 at genai `7bef8ba2` (Machine A Codex: C8 first
doc-doctor pass `verified`; report `94ded63f` covers all 12 checklist items,
2 PASS / 10 FINDING, with no findings fixed. Codex queue empty; manager triage
next).
Prior: 2026-08-05 (Machine B: NH3 follow-ups #3, #2 AND #1 EXECUTED, reviewer
PASS on all. #1 salinity-fixed pilot + matched-N control → PERSIST: fixing
salinity at matched N did NOT recover the Im_k2 update (banded 0.0313 ≈ control
0.0321, +0.02σ_obs, wrong sign for salinity); size effect small/monotone/toward-obs
→ salinity axis ELIMINATED as the driver of gap (a); remaining candidate is
capacity/embedding (#4). STOPPED per protocol — remedy selection is manager +
reviewer + user. MgSO4/NaCl stay HELD. See below).
Earlier: 2026-08-05 at genai `5ae7d17c`. Manager judgement §0.9 issued:
B3 accepted — the 1.06 km v5–v7 reference disagreement was an n_eff=500
resolution artifact; both old references wandered (fresh truth ~61.4–61.8 km);
§0.7 crosschecks retarget the pooled fresh n_eff=2000 references with a
preregistered 0.36 km shape-excess floor, and the earlier "v7 flow 22/22 vs
v5 reference" claim is corrected (void — wandered target). NH3 separators
accepted: joint mixture and noise-swamping ruled out for the Im_k2 miss;
follow-ups reordered #3 (pushforward plot) → #2 (matched-resolution NH3
reference — the 0.100 MCMC-pp ceiling inherits the discredited resolution
class) → #1 only if a material gap survives; MgSO4/NaCl release criteria
preregistered on the #2 outcome (details MACHINE-B-HANDOFF §0.9).
**#3 EXECUTED 2026-08-05 (reviewer PASS):** pushforward-shape diagnostic on
BOTH the capped anchor AND the deployed 1M flow confirms concentration-failure
and EXCLUDES wrong-mode — deployed concentration ratio 1.215 (>1: posterior
not narrower than prior), high-k2 tail retained (frac≥ceiling 0.233), but
frac≥obs moved only 0.128→0.161 (ceiling-independent under-update evidence).
Nonfinite-drop bias closed (<0.25% of draws, conservative direction).
**#2 EXECUTED 2026-08-05 (reviewer PASS):** matched-resolution NH3 reference MCMC
(n_eff=2000/n_active=1024, seeds 72+172, config unmutated) re-measured the
weighted |Im_k2| ceiling = **0.1037** (n_eff=500 was 0.0999; +0.0038 = +0.11σ_obs
and *away* from SBI — the ceiling is NOT an n_eff=500 artifact, the decisive
contrast with B3). SBI-pp 0.0423 vs matched MCMC-pp 0.1037 → **gap +1.76σ_obs ≫
0.5σ → SURVIVES**. With #3 (concentration-failure, wrong-mode excluded) the
under-update is a genuine flow-training deficiency. **Decision: run #1
(salinity-fixed retrain) with reviewer+user sign-off; MgSO4/NaCl PROCEED
falsified, stay HELD.** Reviewer note: matched MCMC-pp itself sits 0.89σ below obs
(ordinary model/data tension, NOT a remediation target); #1 targets only the
1.76σ SBI-vs-MCMC gap. See `validation_reports/FLOW_UNDERUPDATE_DIAGNOSIS.md`
(Follow-ups #3 + #2) and `validation_reports/nh3_diagnosis/matched_reference/`.
**#1 EXECUTED 2026-08-05 (reviewer PASS): PERSIST → salinity ELIMINATED as the
driver.** Design-review reviewer (PASS WITH CONCERNS) required a sample-size-matched,
salinity-VARYING control (the fixed-salinity band keeps only ~9%/60,039 rows vs
~690k — an ~11× cut biases toward under-concentration on its own); folded in +
band re-documented as a symmetric ±0.084-dex slice about 12.6 ppt. Results (Im_k2
pp-median): banded(N=60k,fixed) 0.0313 ≈ control(N=60k,varying) 0.0321 → salinity
axis IN vs OUT = +0.02σ_obs, wrong sign for the salinity hypothesis; control(60k)
vs S1(643k varying) 0.0392 = −0.20σ (pure size effect, toward obs, discharges the
confound); banded−ceiling(0.1037) = +2.07σ still-open flow gap. Reviewer PASS:
salinity cleared for gap (a) [SBI-pp vs MCMC-pp ceiling], remaining candidate
capacity/embedding (#4); size gain strongly sublinear (11× data → +0.007; ceiling
needs +0.065) so more ocean-only data won't close it; banded pp 5–95 [0.0023,0.340]
spans past obs → expressiveness/concentration signature, not support failure.
STOPPED per protocol; #4 remedy selection is manager + reviewer + user.
`validation_reports/nh3_diagnosis/f1_salinity_fixed/`. Also
2026-08-05: e-notation display fix for C20/C22-scale parameters in the
corner plot + summary table (`12bdbf56`, verified). Earlier-entry history:
`git log -p -- plans/STATUS.md`.
Freshness rule: any session that pushes commits, integrates Machine B
artifacts, or changes a queue must refresh this file's `Updated:` line and the
affected sections in the same session. If this file is >7 days older than
`git log -1`, refresh it before relying on it.

## Roster and lanes (2026-08-01)

| Lane | Models | Role | Queue |
|---|---|---|---|
| Machine A — Claude Code | Fable 5 (manager) + Opus 5 / Sonnet 5 subagents | planning, scientific review/adjudication, integration, GUI, smokes, cross-machine coordination | this file |
| Machine A — Codex 5.6 | Codex 5.6 | delegated implementation + reconnaissance tasks | `plans/CODEX-QUEUE.md` (instructions: `AGENTS.md`) |
| Machine B — Claude Code | Opus 4.8 / Sonnet 4.6 / Haiku 4.5 | heavy compute: cache builds, reference MCMC, SBI training, gates | `plans/MACHINE-B-HANDOFF.md` |

Claude Fable 5 on Machine A is the model manager: it owns adjudication and
decides what enters the Codex and Machine B queues. Scientific verdicts are
never delegated to Codex.

## SBI slot states (GUI `_SBI_ARTIFACT_SLOTS`)

| Slot | State | Notes |
|---|---|---|
| Titan (Andrade, no-ocean) Test50 8D | DEPLOYED | unchanged |
| Titan freegrav no-ocean | delivered | artifact + gate reports in repo |
| Titan NH₃ joint no-ocean+ocean | **SPLIT ratification 2026-08-03: gravity/structure VERIFIED, tidal (k₂) sector NOT verified**; GUI slot wired (in-browser render pending) | `titan_freegrav_nh3_posterior_1m.pt`; SBC pass; PPC shows flow under-updates k₂ (SBI pushforward at prior-pred median vs MCMC/data) — MCMC reference authoritative for k₂/dissipation, do NOT quote SBI Re_k₂/Im_k₂/ζ/η. See amended `validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md`. Real slot replaced placeholder; headless load+sample verified |
| Titan MgSO₄ ocean | **awaiting artifact** | GUI placeholder; paths/centrals/gates unset |
| Titan NaCl ocean | **awaiting artifact** | GUI placeholder; paths/centrals/gates unset |
| Galileo–Europa v1.1 8D | DEPLOYED 2026-07-20 | gates 8/8; PPC 2026-08-04: 0/5 flagged, no tidal warning (tidal channels non-informative by construction) |
| Clipper–Europa v4 geodesy 11D | DEPLOYED 2026-07-21 | user-ratified; PPC 2026-08-04: 0/21 flagged, no tidal warning (Re_k2 informative extreme-tail datum, SBI tracks MCMC 1.43σ update to 0.04σ — decisive clean evidence) |
| Clipper–Europa v5 thick-ice ablation trio | **delivered, NOT ratified — GUI-gated** | baseline gates: sbc 0, limits 2, crosscheck 2; nok2 arm sbc 2; slots render hold warning, conditioning disabled (`cb597490`). PPC 2026-08-04 (baseline): 0/21 flagged, non-trivially clean (C20 22.6σ update tracked 0.12σ); flow-fidelity only, does NOT clear ratification blocks |
| Clipper–Europa v6 freegrav trio | **delivered, NOT ratified — GUI-gated** | baseline gates: sbc 0, limits 2, crosscheck 2; noinduction arm sbc 2; was silently the Clipper default — now gated, v4 default restored (`cb597490`). PPC 2026-08-04 (baseline): 0/20 flagged (gravity loosened by design 0.4–0.5σ, Re_k2 1.66σ + synodic 2.8–3.1σ tracked ≤0.07σ); flow-fidelity only |
| Clipper–Europa v7 open-\|Ae\| 11D | **delivered, NOT ratified** | gates: sbc 0, limits 2, crosscheck 2; never wired into GUI. PPC 2026-08-04: 0/21 flagged (Bind_synodic_x_real 52.9σ update tracked 0.18σ); gaps batch-largest 0.10–0.27σ carried to B5; PPC orthogonal to the B3 parameter-space reference-wander block |
| v1 / v2 | RETIRED | replaced by v1.1 / v4 |
| v3 | VETOED | wrong k2 bounds; 2D cache reused by v4+ |

## ADJUDICATED 2026-08-02 — v5/v6/v7 gates (record: plans/active/europa-v5v6v7-gate-adjudication.md)

Manager + independent Opus-5 adversarial review. Outcomes:

- **Limits FAILs**: superseded-gate artifact (pre-`3fa5a3fc` code demanded
  containment==1.0 pooled over off-manifold sweep points; unreachable for
  NSF flows, deploy path rejects out-of-box). Provisionally overturned —
  official regeneration at HEAD with the constructed sweep grid + anchor
  mode required (B1/B6) before any PASS is recorded.
- **v5: NOT ratifiable.** D_iceIh_km shape-EXCESS failure survives HEAD's
  own decomposition (1.12x tol, headline science param) + zeta_Ih mean
  1.074x + SBC underpowered (n=108 < 200 spec minimum).
- **v6: conditionally ratifiable.** Clean pass at HEAD (0 shape-excess,
  0 mean fails); blocked only on adequately powered SBC (n=102) + report
  regeneration (B1/B2/B6).
- **v7: crosscheck FAIL adjudicated a REFERENCE-side artifact.** True v7
  posterior provably identical to v5's at the fiducial (support strictly
  enlarged; newly admitted region carries 0.000000 posterior mass in both
  references); the 1.06 km reference disagreement has the wrong sign for a
  support effect and the rigid-translation shape of nested-sampling
  log-volume error. v7 flow passes 22/22 mean+sigma gates against the v5
  reference. Which reference wandered is OPEN (v5's log_Z_err 0.281 makes
  it the more suspect run). Blocked on B3/B4/B5.
- **Top open scientific risk:** every campaign with statistical power flags
  the same D–w degeneracy pair (v4 SBC FAIL on Tb+log10_w at n=1000; v5
  shape-excess on D_iceIh; v7 means on D_iceIh+log10_w). B2/B5 test it.
- **Governance:** deployed v4's ledger row corrected — its own reports
  contain undisclosed crosscheck mean breaches (zeta_Ih 1.23x, rho_core
  1.18x) and an SBC FAIL (Tb p=4e-4, log10_w p=0.034, n=1000). Whether v4
  needs re-adjudication is a USER decision.

GUI gating (`cb597490`) stays for all v5/v6 slots; v7 remains unwired.
Machine B follow-up queue B1–B7 in `plans/MACHINE-B-HANDOFF.md` §0.

## NH3 activity model corrected (2026-07-28, commit 6c5ee2af)

CoolProp mixture excess had the WRONG SIGN (under-depressed liquidus 9–37%).
Replaced by Melinder-anchored Redlich-Kister on the Choukroun & Grasset P,T
shape (`NH3_ACTIVITY_MODE='melinder-CG'`, default). "L-K 2002" polynomial
DELETED (unattributable; dilute slope 2x too shallow). Consequence: Titan NH3
campaign Phase 0 MUST re-run under the corrected model; provisional Phase-1
rectangle Tb [248,257] K x w [30,100] ppt. Spec:
`plans/active/titan-nh3-ocean-campaign-spec.md`.

## Titan NH3 JOINT no-ocean+ocean artifact — SPLIT ratification (2026-08-03)

Task #68 delivered: frozen (Tb,w) grid nodes build as REAL no-ocean interiors
and appear in the posterior (not rejected), per the user's 2026-07-30 decision.
Full campaign complete; reviewer status is SPLIT (gravity/structure verified,
tidal sector not — see the PPC subsection above):

- Cache `Test54_nh3_ocean/titan_nh3_joint_structure_grid_2d.pkl` (sha256
  `3d837cf8…`, `retry_frozen_as_no_ocean=True`); 1M gen 689,845 kept; reference
  MCMC 4461 samples; 1M nsf flow `titan_freegrav_nh3_posterior_1m.pt`.
- Gates: **SBC PASS** (353 pairs, min BH-adj KS p 0.772). limits FAIL is N/A for
  a joint mixture (HP ices carry tidal response in the no-ocean branch;
  documented, NOT tuned). crosscheck FAIL is **NOT benign-nuisance-only** — the
  PPC proved the per-param gate is blind to a joint tidal-sector under-update.
- **Split status:** gravity/structure (C20, C22, R_core, rho_core, salinity, Tb,
  dC20/dC22 null) VERIFIED; tidal (k2/zeta/eta) NOT verified — MCMC reference
  authoritative. Amended verdict + PPC finding in
  `validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md`.
- **MgSO4/NaCl inherit the split-status discipline + the new PPC/pushforward-gate
  requirements** (NOT the old benign-nuisance framing). Committed 2026-08-03.

GUI: NH3 joint slot wired into `_SBI_ARTIFACT_SLOTS`
(`PlanetProfileApp/pages/Inference.py`) 2026-08-02 as the newest Cassini–Titan
version; widget keys already namespaced by artifact filename (no shared-key
collision). Headless load+sample verified (artifact loads at deployed path,
2000 finite draws at the slot fiducial, salinity+Tb span the full joint plane).
In-browser render (`implemented, unverified`) still pending.

**BLOCKER found 2026-08-03 (Machine A, user-reported):** the structure cache
`Test54_nh3_ocean/titan_nh3_joint_structure_grid_2d.pkl` was NOT committed —
the slot errors ("Structure cache not found") on every machine but B, and as
the newest Titan version it had become the DEFAULT, breaking the Titan model
selector for users. Repairs on `255585c0` (verified, AppTest 8/8): version
default skips slots whose artifact/cache is missing locally (NH3 stays
selectable with awaiting-cache guidance, Generate disabled); NH3 tag added to
the dropdown-visible label tail; wedge now labels the ocean with the run's
asserted composition (NH3) instead of Titan's MgSO4 body default, threaded
via the slot -> result config (training JSON/config hash untouched).
**Machine B: commit + push the Test54 cache pkl** (v5 precedent, 21 MB) —
until then the NH3 slot is non-functional on Machine A and any HF deploy.
Manager countersign of the B-side ratification + in-browser verify follow
once the cache lands.

Next in Phase B: MgSO4 (seeds 73/7373/73, salinity log10[0,2.288]), then NaCl
(74/7474/74, log10[0,2.477]) — each own config + 2D cache + artifact.

**Posterior-predictive + interior diagnostic (2026-08-03, in gate runner).**
`plans/scripts/titanG_ppc_interior_check.py` (wired into `titanG_nh3_run_gates.py`
→ inherited by MgSO4/NaCl gate scripts) pushes SBI posterior draws back through
the identical theta→x forward loop (new `theta_override` kwarg on
`generate_sbi_dataset`; forward path validated exact, max|fwd−stored|=0) and
reports posterior-predictive coverage vs x_obs + the prior-predictive percentile
of x_obs. NH3 result: x_obs interior on all 4 channels (C20 42/C22 60/Re_k2 80/
Im_k2 86 pctile) — the flow interpolates, confirming the prior verdict. BUT a
NEW finding surfaced: SBI posterior-predictive |Im_k2| median 0.042 vs the MCMC
reference 0.093 (obs 0.135) — the SBI flow under-updates the tidal sector,
traced to log10_zeta/eta_Ih/eta_III shifts that pass the per-param crosscheck by
the thinnest margins (90–92% of tol) and compound in the nonlinear k2
pushforward. The per-param crosscheck gate provably cannot see this (it never
pushes forward to k2).

**Reviewer verdict (2026-08-03): ratification SPLIT.** Gravity/structure sector
VERIFIED; tidal (k₂/dissipation) sector NOT verified — MCMC reference
authoritative, do not quote SBI Re_k2/Im_k2/zeta/eta at the Titan datum. SBC
(prior-averaged) PASS coexists with local miscalibration at the informative
datum (no paradox). RATIFICATION.md, gate_summary.json, GUI scope_note, and the
`project_joint_mixture_gate_scope` memory all amended to the split status (the
old "benign eta nuisance" framing is superseded). **MgSO4/NaCl requirements
before ratification:** run the PPC + prior-pred interior recheck; promote a
pushforward-observable crosscheck to a FIRST-CLASS gate (must flag the known NH3
k2 miss); diagnose the flow under-update before more compute (don't assume "more
epochs" fixes it). **Do not start MgSO4 compute until the flow-diagnosis
approach is chosen.** Report: `validation_reports/titan_freegrav_nh3_1m/ppc/`.
Deployed sampling path confirmed `reject_outside_prior=True` (reviewer MINOR
closed).

## Landed on genai since the 8e91d440 HF snapshot (selection)

- NH3: CoolProp NH3-H2O ocean EOS; Melinder-anchored liquidus correction +
  rewritten `tests/coolprop_nh3_test.py` (15 green).
- GUI: heating-tab radiogenic inventory selector (BSE/CI + age slider);
  Mineralogy tab (post-hoc Perple_X grain-density consistency, tolerance band
  + porosity headroom + cold-edge flags) — both AppTest-verified. The
  Perple_X native-domain density heatmap + selected-draw geotherm overlay
  (Codex C3) is `verified` 2026-08-02: AppTest + manager visual inspection
  of rendered CV3/CI figures. Titan amortized Phase-A/NH₃/MgSO₄/NaCl slot
  scaffolding (Codex C4) is `verified` 2026-08-02 (manager AppTest: Phase-A
  default, placeholder states, disjoint per-slot widget keys). Follow-on
  manager fix `cb597490` gates the unratified v5/v6 slots (adjudication
  hold enforcement; suite 7/7). Perple_X cold-edge band/legend polish (C6)
  is `implemented, unverified`; CV3/CI PNGs and AppTest evidence are recorded
  in `plans/CODEX-QUEUE.md`.
- Fixes: Titan "C/MR^2 = None" message (SetCMR2strings guard + SetupInit
  call); PREM porosity-table path; NaCl `_warn_once` broadcast crash;
  `bulk_overrides` threading in `build_tbw_grid_cache`.
- Repro: frozen PPcl environment.yml; scipy 1.17.1 / numpy 2.4.6 bumps;
  v5/v6/v7 artifacts, caches, reference results, gate reports.

## Public deployment boundary

Hugging Face Space (vsteven-planetprofile.hf.space) is CURRENT at genai
`d53385f1`, deployed 2026-08-04 (user ran the upload one-liner; USER
confirmed the Space running). The batch replaces the 2026-07-21 `8e91d440`
snapshot and adds: Titan NH3 joint slot (+ its cache, sector warning, wedge
NH3 composition), v5/v6 not-ratified GUI gating with v4 as Clipper default,
mineralogy tab + Perple_X geotherm overlay, heating radiogenic-inventory
selector, ready-slot default fallback, NaCl `_warn_once` fix, and the NH3
activity-model correction underneath. Deploy-script cache list now includes
Test54 (`d53385f1`). C5 subsequently made future deploy snapshots derive and
exact-diff the full cache set from `_SBI_ARTIFACT_SLOTS`; build-only evidence
is `verified`, but no new deploy was pushed.

## Awaiting

- Machine A (Claude): v5/v6/v7 gate adjudication (scientific-reviewer pass);
  then wire ratified slots; maintain queues.
- Machine A (Codex 5.6): queue empty — C8 first doc-doctor pass is `verified`
  at `validation_reports/doc_doctor/2026-08-06_first_pass.md` (12/12 items,
  2 PASS / 10 FINDING; no fixes). Findings await manager triage into scoped
  follow-up tasks.
- **Machine A (Claude) — DECISION REQUESTED: NH3 flow under-update diagnosis
  is complete (`711c2fd2`); MgSO4/NaCl compute is HELD pending your call.**
  Diagnosis (reviewer PASS WITH CONCERNS, corrections folded): #3 x-norm scale
  RULED OUT; obs-vector width RULED OUT as causal (no-ocean control assimilates
  informative Re_k2 3.92σ→0.05σ, 0/4 flagged); mechanism localized to the
  ocean-admitting apparatus (joint mixture AND the co-varying salinity axis —
  NOT the mixture alone, since the control also adds a 13th param + removes the
  phase guard); bimodal mode-assignment is a CANDIDATE, not established; #1
  noise swamping cleared for Re_k2 but UNTESTED for the dominant Im_k2 miss.
  Three paths for you/the user to pick, cheapest-first: (1) run two cheap
  SEPARATORS before any build — an ocean-only-with-salinity PPC (splits salinity
  axis from bimodality) + the #1 reduced-noise Im_k2 pilot; (2) proceed with
  MgSO4/NaCl under split-status as-is (tidal sector MCMC-authoritative,
  pushforward gate fires by construction); (3) remediate the conditional first
  (mixture-aware embedding / has_ocean-labelled conditional / sequential round).
  Full report `validation_reports/FLOW_UNDERUPDATE_DIAGNOSIS.md`; rationale in
  `MACHINE-B-HANDOFF.md` §0.6 RESOLVED block. Machine B awaits the decision;
  will not start #1 pilot or any MgSO4/NaCl work unilaterally.
- Machine B: idle after the diagnosis push. HOLD on any v5/v6/v7 retraining
  until adjudication; MgSO4/NaCl HELD pending the decision above.
- USER: live-app click-throughs DONE 2026-08-05 — one change requested and
  shipped (`01649f3f`, verified): Europa induction conditioning defaults
  Re Ae = 0.9 for the synodic band + second harmonic (imaginary parts and
  the orbital band keep config-derived values; inputs stay editable). The
  Space runs the d53385f1 snapshot; the Ae-default change reaches it at the
  next redeploy. v5/v6 ratification (user priority) is executing via
  MACHINE-B-HANDOFF §0.7 — B3 done, B1/B2 against the fresh references
  next.

## Reference docs

- Machine B queue: `plans/MACHINE-B-HANDOFF.md`
- Codex lane: `plans/CODEX-QUEUE.md`, `AGENTS.md`
- Active/archive routing: `plans/active/README.md`, `plans/archive/README.md`
- Constraints per mission-body: `plans/mission-body-constraints.md`
- Gate reports: `validation_reports/` (+ `v5/v6/v7_gate_summary.json`)
- Artifact ledger: `PlanetProfile/Inference/sbi_artifacts/INDEX.md`
- NH3 defect trail: `plans/HANDOFF-2026-07-26-nh3-liquidus-defect.md`
