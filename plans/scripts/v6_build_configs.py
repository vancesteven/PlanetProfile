"""Build the Europa Clipper v6 'free-gravity' config trio from the v5 geodesy config.

v6 vs v5 geodesy channel (ONLY the gravity channel changes; interior priors,
cache, rheology, induction bounds, k2/h2 all identical to v5):
  1. DROP CMR2 from observables entirely (it was the double-count: 0.3547 was
     itself derived from C22 via the hydrostatic Radau relation).
  2. Set C20/C22 to GC21 Table 2 SOL-A (UNCONSTRAINED), unnormalized at 1565 km,
     NO sqrt(5)/sqrt(10/24) conversion (reviewer block-level correction):
        C20 = -J2 = -4.3759e-4 +/- 7.747e-5
        C22 =        1.3862e-4 +/- 2.44e-6
  3. Widen BOTH non-hydrostatic offset boxes so gravity carries ~ZERO interior
     C/MR2 information (USER decision 2026-07-22 "truly agnostic", resolving the
     scientific-reviewer MAJOR finding):
        dC20_nh: [-2e-5, 2e-5] -> [-3.9e-4, 3.9e-4]  (~5 sigma(C20)=7.747e-5)
        dC22_nh: [-2e-5, 2e-5] -> [-5e-5, 5e-5]
     Rationale: with dC22_nh at only +/-2e-5 and the TIGHT C22 likelihood
     (sigma=2.44e-6), C22_h+dC22_nh~=C22_obs pins the interior C22_h -> k_f ->
     C/MR2 into a box ~[0.335,0.373] centered on the hydrostatic reference
     (reviewer: interiors at C/MR2=0.30 need dC22_nh=+4.8e-5, C/MR2=0.40 need
     -4.7e-5, both rejected by +/-2e-5). That silently re-imposed a
     hydrostaticity-centered MoI constraint -- the very coupling the redesign
     removes. +/-5e-5 lets C/MR2 span ~[0.30,0.40] freely so gravity is genuinely
     agnostic and the interior comes from k2/h2/induction. The old justification
     conflated C22 MEASUREMENT PRECISION (2.44e-6) with the PRIOR WIDTH on the
     physical non-hydrostatic C22 anomaly (unrelated; can be several e-5).
     CONSEQUENCE: the published hydrostatic C/MR2 (GC21 0.3547) now returns ONLY
     as the optional inference-time reweight (plan task D) -- which is therefore
     REQUIRED before any science claim on the interior C/MR2 (D is inference-time,
     needs no retraining).

Ablation siblings mirror the v5 trio -- identical everywhere except the
observable subset:
  - geodesy_11D  : C20,C22 + k2/h2 + 14 Bind_  (20 obs)
  - noinduction  : C20,C22 + k2/h2             ( 6 obs)
  - nok2         : C20,C22 + 14 Bind_          (16 obs)

Writes directly to PlanetProfile/Inference/configs/ (NEW files, additive; v5
configs are NOT mutated). No compute here -- reviewable JSON only.
"""
import json, copy, os

CFGDIR = "PlanetProfile/Inference/configs"
V5 = os.path.join(CFGDIR, "europa_clipper_v5_geodesy_11D.json")

# --- GC21 Table 2 SOL-A (unconstrained), unnormalized at 1565 km ---
C20_SOLA = [-4.3759e-4, 7.747e-5]
C22_SOLA = [1.3862e-4, 2.44e-6]
# widened non-hydrostatic offset boxes (see module docstring; USER "truly agnostic")
DC20_BOX = [-3.9e-4, 3.9e-4]   # ~5 sigma(C20)
DC22_BOX = [-5.0e-5, 5.0e-5]   # lets interior C/MR2 span ~[0.30,0.40]; gravity agnostic

# metadata fields inherited from v5 that describe the OLD (CMR2-retained,
# narrow-box, synthetic-null) design and are FALSE for v6. We do not delete the
# history but prefix each with an explicit SUPERSEDED marker so a provenance
# reader / automated check is not misled (reviewer MODERATE finding).
_STALE_V5_FIELDS = [
    "description",
    "gravity_forward_model_2026_07_19",
    "gravity_fiducial_2026_07_19",
    "primary_constraint",
    "non_hydrostaticity_deliverable_2026_07_19",
    "fiducial_recompute_2026_07_20",
]
_SUPERSEDE = ("[SUPERSEDED FOR v6 2026-07-22 -- describes the v5 design; see "
              "v6_freegrav_redesign_2026_07_22. In v6: CMR2 is NOT an observable; "
              "gravity is agnostic (widened dC20_nh/dC22_nh boxes -> ~zero interior "
              "C/MR2 info); 11 sampled parameters; gravity is conditioned on REAL "
              "GC21 SOL-A data, not a synthetic hydrostatic null.] ")

with open(V5) as f:
    v5 = json.load(f)

v6 = copy.deepcopy(v5)

# 1. remove CMR2 observable; set C20/C22 to SOL-A. Preserve remaining order.
obs = v6["observables"]
obs.pop("CMR2", None)
obs["C20"] = C20_SOLA
obs["C22"] = C22_SOLA

# 3. widen BOTH offset boxes
v6["param_space"]["dC20_nh"]["bounds"] = DC20_BOX
v6["param_space"]["dC22_nh"]["bounds"] = DC22_BOX

md = v6["metadata"]

# mark stale v5-design fields as superseded (reviewer MODERATE)
for k in _STALE_V5_FIELDS:
    if k in md and isinstance(md[k], str):
        md[k] = _SUPERSEDE + md[k]

md["schema_version"] = ("v6 (v5 interior/induction/k2 + FREE-GRAVITY geodesy: "
                        "CMR2 dropped, C20/C22 = GC21 SOL-A unconstrained, "
                        "agnostic gravity) — clipper_v6_freegrav_11D")
md["v6_freegrav_redesign_2026_07_22"] = (
    "USER-DIRECTED (2026-07-22): eliminate the geodesy double-count. v5 imposed "
    "CMR2=[0.3547,0.0024] as a TIGHT independent observable AT THE SAME TIME as "
    "tight C20/C22 -- but 0.3547 was itself derived FROM C22 via the hydrostatic "
    "Radau relation (GC21 Table 3), so the same physical relation was used twice, "
    "biasing the posterior toward hydrostaticity (the very quantity the campaign "
    "measures). v6 fix: (1) CMR2 REMOVED from observables (no longer a likelihood "
    "term); (2) C20/C22 set to GC21 Table 2 SOL-A UNCONSTRAINED (J2=437.59+/-77.47, "
    "C22=138.62+/-2.44, x1e-6, UNNORMALIZED at R_ref=1565 km -- NO sqrt(5)/sqrt(10/24) "
    "conversion per reviewer block-level correction, since Table 2 is already "
    "unnormalized). SOL-A -C20/C22 = 3.157 (NOT the hydrostatic Tricarico 3.324): "
    "the unconstrained solution does not obey the hydrostatic ratio -- that ~0.3-sigma "
    "gap IS the (statistically insignificant) non-hydrostatic signal, consistent with "
    "GC21's own J2/C22=3.16+/-0.57. RD hydrostatic-reference C/MR2 from SOL-A C22 = "
    "0.3559 (consistent with the v5 fiducial 0.3550). (3) BOTH non-hydrostatic offset "
    "boxes widened so gravity carries ~ZERO interior C/MR2 info (USER 'truly agnostic', "
    "resolving reviewer MAJOR finding): dC20_nh [-2e-5,2e-5]->[-3.9e-4,3.9e-4] "
    "(~5*sigma(C20)); dC22_nh [-2e-5,2e-5]->[-5e-5,5e-5]. The narrow v5 dC22_nh + tight "
    "C22 silently boxed interior C/MR2 to ~[0.335,0.373] centered on the hydrostatic "
    "reference (a residual hydrostaticity-centered MoI constraint); +/-5e-5 lets C/MR2 "
    "span ~[0.30,0.40] freely. The old 'sigma(C22) tiny' rationale conflated C22 "
    "measurement precision with the prior width on the physical non-hydrostatic C22 "
    "anomaly. CONSEQUENCE: the published HYDROSTATIC C/MR2 (GC21 0.3547+/-0.0024) is NO "
    "LONGER an observable; it returns ONLY as the OPTIONAL inference-time reweight on "
    "the derived hydrostatic-reference C/MR2 (plan task D), which is therefore REQUIRED "
    "before any science claim on the interior C/MR2 (D is inference-time, needs no "
    "retraining). The PRIMARY non-hydrostaticity deliverable is UNAFFECTED by the box "
    "widening: the identifiable combination dC20_nh + 3.324*dC22_nh is interior-"
    "independent (C22_h cancels), central estimate unbiased. Dual C/MR2 deliverable: "
    "hydrostatic-reference (RD from C22) vs actual (structure integral); gap = "
    "non-hydrostaticity. Plan: fluffy-snacking-fountain. Interior priors, 2D (Tb x w) "
    "cache, rheology, induction bounds, and k2/h2 conditioning are BYTE-IDENTICAL to v5; "
    "only the gravity channel + the two offset-box widths changed.")
md["v6_hybrid_conditioning_2026_07_22"] = (
    "DESIGN NOTE (reviewer MODERATE): v6 conditions the GRAVITY channel (C20/C22) on "
    "REAL Galileo-era data (GC21 SOL-A, ratio 3.157), while k2/h2 and the 14 Bind_ "
    "induction channels remain the v5 SYNTHETIC fiducial forecast (Re_k2=0.23 Mazarico "
    "projection; Bind from the Tb=264.5 K / w=35.165 ppt / R_core=534.67 km fiducial "
    "interior). v6 is therefore a real-gravity + forecast-tidal/induction HYBRID, NOT "
    "the v5 self-consistent hydrostatic NULL. Verified mutually consistent: the fiducial "
    "interior's C22_h=1.3775e-4 sits 0.36 sigma from real SOL-A C22=1.3862e-4, and "
    "fiducial C/MR2~=0.355 ~= SOL-A hydrostatic-reference 0.356, so the hybrid does not "
    "manufacture tension. The 'hydrostatic NULL' language in superseded v5 fields does "
    "NOT apply to the v6 gravity channel.")
md["v6_gravity_provenance_2026_07_22"] = {
    "source": "Gomez Casajus et al. 2021 (GC21) Table 2, SOL-A (unconstrained)",
    "convention": "unnormalized Stokes coefficients at R_ref = 1565 km; "
                  "NO sqrt(5)/sqrt(10/24) fully-normalized conversion applied "
                  "(Table 2 is already unnormalized -- reviewer block-level fix, "
                  "verified against papers/gomescasajus2021updated.pdf Table 2)",
    "C20": {"value": C20_SOLA[0], "sigma": C20_SOLA[1], "note": "= -J2, J2=437.59+/-77.47 x1e-6"},
    "C22": {"value": C22_SOLA[0], "sigma": C22_SOLA[1], "note": "138.62+/-2.44 x1e-6"},
    "ratio_minus_C20_over_C22": 3.1568,
    "hydrostatic_ref_cmr2_from_C22_rd": 0.35586,
    "cmr2_removed_reason": "double-count: 0.3547 derived from C22 via hydrostatic "
                           "Radau; imposing both CMR2 and C22 uses the relation twice",
    "n_sampled_params": 11,
    "offset_boxes_widened": {
        "dC20_nh": {"old": [-2e-5, 2e-5], "new": DC20_BOX,
                    "reason": "~5*sigma(C20)=7.747e-5; carries honest C20 uncertainty "
                              "and the non-hydrostatic upper-limit tail without a ~3-sigma wall"},
        "dC22_nh": {"old": [-2e-5, 2e-5], "new": DC22_BOX,
                    "reason": "removes the hidden C/MR2 box [0.335,0.373]; lets interior "
                              "C/MR2 span ~[0.30,0.40] so gravity carries ~zero interior info"},
    },
    "task_D_required": ("With gravity agnostic, the published hydrostatic C/MR2 "
                        "(GC21 0.3547+/-0.0024; or Anderson 0.3475; Jacobson 0.3405) is "
                        "applied ONLY as the optional inference-time reweight (plan task "
                        "D). Task D MUST land before any interior-C/MR2 science claim."),
    "S22_note": "S22=-6.21+/-2.90, C21/S21 available but degree-2 tesseral, out of "
                "current C20/C22-only scope",
    "correlation_note": "SOL-A rho(C20,C22)=-0.17 (safe, <1); not yet imposed -- "
                        "diagonal noise retained (effect on the non-hydrostatic-"
                        "combination sigma ~1.8%, negligible). If adopted, cap rho<1 "
                        "(mcmc_runner bivariate branch divides by det; SOL-B rho=1.0 "
                        "would be singular).",
}
# supersede v5-specific 'moi_prior' framing note (CMR2 no longer an observable)
if "moi_prior" in md:
    md["moi_prior"]["v6_status_2026_07_22"] = (
        "SUPERSEDED in v6: CMR2 is NOT an observable in v6, and gravity is agnostic "
        "(widened offset boxes). This Galileo MoI returns ONLY as the OPTIONAL "
        "inference-time reweight (plan task D, required before any interior-C/MR2 "
        "science claim), never a likelihood term.")

# ---- write the three arms (differ ONLY in observables subset) ----
def write_arm(fname, keep):
    arm = copy.deepcopy(v6)
    arm["observables"] = {k: v6["observables"][k] for k in keep}
    arm["metadata"] = copy.deepcopy(md)
    arm["metadata"]["arm_2026_07_22"] = {
        "arm": fname, "n_obs": len(keep), "obs_names": keep,
        "note": "v6 ablation sibling; identical to v6 geodesy except observable subset. "
                "Dataset generated once from the full-obs baseline then NAME-sliced "
                "(mirrors plans/scripts/v6_gen_dataset.py). NOTE: the noinduction arm "
                "still carries induction_bounds (|Ae|>0.7 synodic support cut) applied "
                "at dataset generation, so it is 'no induction LIKELIHOOD', not "
                "induction-free -- the ablation contrast is partially conditioned on the "
                "residual Galileo induction support prior (mirrors v5)."}
    out = os.path.join(CFGDIR, fname)
    with open(out, "w") as f:
        json.dump(arm, f, indent=2)
    print(f"[v6] wrote {out}  ({len(keep)} obs): {keep}")
    return keep

all_names = list(v6["observables"].keys())        # 20: C20,C22,k2/h2,14 Bind
k2h2 = ["C20", "C22", "Re_k2", "Im_k2", "Re_h2", "Im_h2"]
bind = ["C20", "C22"] + [n for n in all_names if n.startswith("Bind_")]

write_arm("europa_clipper_v6_freegrav_11D.json", all_names)
write_arm("europa_clipper_v6_freegrav_noinduction_6obs.json", k2h2)
write_arm("europa_clipper_v6_freegrav_nok2_16obs.json", bind)
print("[v6] done. v5 configs untouched.")
