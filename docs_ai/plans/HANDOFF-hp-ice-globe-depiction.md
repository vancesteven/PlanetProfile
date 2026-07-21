# Handoff: depicting high-pressure ices (III / V / VI) in the interactive globe

**Author:** genai session, 2026-07-20
**For:** machine A (implementer)
**Status:** `implemented, verified` (Machine A, 2026-07-20) — HP-ice shell now
emitted to the globe. Verified on the Titan Test50 slot (median HP thickness
319 km; light-blue shell between the Ih surface and the silicate interior at
r = 2417 km > silicate 2081 km; Europa correctly shows no HP shell).



**UPDATE 2026-07-20 (per-phase):** superseded the single aggregate hp_ice shell — each ice phase (Ih / III / V / VI), the ocean, silicate mantle, and core now render as its OWN labeled shell with finite proportional thickness and a color legend (Ih=textured surface skin, III=cyan, V=teal, VI=violet, ocean=royal blue). Carried via a new InferenceResult.D_icePhase_results dict {III,V,VI}. Verified: Ganymede (all 5 phases nested + legend), Titan (Ih/III/V/VI, no ocean), Europa (ocean+silicate, no HP).

(Original spec below.)

**Status of the depiction:** `not implemented` — the HP-ice *style* exists but no
HP-ice layer is ever emitted to the globe. This document is the findings + a
concrete, grounded implementation spec. It changes NO fundamental scientific
assumption; it only surfaces cache quantities that are already computed.

---

## TL;DR

The globe renderer already knows how to *color* a high-pressure-ice shell
(`LAYER_COLORS['hp_ice']`), but nothing in the data path ever hands it an
`hp_ice` layer. The per-phase HP-ice thicknesses **are computed and stored in the
structure cache** (`D_iceIII_km`, `D_iceV_km`, `D_iceVI_km`) but are **not
propagated to the `InferenceResult`**, so the globe cannot draw them. The fix is
two small, additive steps: (1) package an aggregate HP-ice thickness per posterior
sample onto the result, and (2) insert an `hp_ice` shell between the ocean and the
silicate floor in the two layer-builders.

**Important scope note (Europa):** on the current Europa GSW-seawater cache, HP
ices are **structurally impossible** (basal-ocean P ≈ 156 MPa < the 200 MPa
`PminHPices_MPa` cap — see `cache_builder.py` docstring ~line 812). So for Europa
v3/v4/v5 the HP-ice shell will correctly render as **zero thickness / absent**.
This feature matters for **Ganymede, Callisto, Titan** (thick HP-ice mantles) and
is worth building now so those bodies render correctly. Verify on a Ganymede or
Titan result, not Europa.

---

## Where the data already exists

`PlanetProfile/Inference/cache_builder.py` computes and stores, per structure node:

- line ~515-517: `D_iceIII_km`, `D_iceV_km`, `D_iceVI_km`
- line ~513: `D_hsphere_km` (full hydrosphere: Ih + HP ices + ocean)
- `D_iceIh_km`, `D_ocean_km`

The hydrosphere decomposition is exactly:
`D_hsphere_km = D_iceIh_km + D_iceIII_km + D_iceV_km + D_iceVI_km + D_ocean_km`
(see `cache_builder.py:525`).

These are retrievable at inference time via
`MCMCRunner._get_cache_scalar(theta_dict, key)` (line 1153) — the SAME accessor
already used for `D_ocean_km`, `D_iceIh_km`, `D_hsphere_km`.

## Where the data is dropped (the gap)

Both runners currently package only `D_ocean_results`, `D_iceIh_results`,
`D_hsphere_results` onto the `InferenceResult` — never the HP-ice thicknesses:

- `mcmc_runner.py:1840-1842` (per-sample loop) and `:1890-1892` (result assembly)
- `sbi_runner.py:1211-1213` and `:1258-1260`
- `inference_core.py:264-281` — `InferenceResult` dataclass has
  `D_ocean_results`, `D_iceIh_results`, `D_hsphere_results` but **no**
  `D_iceHP_results` field.

## Where it would be drawn (also missing)

`PlanetProfileApp/Utilities/globe_view.py`:
- `LAYER_COLORS['hp_ice'] = 'rgb(170, 200, 230)'` (line 65) — style ready.
- `posterior_median_layers()` (lines 297-337) builds core → silicate → ocean-top.
  It **never appends an `hp_ice` layer**, and it does not read any HP-ice field.

`PlanetProfileApp/pages/Inference.py`:
- `_sample_layers(i)` (lines 2432-2449) mirrors `posterior_median_layers` for the
  click-to-select-a-sample path — also missing the HP-ice shell.

So there are **two** layer-builders to update in lockstep (median path + sample path).

---

## Radial geometry (how to place the HP-ice shell)

Layers are concentric spheres given by `r_km` (radius from center), drawn
inner→outer. Working from the surface inward, radii are:

```
R_body                              (surface)
R_body - D_iceIh                    = base of Ih  = ocean top      (kind 'ocean')
R_body - D_iceIh - D_ocean          = base of ocean = HP-ice top   (kind 'hp_ice')  <-- NEW
R_body - D_hsphere                  = base of hydrosphere = silicate floor (kind 'silicate')
R_core                              (kind 'core')
```

Because `D_hsphere = D_iceIh + D_ocean + D_iceHP` (HP ices sit between the ocean
and the silicate floor), the HP-ice shell's **outer radius** is
`R_body - D_iceIh - D_ocean` and its **inner radius** is `R_body - D_hsphere`
(= the existing silicate-floor radius). In the concentric-sphere model each layer
is a single `r_km` (its outer radius); the renderer fills inward to the next
smaller shell. So the new layer's `r_km = R_body - D_iceIh - D_ocean`, inserted
**after** the ocean-top layer and **before** the silicate layer in the
inner→outer sorted list (the renderer sorts by `-r_km`, so ordering in the append
list does not matter — only the radius value does).

Aggregate the three HP phases into one shell for display:
`D_iceHP = D_iceIII + D_iceV + D_iceVI`. (A future refinement could draw III/V/VI
as separate tinted shells, but a single aggregate `hp_ice` shell is the right
first step and matches the single `LAYER_COLORS['hp_ice']` entry.)

---

## Implementation steps (all additive, no behavior change for existing bodies)

### Step 1 — carry HP-ice thickness onto the result

**`inference_core.py`** (~line 275, next to `D_hsphere_results`): add field
```python
D_iceHP_results: Optional[np.ndarray] = None  # aggregate III+V+VI thickness (km), per sample
```
Document it in the class docstring alongside `D_iceIh_results`.

**`mcmc_runner.py`**:
- In the per-sample loop (~1820 init, ~1840 append), add a `D_iceHP_results` list
  and append `sum of _get_cache_scalar(theta_dict, k) for k in
  ('D_iceIII_km','D_iceV_km','D_iceVI_km')` (treat missing/NaN as 0.0 so bodies
  without those keys are unaffected).
- Convert to array (~1855) and pass `D_iceHP_results=...` into the
  `InferenceResult(...)` construction (~1892).

**`sbi_runner.py`**: mirror exactly — init `np.full(n_samples, np.nan)` (~1184),
fill in the eval loop (~1212), pass into the result (~1259).

Guard for older caches: if any of the three keys is absent, `_get_cache_scalar`
should yield NaN → treat as 0.0 in the sum (so the aggregate is 0.0, not NaN,
when a body genuinely has no HP ices — e.g. Europa). Confirm `_get_cache_scalar`'s
missing-key behavior and coalesce accordingly.

### Step 2 — emit the hp_ice shell in BOTH layer-builders

**`globe_view.py::posterior_median_layers`** (after the ocean-top append,
lines 333-336): compute `d_hp = _med('D_iceHP_results')`; if `d_hp` is finite and
`> 0.5` km, append
```python
{'name': 'high-pressure ice (top)',
 'r_km': R_km - d_ih - d_oc,
 'kind': 'hp_ice'}
```
(Only when `d_ih` and `d_oc` are both available, since the radius needs them.)

**`Inference.py::_sample_layers`** (after the ocean append, lines 2444-2448):
mirror with the per-sample arrays — pull a `d_hp` array the same way `d_ih`/`d_oc`
are pulled (add `d_hp = getattr(result, 'D_iceHP_results', None)` near where
`d_ih`, `d_oc`, `d_hs` are unpacked, ~line 2353), and append the `hp_ice` dict
guarded by `d_hp.size > i and np.isfinite(d_hp[i]) and d_hp[i] > 0.5`.

### Step 3 — (optional, nice) legend/annotation

The globe has no legend (`showlegend=False`); layers self-describe via
hovertemplate. The HP-ice shell will hover as `high-pressure ice (top): r = NNN km`
automatically. The matplotlib layer-thickness panel in `Inference.py`
(~lines 1889-1952) already notes "hydrosphere includes any HP-ice layers"; consider
adding an explicit "HP ice (III+V+VI) thickness" histogram panel there using the
new `D_iceHP_results` (parallel to the Ice Ih panel at line 1900) so the number is
visible outside the globe too.

---

## Verification (per project UI-verification discipline)

This is a Streamlit/Plotly change → `implemented, unverified` until rendered
in-browser. Required checks:

1. **Do NOT verify on Europa** — HP ices are structurally absent there, so the
   shell correctly won't appear (a valid but uninformative result). Verify on a
   **Ganymede, Callisto, or Titan** inference result that has nonzero
   `D_iceIII/V/VI` in its cache.
2. Launch the app, open that result's "🌐 Interactive Globe" expander, and confirm
   a light-blue `hp_ice` shell renders **between** the (darker-blue) ocean and the
   (brown) silicate floor, at radius `R - D_iceIh - D_ocean`.
3. Rotate into the cutaway wedge; confirm the shell's hover reads
   `high-pressure ice (top): r = … km` and the radius equals
   `R_body − D_iceIh − D_ocean` for that sample.
4. Confirm Europa results are **unchanged** (no spurious HP shell; regression check).
5. Numerical cross-check: for the selected sample,
   `D_iceIh + D_ocean + D_iceHP + (silicate thickness) ≈ D_hsphere` and the shell
   radii are strictly nested (no shell pokes through another) at exaggeration = 1.

## Files to touch (summary)

| File | Change |
|------|--------|
| `PlanetProfile/Inference/inference_core.py` | add `D_iceHP_results` field + docstring |
| `PlanetProfile/Inference/mcmc_runner.py` | populate + pass `D_iceHP_results` (2 spots) |
| `PlanetProfile/Inference/sbi_runner.py` | populate + pass `D_iceHP_results` (2 spots) |
| `PlanetProfileApp/Utilities/globe_view.py` | emit `hp_ice` layer in `posterior_median_layers` |
| `PlanetProfileApp/pages/Inference.py` | emit `hp_ice` layer in `_sample_layers`; unpack `d_hp`; optional HP-ice histogram panel |

## Non-goals / cautions

- Do **not** split III/V/VI into separate shells in this pass — one aggregate
  `hp_ice` shell first (matches the single style entry). Per-phase tinting is a
  clean follow-on if desired (the per-phase cache fields are already there).
- Do **not** change `cache_builder.py` — the data is already computed and stored.
- Do **not** retrofit HP ices onto the Europa seawater path — they are physically
  absent there by construction; the guard (`> 0.5 km`) handles this.
- Keep the change additive: older result pickles lacking `D_iceHP_results` must
  still render (the `getattr(result, 'D_iceHP_results', None)` + finite/`>0.5`
  guards ensure the shell is simply omitted, not a crash).
