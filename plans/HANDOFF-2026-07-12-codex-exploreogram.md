# HANDOFF → codex: Exploreogram fixes (state + remaining work)

Branch `genai`, repo /Users/svance/ppgenai. Env: conda `PP` on this machine
(`source /opt/miniconda3/etc/profile.d/conda.sh && conda run -n PP ...`). App:
`cd PlanetProfileApp && conda run -n PP streamlit run PlanetProfileApp.py`.
House rules: repo CLAUDE.md — status vocabulary (`verified` / `implemented,
unverified` / `not implemented`); matplotlib output must be visually inspected
before claiming a plot fix; no Test/ edits without permission; no push without
approval. Wider campaign context: plans/HANDOFF-2026-07-09-test50-sbi-validation.md.

## User-reported Exploreogram issues (2026-07-12)

1. Run fails: `Previous run failed: No module named 'configPP'`.
2. No planet selected → page should default to Europa.
3. Default Europa preset should be: x = mean ocean conductivity (salinity driver),
   y = ocean thickness (ice-shell-thickness driver), z = induced magnetic field with
   synodic + synodic 2nd (harmonic) + orbital excitations, log axes, Bx component.
4. Matplotlib exploreogram plots show a jagged/oscillatory appearance along the
   y-axis at high salinity.

## DONE this session (commit accompanying this file) — status: verified via
scripted AppTest render (streamlit.testing.v1); user click-through pending

- **(1) configPP import**: bare `from configPP import configAssign` only resolved
  when the CWD held `UserConfigs/` (app CWD is PlanetProfileApp/ → ModuleNotFound).
  Changed to the qualified `from UserConfigs.configPP import configAssign` (repo
  root is on sys.path; UserConfigs is a namespace package) at FOUR sites:
  Exploreogram.py (x2), RunPlanetProfile.py, FigureSettings.py. Same disease as the
  earlier PlanetProfile/__init__ + SPICE + SpacecraftMAGdata CWD fixes.
- **(2) Europa default body**: the hard `st.stop()` when no Planet is selected now
  auto-loads Europa via `Utilities.PlanetLoader.load_planet_module(parent_directory,
  'Europa')` with an info banner (fallback to the old warning if loading fails).
- **(3) Europa default preset**: session-state seeding (only when keys absent, so
  user choices persist): explore_xName='sigmaMean_Sm' (driver default wOcean_ppt is
  automatic), explore_yName='D_km' + y_driver_select='zb_approximate_km',
  explore_zName='Amp_nT', exploreogram_x_log/y_log=True,
  exploreogram_induct_component='Bx', selected_excitations=['synodic',
  'synodic 2nd', 'orbital'], x_min=1.0 (positive salinity floor so the log axis is
  honored). Also fixed a pre-existing bug: the z-variable selectbox hardcoded
  index=0, stomping any seeded/persisted z selection on first render.
- AppTest verification: page renders with exactly the requested defaults and no
  exceptions (xName sigmaMean_Sm / yName D_km / zName Amp_nT / y driver
  zb_approximate_km / logs True/True / component Bx / excitations synodic +
  synodic 2nd + orbital).

## Item (4) jagged y-axis plots — status: implemented, unverified (Machine B verifies)

Machine A implemented (2026-07-13, this commit), per the recon below:

- **Fold-over guard** `_sort_grid_axes` in PlanetProfile/Plotting/ExplorationPlots.py,
  applied in `PlotExploreOgram` right after the VALID/NaN masking: any grid row/column
  whose (derived) axis coordinates are non-monotone is argsorted along that direction
  with z + contour carried along. Rectilinear grids pass through untouched (CLI
  behavior identical for normal axes). Mechanism VERIFIED visually on a synthetic
  folded grid reproducing the reported artifact
  (plans/scripts/foldover_mechanism_demo.py — before = jagged oscillatory bands
  growing with salinity, after = clean).
- **y-branch rectilinearization** in PlanetProfileApp/pages/Exploreogram.py: the
  sigmaMean-as-Y substitution injected the RAW 2D array (the x-branch already used
  the column-mean fix) — now mirrored with a per-row mean, eliminating the
  non-rectilinear pcolormesh rows that produced the reported jagged y-axis at high
  salinity.
- **inductionData CWD fix** in PlanetProfile/MagneticInduction/MagneticInduction.py:
  GetBexc resolved `<Body>/inductionData` relative to the CWD, AND an empty
  auto-created PlanetProfileApp/Europa/inductionData shadowed the repo-root data —
  fallback now keys on the Be1xyz FILE, not the directory. Verified: excitation
  moments (11 Europa frequencies) load from the app CWD.

### MACHINE B — final verification run (user directive: intensive runs on B)

From repo root, latest genai:

```bash
mamba run -n PPcl python plans/scripts/exploreogram_jagged_verify.py
```

This drives the real Streamlit page (streamlit.testing AppTest) through a 5x5
Europa exploreogram with the user's default preset (sigmaMean/salinity 1-100 ppt x,
D_km/ice-thickness y, induction z, Bx, log axes), matplotlib path. ~25 Europa
models, minutes on Machine B. Then:
1. Locate the saved figure (script prints recent figure paths; also check
   Europa/figures/).
2. VISUALLY inspect: the y-axis must not show the jagged/oscillatory banding at
   high salinity; cells must not overlap.
3. Report verified / not-verified with the figure path via a handoff addendum.
(If the AppTest route mis-behaves on PPcl, the equivalent manual route is the
app's export-script feature or a direct ParPlanetExplore driver per the notes
above.)

## Superseded recon (for reference)

Recon already done (ranked suspects, file:line):

1. **PRIMARY: `ax.pcolormesh(x, y, z, shading='auto')` over non-monotone derived-axis
   coordinates** — PlanetProfile/Plotting/ExplorationPlots.py:371 (single-panel
   `PlotExploreOgram`, grid data from
   Utilities/EssentialHelpers.py:276 `extract_and_validate_plot_data`, reshaped at
   ExplorationPlots.py:321-324). When the y-axis is a DERIVED quantity (D_km driven
   by zb/Tb), the per-cell y-coordinates are model outputs which become non-monotone
   at high salinity → pcolormesh draws overlapping/distorted quads → wavy rows.
   No sorting is applied anywhere.
2. NaN'd broken models: runs set ALLOW_BROKEN_MODELS=True (Exploreogram.py:854);
   high salinity produces invalid models → NaN (ExplorationPlots.py:335-336) →
   ragged boundary interacting with (1).
3. Phase wrapping: Planet.Magnetic.phase used raw (ResultsIO.py:316,332), no
   np.unwrap — in amplitude_phase display mode ±180° wraps look oscillatory.
4. Plotly single-panel uses only excitation peak 0 (exploreogram_plotly.py:426) —
   multi-excitation selections silently show the first frequency only (separate,
   worth fixing while in there).

Suggested fix shape for (1): in the matplotlib path, when x/y come from
FigLbl.xCustomAxis/yCustomAxis (derived axes), sort each row/column by the axis
coordinate before pcolormesh (np.argsort per column for y, per row for x), carrying
z along; alternatively detect non-monotonicity and fall back to a regular grid in
the DRIVER variable with the derived value shown via contour labels. Keep the CLI
plotting behavior identical for non-derived axes.

**Verification requirement (house rule)**: reproduce first — run a small Europa
exploreogram (the new defaults; grid ~6x6 to keep runtime minutes) spanning high
salinity (x driver wOcean_ppt up to ~100 ppt), save the matplotlib PNG, LOOK at it,
then fix, re-render, and compare. The page's in-process entry is
`ParPlanetExplore(Planet, Params, xList, yList)` (Exploreogram.py:962) followed by
`GetLayerMeans`/`ExtractResults`; the export-script feature (Exploreogram.py:2097+)
writes a standalone matplotlib script that is convenient for headless reproduction.

NOTE: NaCl ocean conductivity model changed THIS session (Pan et al. 2021,
PlanetProfile/Thermodynamics/NaCl/NaClProps.py, now the NaCl default) — Europa's
default is Seawater (GSW, unchanged), so the jagged-plot repro is unaffected, but
any NaCl exploreograms will differ from older runs by design.

## Repo state snapshot (updated 2026-07-13, commit 3c55b4a6)

- Titan: 500k SBI artifact deployed, gates green within |Im k2| <= 0.20, GUI guard.
  1M clean-gate redeploy blocked on ratification items 1+2 (Machine A).
- Europa: **Galileo-run 1M SBI artifact DEPLOYED** (europa_seawater_andrade_posterior_1m.pt),
  synodic-only induction, with a GUI domain guard (|Im k2|<=0.15) + default Tb>=261.5 K
  truncation + scope note. GUI-verified via AppTest. Do not touch
  configs/europa_seawater_andrade_7D.json. See HANDOFF-2026-07-09 ADDENDUM 2026-07-13b.
- Callisto: Pan-2021 NaCl conductivity landed (bf545b5e); cache rebuild + Ae DROP diagnostic
  DONE; C-B phase (config + real k2/h2 targets) is Machine A's next authoring task.
- MgSO4: Pan et al. (2020) conductivity implemented (opt-in; Vance2018 still default) +
  regression test (Test/TestMgSO4Conductivity.py). verified.
- Known env quirk (Machine A): PPCallisto porous-silicate builds crash at
  PhaseConv phase ID 7 (ice VII pore fluid) — pre-existing, unrelated to the above.
- Exploreogram jagged y-axis plots (item 4 above): still **not implemented**.
