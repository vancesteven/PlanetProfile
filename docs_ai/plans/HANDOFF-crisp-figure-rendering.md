# HANDOFF — Crisp figure rendering for the Inference page (Machine A)

**Owner to complete:** Machine A
**Author:** genai session, 2026-07-22
**File in scope:** `PlanetProfileApp/pages/Inference.py` only
**Baseline commit:** `ec711672` (branch `genai`)
**Status of the pilot:** the **corner plot** is done and `verified` (user confirmed in-browser
"looks good"). This handoff applies the *same* fix pattern to the remaining 11 `st.pyplot(...)`
figures on the page.

---

## Problem (same root cause as the corner plot)

`st.pyplot(fig)` rasterizes the matplotlib figure to a **screen-DPI PNG**, then the browser
downscales that bitmap into the ~700 px column. On large/many-subplot figures the text ends up
blurry and unreadable. Fix = render **inline as vector SVG** (crisp text at any zoom) and offer a
**PDF/PNG download**; for figures with tens of thousands of raw sample points, **rasterize only the
heavy 2D artists** so the SVG doesn't balloon to ~10 MB and choke the browser.

## The pattern to copy (already implemented — reuse, do NOT reinvent)

Two module-level helpers already exist in `PlanetProfileApp/pages/Inference.py`:

- **`_render_corner_figure(samples, labels, *, seed=0)`** — line 2187. Corner-specific; not reused
  here except as a reference for the rasterize-collections idiom.
- **`_display_corner(fig, *, key)`** — line 2228. **This is the reusable display path.** It:
  1. `fig.savefig(buf, format='svg', bbox_inches='tight')` then
     `st.image(buf.getvalue().decode('utf-8'), width='stretch')`.
     ⚠️ **The `.decode('utf-8')` is load-bearing** — `st.image` renders SVG only from a *string*;
     raw bytes are routed to PIL and raise `cannot identify image file <_io.BytesIO...>` (the exact
     bug we already hit and fixed).
  2. Adds side-by-side **Download PDF** + **Download PNG (200 dpi)** buttons, each with a unique
     `key=f'{key}_pdf'` / `key=f'{key}_png'`.

### Step 1 — generalize the display helper (rename, keep corner call working)

Promote `_display_corner` into a general `_display_vector_fig(fig, *, key, download_label='plot')`
(or add a thin general sibling and have `_display_corner` delegate). Same body; only the download
button label/file-name should be parameterized. Keep the existing corner call site working.

### Step 2 — add a rasterize helper for HEAVY figures

```python
def _rasterize_heavy_artists(fig):
    """Rasterize scatter/line-collection/image density so SVG/PDF stay small;
    text, axes, legends, colorbars stay vector."""
    for ax in fig.get_axes():
        for coll in ax.collections:   # scatter PathCollections, LineCollections
            coll.set_rasterized(True)
        for im in ax.images:
            im.set_rasterized(True)
```

Call it on HEAVY figures **before** `_display_vector_fig`. (This is exactly what
`_render_corner_figure` does at its tail.) A per-figure `fig.set_dpi(150)` before save keeps the
rasterized layer sharp without exploding size.

### Step 3 — replace each `st.pyplot(fig)` with the helper

```python
# LIGHT figure:
_display_vector_fig(fig, key='layer_thicknesses', download_label='layer thicknesses')

# HEAVY figure:
_rasterize_heavy_artists(fig)
_display_vector_fig(fig, key='k2_scatter', download_label='k2 posterior')
```

Keys **must be unique per figure**, and for the one looped site, unique per iteration.

---

## The 11 sites (line numbers at baseline `ec711672`; they will shift as you edit top-down —
## edit BOTTOM-UP or re-grep after each change)

### HEAVY — rasterize the 2D artists, then vector-display

| Line | var | What it shows | Heavy primitive | Key note |
|---|---|---|---|---|
| 2555 | `fig` | k₂ posterior: \|Im k₂\| vs Re k₂, all samples, colored by silicate-heating fraction; red 1σ/2σ ellipses | full-posterior `scatter` + colorbar | single |
| 2619 | `fig` | Aₑ complex plane, one subplot per excitation label; per-sample Re/Im scatter + \|Aₑ\|=1 circle | `scatter` per panel | single (Ae-units branch) |
| 2658 | `fig` | Bind = Aₑ·Be surface field, one subplot per component; scatter + conditioned star + 1σ circle | `scatter` per panel | **IN A LOOP** `for lab in labels:` → key must include `lab` |
| 2756 | `figk` | complex-plane signals by model (salinity/v3+): per-sample k₂/h₂/Aₑ paths, colored by salinity, sized by D_ocean | `LineCollection` (≤1200) + `scatter` + colorbar | single; shares name `figk` with 2824 — disambiguate key (e.g. `signals_salinity`) |
| 2824 | `figk` | complex-plane signals by model (fixed-salinity/v1–v2 else): node-grouped paths over faint per-sample k₂/h₂ clouds | faint `scatter` clouds + node `plot` | single; key `signals_fixed` |
| 3581 | `fig` | heating distribution: per-parameter scatter panels (heating vs each param, per phase) + stacked-bar panel | `scatter` across params×phases | single |

### LIGHT — vanilla vector display is enough (no rasterization needed)

| Line | var | What it shows | Primitive |
|---|---|---|---|
| 2470 | `fig` | interior layer-thickness histograms (ice/hydrosphere/mantle/core), 2-col grid, median + 16–84% lines | `hist`, `axvline` |
| 3135 | `figu` | non-hydrostatic identifiable combination u histogram, median + 68% + ±σ_u bands | `hist`, `axvline`, `axvspan` |
| 3316 | `fig` | C/MR² posterior: actual vs hydrostatic-reference histograms + observed Gaussian + 1σ/2σ bands | `hist`, `plot`, `axvline`, `axvspan` |
| 3474 | `fig_stack` | per-model reservoir fractions stacked bar (Silicate/Ice VI/V/III/Ih), one bar per model | `bar` (stacked) |
| 3495 | `fig_pie` | median heating partitioning pie (Silicate / HP Ice / Ice Ih), inside `col_pie`, `if tot_med>0` guard | `pie` |

**Watch item (3474):** ~5 layers × n_models rectangles. If n_models is large (many hundreds+),
the vector rect count grows — if a specific run makes the SVG heavy, rasterize this one too. No
scatter/hist2d, so LIGHT by default.

---

## Verification (per CLAUDE.md UI discipline — do NOT mark `done`/`fixed`)

A figure fix is **not** verified by `py_compile` or code review. Required:
1. Launch the app (`streamlit run PlanetProfileApp/…` or refresh a running instance on `:8501`).
2. Generate a posterior so each figure renders; **visually confirm** text is now crisp in-browser.
3. Click each **Download PDF**; open the PDF and confirm vector text is sharp and the plot is
   complete. Cite the artifact.
4. For HEAVY figures, confirm the inline SVG renders without browser stall and the exported SVG/PDF
   is small (target: single-MB SVG, sub-MB PDF — the corner plot went 10 MB→1.2 MB SVG / 0.6 MB PDF).
5. Status vocabulary: `verified` (with artifact) / `implemented, unverified` / `not implemented`.
   A figure whose in-browser render you have not eyeballed is `implemented, unverified`.

## Gotchas learned on the corner pilot
- `st.image` needs an **SVG string**, not bytes (`.decode('utf-8')`).
- Rasterize **collections/images only**, never the whole figure, or you lose the crisp text.
- corner's built-in Gaussian smoothing (`smooth=`, `smooth1d=`, `plot_datapoints=False`,
  `fill_contours=True`) was chosen over scipy `gaussian_kde` because KDE added ~8 s/plot for a
  marginal smoothness gain. For these non-corner scatter plots there is no smoothing decision —
  just rasterize the scatter.
- `st.download_button` supports `icon=':material/download:'` in the installed Streamlit (1.59.1).
- `width='stretch'` (NOT the deprecated `use_container_width`).
