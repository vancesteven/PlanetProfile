# B8 flag-on visual verification — 2026-10-20

Status: **verified** (the reviewer cites it as satisfying B8; the flag-flip verdict is FLIP AFTER changes, see plans/STATUS.md).

Setup:
- HEAD 4956d58d plus the PP_IS_CORRECTION_FEATURE env override, which was committed later in 6c7e3f62.
- Linux aarch64 container, micromamba PPcl, headless Chromium.
- Icons render as boxes because the container has no emoji font.

| file | condition | observed |
|---|---|---|
| on_default_view5.png / on_default_post.txt | env override on, default Cassini datum, 20000 samples, derived "All" | corrected warning scoped to Re k2, \|Im k2\|, zeta, eta x5, alpha; ESS 1,271; Pareto-k 0.19; ocean 0.929; tension caveat |
| on_changed_view4.png / on_changed_post.txt | override on, Im_k2 0.135 -> 0.145 | fallback text, no corrected label |
| off_default_view5.png / off_default_post.txt | flag off (control) | original split-ratification warning |

At the GUI default sliders (10,000 samples / derived 1000), the app falls back via the B4 n_derived guard. Per the reviewer, this is correct and fail-closed; the message wording is to be fixed.
