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

## Reviewer text fixes — re-verified 2026-10-20

| file | condition | observed |
|---|---|---|
| fix_default_sliders_hint.png | override on, default sliders (10,000 / derived 1000), default datum | "full-N likelihoods absent ... The datum is unchanged: set Posterior samples = 20,000 and Derived subset = All to use the precomputed correction"; no ".)." |
| fix_changed_datum.png | Im_k2 0.145 | "modified datum: x_obs differs ... The correction is validated only at the default Cassini datum"; no ".)." |
| fix_full_post.txt | 20,000 + All | corrected path still active; summary caption ends "all parameters weighted; the quarantine lift applies to the tidal set named in the warning." |

The container drops italic and inline-code glyphs ("Observable are free" is missing an italic word; "(see )" is blank), which supports the reviewer's font-artifact reading of item (d). Confirm once in a Mac browser.

## Flag flip — FEATURE_ENABLED = True, verified 2026-10-21

No env var set, i.e. the deployed default.

| file | condition | observed |
|---|---|---|
| flip_deployed_default_corrected.png | 20,000 + All, default datum | "Tidal sector at the DEFAULT Cassini datum only ..."; ESS 1,271; Pareto-k 0.19; max weight 0.62%; ocean 0.929 |
| flip_deployed_default_sliders_hint.png | default sliders | size-mismatch fallback with the "The datum is unchanged: set Posterior samples = 20,000 and Derived subset = All" hint |
| flip_env0_flag_off.png | PP_IS_CORRECTION_FEATURE=0 | original split-ratification warning, no corrected path |
