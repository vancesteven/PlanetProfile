# Direct answer to §0.25's key question: does the failing region intersect the GUI-reachable box?

§0.25 (`plans/MACHINE-B-HANDOFF.md`) asks: "does the failing region intersect
the GUI-reachable x_obs_limits box, or lie beyond it? If beyond, deploy
posture is UNCHANGED and the finding is a documented boundary, not a new
quarantine." It states the expected disposition as scope-bound: "corrected
at the fiducial, uncorrected-and-quarantined elsewhere."

**Answer, computed directly from the already-collected c12_sweep_report.json
(no new compute): the failing region DOES intersect the GUI-reachable box.**

The NH3 slot's only declared `x_obs_limits` entry is `Im_k2 in (0.0, 0.30)`
(`PlanetProfileApp/pages/Inference.py`, `titan_freegrav_nh3_posterior_1m.pt`
slot). C20, C22, and Re_k2 have no declared GUI box, so "GUI-reachable" is
defined here as `Im_k2 <= 0.30` alone (the widest reading the deploy code
actually enforces), then separately intersected with the training-support
box for a stricter reading.

| region | n | failed | frac_fail |
|---|---|---|---|
| All 149 swept points | 149 | 29 | 19.5% |
| Im_k2 <= 0.30 (the ONLY declared GUI limit) | 114 | 17 | **14.9%** |
| Training-support box (Re_k2 in [-0.1,1.5], Im_k2 in [0,1.0]) | 123 | 26 | 21.1% |
| Training-support box AND Im_k2 <= 0.30 (narrowest reachable set) | 113 | 16 | **14.2%** |

Even under the narrowest reading of "reachable" — a point a GUI user could
actually submit AND that sits inside the flow's own training support — 16 of
113 points (14.2%) still fail Pareto-k <= 0.7. This is not a boundary effect
concentrated just outside the box; it persists well inside it.

## What the in-box failures look like

Of the 17 failures with Im_k2 <= 0.30: 1 is a hard C4 reject (75.6%
sentinel-rejected, at the fiducial's own C22 but Re_k2 pushed to the
training-box edge 1.501 and Im_k2=0.048 — i.e. a point the box formally
admits but the guard cannot actually condition on), 4 are "wild" (Pareto-k >
1.5, up to 16.27), and 12 are "marginal" (Pareto-k in 0.7-1.5). So it is not
purely a hard-edge failure mode — most in-box failures are moderate
degradations of reliability, not catastrophic ones, but they are common
enough (1 in 7) that "the correction is validated within the reachable
domain" is not supportable as stated.

## What this does NOT settle

- Whether these failures matter in practice depends on how often actual GUI
  users condition near these regions versus near the fiducial — this sweep
  draws from the PRIOR-predictive distribution, which need not match the
  distribution of x's real users submit.
- This does not re-open the fiducial-point validations (C3/C5.3/C10/C11/
  C13/C16), which are unaffected.
- Re_k2 has no declared GUI limit at all, so any "GUI-reachable" claim about
  it is a description of the training-support box only, not an enforced
  deploy-time guard — worth the reviewer/manager flagging if a Re_k2 box is
  intended to exist and is simply undeclared (parallel to r5's D1 finding
  about an unenforced ocean-branch MoI window on the Enceladus campaign).

Not self-adjudicated. This resolves the FACTUAL question posed in §0.25;
the DISPOSITION (does 14-15% failure inside the reachable box change deploy
posture, warrant a narrower `x_obs_limits`, or get accepted as a documented
caveat) is the reviewer's/manager's call, not mine.
