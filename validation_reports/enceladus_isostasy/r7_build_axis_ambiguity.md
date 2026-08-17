# r7 build blocked: zb-axis endpoint convention is self-contradictory — no compute spent

Before launching the ~19.5h production build authorized by r7
(`r7_ADJUDICATION.md`, §0.29), I reconstructed the zb axis from
`enceladus_cassini_isostasy_7D.json`'s `metadata.structure_cache_spec
.zb_km_grid` to verify it against the declared `n_total: 87` before
spending any compute. It does not reconcile, and the discrepancy is
exactly one node at the 45 km ocean/no-ocean support edge — the same
edge the config's own commentary flags as safety-critical.

## The two readings

The five segments are declared as:

| lo | hi | step | declared n |
|----|----|------|-----------:|
| 5.0 | 20.0 | 1.0 | 15 |
| 20.0 | 22.0 | 0.5 | 4 |
| 22.0 | 30.0 | 0.25 | 32 |
| 30.0 | 42.0 | 0.5 | 24 |
| 42.0 | 45.0 | 0.25 | 12 |

Sum of declared `n` = 87, matching `zb_km_grid.n_total: 87`.

**Reading A (half-open per segment, shared boundary counted once as
each segment's `lo` only):** `[5.0, 6.0, ..., 19.0]` (15) +
`[20.0, 20.5, 21.5]`... i.e. each segment contributes floor((hi-lo)/step)
nodes starting at its own `lo`, and the final segment's `hi=45.0` is
NEVER emitted because it would be the 13th node of a 0.25 km/3 km
segment (12 is `(45-42)/0.25 = 12`, i.e. floor, not floor+1). This
reading reproduces 87 nodes exactly — but the last node is 44.75 km,
**not 45.0 km**. Verified directly:

```
segs = [(5.0,1.0,15),(20.0,0.5,4),(22.0,0.25,32),(30.0,0.5,24),(42.0,0.25,12)]
zb = [lo+step*i for lo,step,n in segs for i in range(n)]
len(zb) -> 87; zb[-1] -> 44.75
```

**Reading B (inclusive of both endpoints, per the config's own
`endpoint_convention` field):**

> "Segments are inclusive of BOTH endpoints... This fixes reviewer
> NEW-MODERATE-1: the old half-open reading left the box edge with no
> node, so 1.5% of the zb prior had no bracketing pair for the
> bilinear interpolant."

Under this reading each segment is `linspace(lo, hi, round((hi-lo)/step)+1)`
and adjacent segments share their boundary node (counted once). The
last segment (42→45, step 0.25) has `(45-42)/0.25+1 = 13` nodes, not
12. Full inclusive count across all five segments (de-duplicating the
4 shared boundaries at 20, 22, 30, 42): **88 nodes**, span
`[5.0, 45.0]` km — the 45.0 km edge IS present.

## Why this is not a rounding nitpick

The 42-45 km segment's own stated rationale is:

> "OCEAN/NO-OCEAN SUPPORT EDGE at ~45 km. THIS is where the Titan
> gate-3 absolute-onset-placement lesson applies... The edge is
> w-DEPENDENT..., so it is a curve in (zb,w) and the fine step must
> hold across the w axis."

Reading A — the one that reproduces the declared `n_total: 87` — has
**no node at 45.0 km at all**. Every (zb, w) sample whose true zb lands
between 44.75 and 45.0 would be extrapolated past the last grid node
rather than bilinearly interpolated between bracketing nodes, at
precisely the edge the segment exists to resolve. This is the same
defect class E1 was about (a placement/bracketing mismatch at a
physically load-bearing boundary), just one node over rather than a
root-find tolerance bug.

The config's own `endpoint_convention` text explicitly names this
exact failure mode ("the old half-open reading left the box edge with
no node") as something already fixed. Reading A silently reintroduces
it; the declared `n_total: 87` is consistent ONLY with Reading A.
**The document is self-contradictory: the endpoint-convention prose
promises Reading B's behavior, but the segment `n` values (which sum
to exactly 87, matching the declared total) are Reading A's.**

## What I did NOT do

I did not pick a reading and build. Both are defensible mechanically,
they disagree by exactly one node at the safety-critical edge, and r7
represents 44 tests + an independent 12-node production verification
already invested in this axis — not a shape I should silently resolve
by picking whichever arithmetic is more convenient. No compute spent;
`plans/scripts/enceladus_zbw_production_build.py --dry-run` catches
this before any `build_single_structure` call and refuses to proceed.

## Decision needed

Which reading is correct — 87 nodes ending at 44.75 km (Reading A,
matches the declared total as-is), or 88 nodes ending at 45.0 km
(Reading B, matches the endpoint-convention prose, requires
correcting `n_total` to 88 and the last segment's `n` to 13)? Once
ruled, `plans/scripts/enceladus_zbw_production_build.py` needs at most
a one-line change to `_zb_axis()` and `--dry-run` re-verifies before
the real ~19.5h build launches.
