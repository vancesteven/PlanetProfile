"""Query-time scope guard for the Enceladus ocean-branch SBI artifact.

Reviewer adjudication 2026-10-27 (validation_reports/enceladus_isostasy/sbi/
RESULT.md): the artifact is deployed SCOPE-BOUND. Preregistered gates G3 and
G5 FAILED; the thin-shell regime (zb < 15 km) is mis-shaped (over-dispersed zb,
mis-calibrated C2) as well as leaking into the cache hole. So:

  * REFUSAL tier: draw ``n_draws`` flow samples; if MORE than ``threshold``
    (preregistered 1 %) fall at zb < ``zb_thin_km`` (15 km) OR in a cache-hole
    cell, return NO posterior -- a ``Refused`` with the measured fractions and
    the instruction to use pocoMC.
  * Below the threshold: draws that fall in hole cells are DROPPED (support
    enforcement: the target pi_eff * L is identically zero there, the same
    logic as reject_outside_prior on the box) and the dropped fraction is
    ALWAYS reported. This does NOT change G5, which stays FAIL and must never
    be re-evaluated on truncated draws.
  * ``validated`` is True ONLY for x equal to the Cassini datum or one of the
    six G4 injection data sets; every other accepted x carries the
    "not individually validated" label.

A hole cell is any theta for which ``apply_parameters`` raises
``UnservableSampleError`` or the ocean (zb, w) struct is None -- the same two
tests, in the same order, as the likelihood (and
plans/scripts/enceladus_prior_pushforward.classify_and_derive).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from .forward_models import apply_parameters, UnservableSampleError

THRESHOLD = 0.01
ZB_THIN_KM = 15.0
N_DRAWS = 10000

SCOPE_TEXT = (
    'Enceladus ocean-branch SBI (scope-bound). Validated against pocoMC at the '
    'Cassini datum and on six off-nominal injection data sets. It is not '
    'globally validated: preregistered gates G3 (log10_w rank uniformity, '
    'diagnosed as held-out sampling fluctuation) and G5 (thin-shell cache-hole '
    'leakage) FAILED. Queries whose posterior puts more than 1 % mass at '
    'zb < 15 km or in a cache-hole cell are refused; use pocoMC for those. '
    'Below that threshold, draws falling in hole cells are dropped and the '
    'fraction is shown. Results at other x are "effective-prior posterior; '
    'servable cells only; not individually validated". Cite pocoMC run (a), '
    'not this artifact, for published values.')
NOT_VALIDATED_LABEL = ('effective-prior posterior; servable cells only; '
                       'not individually validated')


@dataclass
class Refused:
    reason: str
    frac_thin: float
    frac_hole: float
    frac_union: float
    n_draws: int
    scope_text: str = SCOPE_TEXT

    def __bool__(self):
        return False


@dataclass
class Accepted:
    samples: np.ndarray
    dropped_fraction: float
    frac_thin: float
    frac_hole: float
    n_draws: int
    validated: bool
    label: str
    scope_text: str = SCOPE_TEXT
    meta: Dict = field(default_factory=dict)

    def __bool__(self):
        return True


def is_hole(runner, theta) -> bool:
    td = runner._expand_theta(np.asarray(theta, float))
    try:
        apply_parameters(td, runner.structure_data)
    except UnservableSampleError:
        return True
    if 'zb_km' in td and runner._ocean_zbw_struct_for_theta(td) is None:
        return True
    return False


def _same_x(x: Dict[str, float], ref: Dict[str, float], rtol=1e-12) -> bool:
    return (set(x) == set(ref)
            and all(np.isclose(float(x[k]), float(ref[k]), rtol=rtol, atol=0.0)
                    for k in ref))


def guarded_posterior(sbi_runner, mcmc_runner, x_obs: Dict[str, float], *,
                      validated_x: Sequence[Dict[str, float]] = (),
                      n_draws: int = N_DRAWS, seed: int = 75,
                      threshold: float = THRESHOLD,
                      zb_thin_km: float = ZB_THIN_KM):
    """Return ``Accepted`` or ``Refused`` (falsy) for one conditioning ``x_obs``.

    ``sbi_runner``: an SBIRunner with the artifact loaded (``sample_posterior``).
    ``mcmc_runner``: an MCMCRunner on the same config + production cache, used
    only for the hole classification.
    ``validated_x``: the x dicts the artifact was validated at (the Cassini
    datum and the six G4 injection data sets).
    """
    th = np.asarray(sbi_runner.sample_posterior(
        x_obs, n_samples=n_draws, seed=seed, reject_outside_prior=True), float)
    names = list(sbi_runner.param_names)
    zb = th[:, names.index('zb_km')]
    thin = zb < zb_thin_km
    hole = np.array([is_hole(mcmc_runner, t) for t in th], bool)
    union = thin | hole
    f_thin, f_hole, f_union = (float(thin.mean()), float(hole.mean()),
                               float(union.mean()))
    if f_union > threshold:
        return Refused(
            reason=(f'{100 * f_union:.2f} % of {len(th)} flow draws fall at '
                    f'zb < {zb_thin_km:g} km or in a cache-hole cell '
                    f'(threshold {100 * threshold:g} %): the thin-shell regime '
                    'is outside this artifact\'s validated scope -- use pocoMC.'),
            frac_thin=f_thin, frac_hole=f_hole, frac_union=f_union,
            n_draws=len(th))
    validated = any(_same_x(x_obs, v) for v in validated_x)
    return Accepted(
        samples=th[~hole], dropped_fraction=f_hole, frac_thin=f_thin,
        frac_hole=f_hole, n_draws=len(th), validated=validated,
        label=('validated against pocoMC at this datum (scope-bound artifact)'
               if validated else NOT_VALIDATED_LABEL),
        meta={'threshold': threshold, 'zb_thin_km': zb_thin_km, 'seed': seed})
