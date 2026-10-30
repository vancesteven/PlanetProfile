"""Two-run branch comparison by post-hoc Bayes factor (task A6).

Frozen-branch design ruling (validation_reports/enceladus_isostasy/
frozen_branch_DESIGN_RULING.md): the ocean/no-ocean branch is drawn by TWO
SEPARATE nested-sampling runs, each integrating over ITS OWN support, and
combined afterwards at declared 0.5/0.5 prior odds. A single-theta branch
indicator was REJECTED.

Why the discrete coordinate was rejected
----------------------------------------
The two branches are NOT NESTED. The ocean branch samples five parameters
(zb_km, log10_wOcean_ppt, rho_ice_kgm3, compensation_C2, libration_sys_frac);
the frozen branch samples three (rho_rock_kgm3, rho_ice_kgm3,
libration_sys_frac), with zb DERIVED by mass conservation and both
log10_wOcean_ppt and compensation_C2 undefined rather than merely
unconstrained. A discrete model index inside one sampler would therefore need
a PSEUDO-PRIOR over the coordinates that exist only on one side -- and the
reported branch mass would depend on that arbitrary choice. That is exactly
the fragile statistic the NH3 C16 work identified, relocated rather than
cured. Two independent evidence integrals have no such freedom.

What this module does and does not do
-------------------------------------
It does the ARITHMETIC and the VALIDITY CHECKS, not the sampling. The two
nested-sampling runs are compute (Machine B); each writes ``log_Z`` and
``log_Z_err`` into its result metrics (``mcmc_runner`` already records both
from the sampler). This module combines them and refuses combinations that
are not comparable.

The single most important guardrail: a Bayes factor between two runs is only
meaningful if BOTH ran on IDENTICAL DATA with an IDENTICALLY NORMALIZED
likelihood. Two runs that differ in their observable set, their sigmas, or
their Sigma_model inflation produce log_Z values whose difference contains
the ratio of two different normalizing constants, and the resulting "Bayes
factor" is an artifact. :func:`branch_bayes_factor` requires the caller to
pass both observable dicts and compares them exactly.
"""
from __future__ import annotations

import math
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

__all__ = ['branch_bayes_factor', 'branch_comparison_from_results']

# Declared prior odds for the branch comparison (ruling; config
# metadata.branch_model.prior_odds). Equal odds is a DECLARATION, not a
# measurement -- it is what makes the reported posterior odds readable as the
# data's contribution alone.
DEFAULT_PRIOR_ODDS: Dict[str, float] = {'frozen': 0.5, 'ocean': 0.5}

# Above this |log B| the posterior odds underflow double precision and the
# NUMBER stops being interpretable even though the VERDICT is unambiguous.
# The frozen branch is discriminated at tens of sigma, so this is the
# expected regime and must be reported as "one-sided", never as a ratio.
_ONE_SIDED_LOG_B = 500.0


def _observables_identical(a: Mapping[str, Any],
                           b: Mapping[str, Any]) -> Tuple[bool, str]:
    """Exact comparison of two observable dicts {name: (value, sigma)}."""
    ka, kb = set(a), set(b)
    if ka != kb:
        return False, (f'observable NAMES differ: only in frozen '
                       f'{sorted(ka - kb)}, only in ocean {sorted(kb - ka)}')
    for k in sorted(ka):
        va, vb = a[k], b[k]
        try:
            pa = tuple(float(x) for x in va)
            pb = tuple(float(x) for x in vb)
        except (TypeError, ValueError):
            if va != vb:
                return False, f'observable {k!r} differs: {va!r} vs {vb!r}'
            continue
        if pa != pb:
            return False, (f'observable {k!r} differs: {pa} vs {pb} -- the '
                           'two evidences carry different likelihood '
                           'normalizations and cannot be differenced')
    return True, ''


def branch_bayes_factor(
    log_Z_frozen: float,
    log_Z_ocean: float,
    observables_frozen: Mapping[str, Any],
    observables_ocean: Mapping[str, Any],
    log_Z_frozen_err: Optional[float] = None,
    log_Z_ocean_err: Optional[float] = None,
    prior_odds: Optional[Mapping[str, float]] = None,
) -> Dict[str, Any]:
    """Post-hoc branch comparison from two independent evidence integrals.

    ``log B = log Z_frozen - log Z_ocean`` (ruling's sign convention), then
    ``log posterior odds = log B + log(pi_frozen / pi_ocean)``.

    Args:
        log_Z_frozen, log_Z_ocean: log evidence from each branch's own
            nested-sampling run, each integrated over ITS OWN support.
        observables_frozen, observables_ocean: the conditioned observable
            dicts ``{name: (value, sigma)}`` the two runs used. These MUST
            be identical; see the module docstring.
        log_Z_frozen_err, log_Z_ocean_err: sampler-reported numerical
            uncertainties. Optional, but their absence is recorded in the
            result so a reader can see the quoted log B is unbudgeted.
        prior_odds: declared prior model probabilities, default 0.5/0.5.
            Need not be normalized; only the ratio is used.

    Returns a dict carrying ``log_bayes_factor``, ``log_bayes_factor_err``
    (or None), ``log_posterior_odds``, ``P_frozen`` / ``P_ocean``,
    ``one_sided`` (True when the odds underflow and only the direction is
    reportable), ``prior_odds``, ``favoured`` and ``notes``.

    Raises:
        ValueError: if either evidence is missing/non-finite, if the two
            observable sets are not identical, or if a prior odd is
            non-positive.
    """
    ok, why = _observables_identical(observables_frozen, observables_ocean)
    if not ok:
        raise ValueError(
            'branch evidences are not comparable -- ' + why + '. A Bayes '
            'factor requires both runs to condition on IDENTICAL data with '
            'an identically normalized likelihood.')

    for name, val in (('log_Z_frozen', log_Z_frozen),
                      ('log_Z_ocean', log_Z_ocean)):
        if val is None or not math.isfinite(float(val)):
            raise ValueError(
                f'{name} is {val!r}; both branches must report a finite log '
                'evidence. A branch whose sampler did not produce one has '
                'not been compared -- do not substitute a bound.')

    odds = dict(prior_odds) if prior_odds else dict(DEFAULT_PRIOR_ODDS)
    pi_f = float(odds.get('frozen', 0.0))
    pi_o = float(odds.get('ocean', 0.0))
    if pi_f <= 0.0 or pi_o <= 0.0:
        raise ValueError(
            f'prior odds must both be positive, got frozen={pi_f}, '
            f'ocean={pi_o}. A zero prior odd removes the branch from the '
            'comparison rather than letting the data speak.')

    log_B = float(log_Z_frozen) - float(log_Z_ocean)
    log_B_err = None
    if log_Z_frozen_err is not None and log_Z_ocean_err is not None:
        ef, eo = float(log_Z_frozen_err), float(log_Z_ocean_err)
        if math.isfinite(ef) and math.isfinite(eo):
            # Independent runs, so the numerical errors add in quadrature.
            log_B_err = math.sqrt(ef ** 2 + eo ** 2)

    log_post_odds = log_B + math.log(pi_f / pi_o)
    one_sided = abs(log_post_odds) > _ONE_SIDED_LOG_B
    if one_sided:
        P_frozen = 1.0 if log_post_odds > 0 else 0.0
        P_ocean = 1.0 - P_frozen
    else:
        # P_frozen = 1 / (1 + exp(-log_post_odds)), computed stably.
        if log_post_odds >= 0:
            P_frozen = 1.0 / (1.0 + math.exp(-log_post_odds))
        else:
            e = math.exp(log_post_odds)
            P_frozen = e / (1.0 + e)
        P_ocean = 1.0 - P_frozen

    notes = []
    if log_B_err is None:
        notes.append(
            'No sampler log_Z uncertainty was supplied for at least one '
            'branch; the quoted log B carries NO numerical error budget.')
    elif abs(log_B) <= 3.0 * log_B_err:
        notes.append(
            f'|log B| = {abs(log_B):.3f} is within 3x its numerical error '
            f'{log_B_err:.3f}: the comparison is NOT resolved by these runs. '
            'Increase n_effective rather than quoting the sign.')
    if one_sided:
        notes.append(
            f'|log posterior odds| = {abs(log_post_odds):.1f} exceeds '
            f'{_ONE_SIDED_LOG_B:.0f}; the odds underflow double precision. '
            'Report the verdict as ONE-SIDED, not as a ratio or a '
            'probability -- the magnitude is not interpretable.')
    notes.append(
        'Each evidence is integrated over its OWN support (frozen 3 free '
        'parameters, ocean 5). The branches are not nested and no discrete '
        'theta coordinate exists; the prior odds below are DECLARED, not '
        'inferred.')

    return {
        'log_bayes_factor': log_B,
        'log_bayes_factor_err': log_B_err,
        'log_bayes_factor_convention': 'log Z_frozen - log Z_ocean',
        'log_posterior_odds': log_post_odds,
        'P_frozen': P_frozen,
        'P_ocean': P_ocean,
        'one_sided': one_sided,
        'prior_odds': {'frozen': pi_f, 'ocean': pi_o},
        'log_Z': {'frozen': float(log_Z_frozen),
                  'ocean': float(log_Z_ocean)},
        'log_Z_err': {'frozen': log_Z_frozen_err,
                      'ocean': log_Z_ocean_err},
        'favoured': ('frozen' if log_post_odds > 0
                     else 'ocean' if log_post_odds < 0 else 'neither'),
        'notes': notes,
    }


def _extract(result: Any, key: str, default=None):
    """Read ``key`` from a runner result object, its ``.metrics``, or a dict."""
    for holder in (result, getattr(result, 'metrics', None)):
        if holder is None:
            continue
        if isinstance(holder, Mapping) and key in holder:
            return holder[key]
        if hasattr(holder, key):
            return getattr(holder, key)
    return default


def branch_comparison_from_results(
    frozen_result: Any,
    ocean_result: Any,
    observables_frozen: Optional[Mapping[str, Any]] = None,
    observables_ocean: Optional[Mapping[str, Any]] = None,
    prior_odds: Optional[Mapping[str, float]] = None,
) -> Dict[str, Any]:
    """:func:`branch_bayes_factor` driven straight off two runner results.

    Accepts anything that exposes ``log_Z`` / ``log_Z_err`` (and optionally
    ``observables``) directly, under ``.metrics``, or as mapping keys --
    which covers ``InferenceResult``, its ``metrics`` dict, and a JSON
    manifest reloaded from disk.

    ``observables_*`` override what is found on the results; pass them
    explicitly when the results do not carry their observable dicts, since
    the identity check is what makes the comparison legitimate.
    """
    obs_f = (observables_frozen if observables_frozen is not None
             else _extract(frozen_result, 'observables'))
    obs_o = (observables_ocean if observables_ocean is not None
             else _extract(ocean_result, 'observables'))
    if obs_f is None or obs_o is None:
        raise ValueError(
            'observable dicts for BOTH branches are required: the Bayes '
            'factor is only defined when the two runs conditioned on the '
            'same data. Pass observables_frozen / observables_ocean '
            'explicitly when the results do not carry them.')
    return branch_bayes_factor(
        log_Z_frozen=_extract(frozen_result, 'log_Z'),
        log_Z_ocean=_extract(ocean_result, 'log_Z'),
        observables_frozen=obs_f,
        observables_ocean=obs_o,
        log_Z_frozen_err=_extract(frozen_result, 'log_Z_err'),
        log_Z_ocean_err=_extract(ocean_result, 'log_Z_err'),
        prior_odds=prior_odds,
    )
