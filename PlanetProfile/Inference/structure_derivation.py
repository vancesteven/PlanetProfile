"""
Structure-derivation helpers for v2 BodyConfigs (Phase C1 Stage 2).

Provides three pure-math helpers used by the MCMC runner when a config declares
``derived_params.rho_sil_kgm3.derivation == 'mass_conservation'``:

  - ``compute_cmr2(layers, R_body_m, M_body_kg)``
        Exact analytic moment-of-inertia coefficient C/MR² for a body modeled
        as a stack of piecewise-constant-density spherical shells.

  - ``derive_silicate_density(...)``
        Solve for silicate density given total mass, hydrosphere mass, the
        silicate-shell geometry, and a sampled core. Returns a hard-bound
        rejection flag.

  - ``extract_hydrosphere_layers(structure)``
        Convert a PP-style cache structure dict into ``(layers, R_oceanbot_m,
        M_hydrosphere_kg)``, masking out cells with PP phase >= 50 (rocky /
        metallic interior).

Conventions
-----------
All quantities are SI (m, kg, kg/m³). A ``layer`` is a ``(r_inner_m,
r_outer_m, rho_kgm3)`` tuple, non-overlapping and with ``r_outer > r_inner >=
0``. PP phase codes (see ``PlanetProfile.Utilities.defineStructs``):
0=ocean, 1–7=ice phases, 30=clathrate, 50–59=silicate, 100/105=Fe,
110/115=FeS. The hydrosphere mask is ``phase < 50``.
"""
from typing import Iterable, List, Optional, Tuple, Dict, Any

import numpy as np


# Constants.phaseSil = 50 in PP. Anything at or above this is rocky/metallic.
PHASE_SIL_THRESHOLD: int = 50


def _shell_mass_inertia(r_in: float, r_out: float, rho: float) -> Tuple[float, float]:
    """Mass and polar moment-of-inertia of a single constant-density shell.

    Closed form:
        m_shell = (4π/3) (r_out³ - r_in³) ρ
        I_shell = (8π/15) (r_out⁵ - r_in⁵) ρ
    """
    M = (4.0 / 3.0) * np.pi * (r_out ** 3 - r_in ** 3) * rho
    I = (8.0 / 15.0) * np.pi * (r_out ** 5 - r_in ** 5) * rho
    return M, I


def compute_cmr2(
    layers: Iterable[Tuple[float, float, float]],
    R_body_m: Optional[float] = None,
    M_body_kg: Optional[float] = None,
) -> float:
    """Compute CMR² = I / (M · R²) for piecewise-constant spherical shells.

    Parameters
    ----------
    layers
        Iterable of ``(r_inner_m, r_outer_m, rho_kgm3)`` tuples. Must be
        non-overlapping with ``r_outer > r_inner >= 0`` per layer.
    R_body_m
        Surface radius. If ``None``, defaults to ``max(r_outer)`` across
        layers. Pass explicitly when the body's surface lies above the
        outermost cached shell.
    M_body_kg
        Total mass. If ``None``, defaults to the integrated mass of the
        provided layers. Pass explicitly when computing CMR² for a body with
        a known observed mass that is being held fixed by mass-conservation
        elsewhere.

    Returns
    -------
    cmr2 : float (dimensionless)
    """
    M_int = 0.0
    I_int = 0.0
    R_max = 0.0
    for r_in, r_out, rho in layers:
        if not (r_out > r_in >= 0.0):
            raise ValueError(
                f"Bad layer geometry: r_inner={r_in!r}, r_outer={r_out!r}"
            )
        m, i = _shell_mass_inertia(r_in, r_out, rho)
        M_int += m
        I_int += i
        if r_out > R_max:
            R_max = r_out
    R = R_body_m if R_body_m is not None else R_max
    M = M_body_kg if M_body_kg is not None else M_int
    if M <= 0.0 or R <= 0.0:
        raise ValueError(f"Invalid mass/radius: M_body_kg={M!r}, R_body_m={R!r}")
    return I_int / (M * R * R)


def derive_silicate_density(
    M_total_kg: float,
    M_hydrosphere_kg: float,
    R_oceanbot_m: float,
    R_core_m: float,
    rho_core_kgm3: float,
    bounds: Tuple[float, float] = (2200.0, 3500.0),
) -> Tuple[float, bool]:
    """Solve for silicate density via mass conservation, with hard-bound check.

        rho_sil = (M_total - M_hydrosphere - M_core) / V_silicate

    where ::

        M_core   = (4π/3) R_core³ ρ_core
        V_sil    = (4π/3) (R_oceanbot³ - R_core³)

    Parameters
    ----------
    M_total_kg
        Body's observed total mass (held fixed; not sampled).
    M_hydrosphere_kg
        Mass of the cached hydrosphere (ice + ocean + clathrate). Comes from
        ``extract_hydrosphere_layers``.
    R_oceanbot_m
        Inner radius of the hydrosphere = outer radius of the silicate shell.
    R_core_m
        Sampled core radius. ``0.0`` admits the no-core limit.
    rho_core_kgm3
        Sampled core density (typically in [5000, 8000] kg/m³ for metallic
        cores; the bound is enforced by the prior, not here).
    bounds
        ``(rho_min, rho_max)`` in kg/m³ for ρ_sil. Sample is rejected if the
        derived value falls outside.

    Returns
    -------
    (rho_sil_kgm3, accepted)
        ``rho_sil_kgm3`` is the derived density (NaN on degenerate geometry).
        ``accepted`` is True iff ``bounds[0] <= rho_sil <= bounds[1]`` AND
        the silicate volume is strictly positive.
    """
    V_sil = (4.0 / 3.0) * np.pi * (R_oceanbot_m ** 3 - R_core_m ** 3)
    if V_sil <= 0.0:
        return float("nan"), False
    M_core = (4.0 / 3.0) * np.pi * (R_core_m ** 3) * rho_core_kgm3
    M_sil = M_total_kg - M_hydrosphere_kg - M_core
    rho_sil = M_sil / V_sil
    rho_min, rho_max = bounds
    accepted = bool(rho_min <= rho_sil <= rho_max)
    return float(rho_sil), accepted


def extract_hydrosphere_layers(
    structure: Dict[str, Any],
) -> Tuple[List[Tuple[float, float, float]], float, float]:
    """Extract hydrosphere shells from a PP-style cache structure dict.

    The cache layout follows ``Test50``'s ``_build_structure`` output:
    arrays ``r_m``, ``rho``, ``phases`` are per-cell, with ``r_m[i]`` the
    OUTER face of cell ``i``. The bottom of cell ``i`` is ``r_m[i-1]``, and
    the bottom of cell 0 is treated as 0 m.

    Parameters
    ----------
    structure
        Dict with keys ``r_m``, ``rho``, ``phases`` (numpy arrays, same
        length).

    Returns
    -------
    (layers, R_oceanbot_m, M_hydrosphere_kg)
        ``layers`` is a list of ``(r_in, r_out, rho)`` tuples for cells with
        phase < 50. ``R_oceanbot_m`` is the inner radius of the lowest
        hydrosphere cell. ``M_hydrosphere_kg`` is the integrated mass.

    Notes
    -----
    Returns ``([], 0.0, 0.0)`` if no hydrosphere cells are present (caller
    is expected to handle this; mass-conservation derivation will reject).

    The cache is assumed to have contiguous hydrosphere cells (no rocky
    interruptions in the middle of the water column). PP structures always
    satisfy this — a brief sanity check is included.
    """
    r_m = np.asarray(structure["r_m"], dtype=np.float64)
    rho = np.asarray(structure["rho"], dtype=np.float64)
    phases = np.asarray(structure["phases"])

    if r_m.size != rho.size or r_m.size != phases.size:
        raise ValueError(
            "structure['r_m'], 'rho', 'phases' must all be the same length"
        )

    # Normalize to ascending r order so cell i runs from r_m[i-1] to r_m[i].
    sort_idx = np.argsort(r_m)
    r_m = r_m[sort_idx]
    rho = rho[sort_idx]
    phases = phases[sort_idx]

    r_inner = np.concatenate(([0.0], r_m[:-1]))
    r_outer = r_m

    hydro_mask = phases < PHASE_SIL_THRESHOLD
    if not np.any(hydro_mask):
        return [], 0.0, 0.0

    # Sanity: hydrosphere must be a contiguous block at the top of the body
    # (no rocky interruptions). Find the first hydrosphere index from below
    # and assert all cells above it are also hydrosphere.
    first_hydro = int(np.argmax(hydro_mask))
    if not np.all(hydro_mask[first_hydro:]):
        raise ValueError(
            "Hydrosphere cells (phase < 50) are not contiguous at the top "
            "of the cached structure. PP caches are expected to satisfy "
            "this; got phase array with non-contiguous water column."
        )

    R_oceanbot_m = float(r_inner[first_hydro])

    # Skip zero-thickness shells. PP caches legitimately duplicate radii at
    # layer interfaces (TidalPy's radial_solver requires the boundary radius
    # to be present twice — once for the layer below, once for the layer
    # above). Those duplicates manifest here as cells with r_in == r_out,
    # which compute_cmr2 rejects with "Bad layer geometry". A zero-thickness
    # shell contributes zero mass and zero inertia, so dropping it is the
    # correct behaviour for the mass-conservation derivation. Surfaced as
    # 0/4096 non-rejected samples for Callisto smoke MCMC on 2026-05-23.
    layers: List[Tuple[float, float, float]] = []
    M_hydro = 0.0
    for ri, ro, rh in zip(
        r_inner[first_hydro:], r_outer[first_hydro:], rho[first_hydro:]
    ):
        if ro <= ri:
            continue
        layers.append((float(ri), float(ro), float(rh)))
        m, _ = _shell_mass_inertia(ri, ro, rh)
        M_hydro += m
    return layers, R_oceanbot_m, float(M_hydro)
