"""Human-readable statement of the modeling assumptions behind an
inference run (user 2026-07-24: "every run needs to clearly describe
assumptions").

`describe_assumptions(result)` builds markdown from the run's own
config — sampled priors, fixed parameters, rheology, the viscosity
law actually applied in the tidal forward model, and the explicit
list of physics NOT modeled (porosity, convective shell partition,
silicate/core thermodynamics) — so the statement always matches the
run being displayed rather than a hand-maintained doc.
"""
from __future__ import annotations

from pathlib import Path


def _prior_line(name: str, spec: dict) -> str:
    kind = (spec or {}).get('prior_type', spec.get('type', '?')
            if isinstance(spec, dict) else '?')
    if isinstance(spec, dict) and 'bounds' in spec:
        lo, hi = spec['bounds']
        rng = f"[{lo:g}, {hi:g}]"
    elif isinstance(spec, dict) and 'mean' in spec:
        rng = f"mean {spec['mean']:g} ± {spec.get('std', float('nan')):g}"
    else:
        rng = ''
    return f"`{name}`: {kind} {rng}".rstrip()


def describe_assumptions(result) -> str:
    cfg = result.config
    param_space = getattr(cfg, 'param_space', {}) or {}
    fixed = getattr(cfg, 'fixed_params', {}) or {}
    groups = getattr(cfg, 'param_groups', {}) or {}
    derived = getattr(cfg, 'derived_params', {}) or {}
    obs = getattr(cfg, 'observables', {}) or {}
    ap = (getattr(cfg, 'arrhenius_params', None)
          or (getattr(cfg, 'sampler_settings', {}) or {}).get(
              'arrhenius_params') or {})

    names = set(param_space) | set(fixed)
    andrade = 'alpha' in names or any(n.startswith('log10_zeta')
                                      for n in names)
    maxwell_mu = 'log10_mu_Ih' in names

    lines = []
    body = getattr(cfg, 'bodyname', '') or 'body'
    mode = ('amortized neural posterior (SBI/NPE)'
            if (result.metadata or {}).get('sampler') == 'sbi'
            else 'MCMC (emcee ensemble sampler)')
    lines.append(f"**Inference:** {body}, {mode}. The posterior is over "
                 "the parameters below; everything else is held at the "
                 "values baked into the structure cache.")

    if param_space:
        lines.append("**Sampled parameters and priors:** "
                     + "; ".join(_prior_line(n, s)
                                 for n, s in param_space.items()) + ".")
    if groups:
        lines.append("**Grouped:** "
                     + "; ".join(f"`{k}` also sets "
                                 + ", ".join(f"`{m}`" for m in v)
                                 for k, v in groups.items()) + ".")
    if fixed:
        lines.append("**Fixed (not sampled):** "
                     + ", ".join(f"`{k}` = {v:g}"
                                 for k, v in fixed.items()) + ".")
    if derived:
        lines.append("**Derived parameters:** "
                     + "; ".join(f"`{k}` via {v.get('derivation', '?')}"
                                 for k, v in derived.items()) + ".")
    if obs:
        lines.append("**Conditioned on:** "
                     + ", ".join(f"{k} = {v[0]:.4g} ± {v[1]:.2g}"
                                 for k, v in obs.items()) + ".")

    cache = Path(str(getattr(cfg, 'structure_cache_path', '') or '')).name
    lines.append(
        "**Interior structure:** each draw's radial profile "
        "(ρ, T, P, elastic moduli, layer boundaries) is interpolated from "
        f"a precomputed PlanetProfile structure cache (`{cache}`) — "
        "hydrostatic, spherically symmetric, self-consistent profiles "
        "built by the full PlanetProfile forward model at grid nodes, "
        "blended per-layer between nodes. Only the sampled parameters "
        "vary; ocean composition, silicate/core EOS, and heat-flux "
        "assumptions are those of the cache build.")

    rheo = ('Andrade (sampled α, per-phase ζ)' if andrade
            else ('Maxwell with sampled ice-Ih shear modulus'
                  if maxwell_mu else 'Maxwell'))
    lines.append(
        f"**Tidal response (k₂, h₂, heating):** TidalPy radial solver on "
        f"the layered viscoelastic profile. Solid ice and silicate use "
        f"{rheo} rheology; the liquid ocean is a Newtonian fluid; the "
        "iron core and any clathrate layer are treated as elastic.")

    visc = ("**Viscosity:** the base η(r) profile comes from the "
            "structure cache; any sampled `log10_eta_*` replaces η "
            "uniformly within that phase's layers (Ih, III/V/VI or "
            "lumped HP, silicate).")
    if ap:
        Es = (ap.get('activation_energy_kJ_mol', {}) or {})
        e_txt = (", ".join(f"{k}: {v:g} kJ/mol" for k, v in Es.items())
                 or "config values")
        visc += (" An Arrhenius temperature dependence is then applied "
                 "layer-by-layer: η(T) = η_ref · exp[(E/R)(1/T − "
                 "1/T_ref)], with T_ref set to the draw's basal ice "
                 "temperature — so the sampled η is the basal (melting-"
                 f"point) viscosity. Activation energies: {e_txt}.")
    else:
        visc += (" No Arrhenius temperature scaling is configured for "
                 "this run — η is uniform within each phase layer.")
    lines.append(visc)

    lines.append(
        "**Moment of inertia:** C/MR² is the structure integral of the "
        "draw's density profile (core-aware when core parameters are "
        "sampled); gravity-pair configs derive the hydrostatic C₂₀/C₂₂ "
        "reference from it via Radau–Darwin.")

    lines.append(
        "**Not modeled:** porosity is OFF everywhere — structure caches "
        "are built with `Do.POROUS_ROCK = False`, so φ, Ppore, "
        "rhoMatrix/rhoPore (and thermal conductivity k, Q_S) have no "
        "values in the data table. The ice shell carries no "
        "conductive/convective partition (wedge drawn fully conductive). "
        "Silicate and core thermodynamics (Cp, α) are not in the cache; "
        "in the data table Cp and α for H₂O phases come from SeaFreeze, "
        "with the ocean approximated by the pure-water spline "
        "(few-percent error at seawater salinities).")

    return "\n\n".join(lines)
