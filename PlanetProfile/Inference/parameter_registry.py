"""
Parameter registry for Bayesian inference.

Centralizes parameter definitions for the inference UI. Adding new parameters
requires only:
1. Adding entry to PARAMETER_REGISTRY below
2. Implementing @parameter_hook in forward_models.py

No UI code changes needed - the registry drives dynamic UI generation.

Author: PlanetProfile Team
Date: 2026-04-29
"""
from dataclasses import dataclass
from typing import List, Optional, Dict, Any


@dataclass
class ParameterDef:
    """
    Definition of a single inference parameter.

    Attributes:
        id: Unique identifier (e.g., 'alpha', 'log10_eta_Ih')
        label: Human-readable name for UI
        latex_label: LaTeX string for plotting (e.g., r'$\alpha$')
        description: Help text explaining the parameter
        category: Grouping category ('rheology', 'structure', 'ocean', 'magnetic')
        default_prior: Default prior type ('uniform', 'normal', 'log-uniform')
        default_bounds: [min, max] for uniform/log-uniform priors
        default_mean: Mean for normal prior (optional)
        default_std: Standard deviation for normal prior (optional)
        units: Physical units (e.g., 'Pa·s', 'K', 'S/m')
        requires_structure_rebuild: If True, varying this parameter requires
            regenerating structure cache (e.g., Tb_K, ocean salinity)
        rheology_constraint: Required rheology ('andrade', 'maxwell', None=any)
    """
    id: str
    label: str
    latex_label: str
    description: str
    category: str
    default_prior: str
    default_bounds: List[float]
    default_mean: Optional[float] = None
    default_std: Optional[float] = None
    units: Optional[str] = None
    requires_structure_rebuild: bool = False
    rheology_constraint: Optional[str] = None


# ============================================================================
# Parameter Registry
# ============================================================================

PARAMETER_REGISTRY: Dict[str, ParameterDef] = {
    # -------------------------------------------------------------------------
    # Rheology Parameters (Andrade)
    # -------------------------------------------------------------------------
    'log10_wOcean_ppt': ParameterDef(
        id='log10_wOcean_ppt',
        label='Log₁₀(Ocean Salinity)',
        latex_label=r'$\log_{10}(w_{ocean}/\mathrm{ppt})$',
        description='Ocean salinity (log10 of parts per thousand). '
                    'Drives ocean conductivity (all induction channels), '
                    'density, freezing point, and layer thicknesses. '
                    'Requires a salinity-resolved (2D Tb x w) structure '
                    'cache; range 0.1-100 ppt (user 2026-07-18).',
        category='ocean',
        default_prior='uniform',
        default_bounds=[-1.0, 2.0],
        units='log10(ppt)',
        requires_structure_rebuild=True,
        rheology_constraint=None
    ),

    # -------------------------------------------------------------------------
    # Gravity nuisance parameters (v4 geodesy: sampled non-hydrostaticity)
    # -------------------------------------------------------------------------
    'dC20_nh': ParameterDef(
        id='dC20_nh',
        label='Non-hydrostatic ΔC₂₀',
        latex_label=r'$\Delta C_{20}^{nh}$',
        description='Additive non-hydrostatic offset on the unnormalized '
                    'C20 gravity coefficient. Sampled alongside the '
                    'physical parameters so real non-hydrostatic power '
                    'cannot bias the interior posterior; its posterior IS '
                    'the inferred degree of non-hydrostaticity. Requires '
                    "gravity_forward_model='clairaut_hydrostatic' "
                    '(v4 geodesy plan).',
        category='gravity',
        default_prior='uniform',
        default_bounds=[-2.0e-5, 2.0e-5],
        units='dimensionless',
        requires_structure_rebuild=False,
        rheology_constraint=None
    ),

    'dC22_nh': ParameterDef(
        id='dC22_nh',
        label='Non-hydrostatic ΔC₂₂',
        latex_label=r'$\Delta C_{22}^{nh}$',
        description='Additive non-hydrostatic offset on the unnormalized '
                    'C22 gravity coefficient (see dC20_nh). Only the '
                    'ratio-breaking combination is identifiable; the '
                    'ratio-preserving one is degenerate with a C/MR² '
                    'shift (v4 plan, review R4).',
        category='gravity',
        default_prior='uniform',
        default_bounds=[-2.0e-5, 2.0e-5],
        units='dimensionless',
        requires_structure_rebuild=False,
        rheology_constraint=None
    ),

    'alpha': ParameterDef(
        id='alpha',
        label='Andrade Exponent',
        latex_label=r'$\alpha$',
        description='Andrade rheology frequency exponent (dimensionless). '
                    'Controls viscoelastic dissipation spectrum width.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[0.2, 0.4],
        units=None,
        requires_structure_rebuild=False,
        rheology_constraint='andrade'
    ),

    'log10_zeta': ParameterDef(
        id='log10_zeta',
        label='Log₁₀(ζ Andrade Compliance)',
        latex_label=r'$\log_{10}(\zeta)$',
        description='Universal Andrade compliance applied to all solid layers. '
                    'Controls the strength of the viscoelastic dissipation.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[-4, 2],
        units=None,
        requires_structure_rebuild=False,
        rheology_constraint='andrade'
    ),

    'log10_zeta_Ih': ParameterDef(
        id='log10_zeta_Ih',
        label='Log₁₀(ζ Ice Ih)',
        latex_label=r'$\log_{10}(\zeta_{\rm Ih})$',
        description='Andrade compliance for Ice Ih. PP default: 0.005 (log₁₀ ≈ -2.3). '
                    'Controls viscoelastic dissipation in the ice shell.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[-4, 0],
        units=None,
        requires_structure_rebuild=False,
        rheology_constraint='andrade'
    ),

    'log10_zeta_HP': ParameterDef(
        id='log10_zeta_HP',
        label='Log₁₀(ζ HP Ice)',
        latex_label=r'$\log_{10}(\zeta_{\rm HP})$',
        description='Andrade compliance for high-pressure ice (III, V, VI). '
                    'PP default: 0.05 (log₁₀ ≈ -1.3). HP ice is typically '
                    '~10× more compliant than Ice Ih.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[-4, 2],
        units=None,
        requires_structure_rebuild=False,
        rheology_constraint='andrade'
    ),

    'log10_zeta_sil': ParameterDef(
        id='log10_zeta_sil',
        label='Log₁₀(ζ Silicate)',
        latex_label=r'$\log_{10}(\zeta_{\rm sil})$',
        description='Andrade compliance for silicate mantle. PP default: 0.005 '
                    '(log₁₀ ≈ -2.3). Typically similar to Ice Ih.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[-4, 0],
        units=None,
        requires_structure_rebuild=False,
        rheology_constraint='andrade'
    ),

    # -------------------------------------------------------------------------
    # Rheology Parameters (Viscosities - both Andrade and Maxwell)
    # -------------------------------------------------------------------------
    'log10_eta_Ih': ParameterDef(
        id='log10_eta_Ih',
        label='Log₁₀(Ice Ih Viscosity)',
        latex_label=r'$\log_{10}(\eta_{\rm Ih})$',
        description='Ice Ih shell viscosity. Dominant control on k₂ response '
                    'at Saturnian tidal periods.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[12, 16],
        units='Pa·s',
        requires_structure_rebuild=False,
        rheology_constraint=None  # Available for both Andrade and Maxwell
    ),

    'log10_eta_HP': ParameterDef(
        id='log10_eta_HP',
        label='Log₁₀(HP Ice Viscosity)',
        latex_label=r'$\log_{10}(\eta_{\rm HP})$',
        description='High-pressure ice (III, V, VI) viscosity. Affects deep '
                    'layer contribution to tidal Love number.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[10, 18],
        units='Pa·s',
        requires_structure_rebuild=False,
        rheology_constraint=None
    ),

    'log10_eta_III': ParameterDef(
        id='log10_eta_III',
        label='Log₁₀(Ice III Viscosity)',
        latex_label=r'$\log_{10}(\eta_{\rm III})$',
        description='Ice III viscosity.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[10, 18],
        units='Pa·s',
        requires_structure_rebuild=False,
        rheology_constraint=None
    ),

    'log10_eta_V': ParameterDef(
        id='log10_eta_V',
        label='Log₁₀(Ice V Viscosity)',
        latex_label=r'$\log_{10}(\eta_{\rm V})$',
        description='Ice V viscosity.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[10, 18],
        units='Pa·s',
        requires_structure_rebuild=False,
        rheology_constraint=None
    ),

    'log10_eta_VI': ParameterDef(
        id='log10_eta_VI',
        label='Log₁₀(Ice VI Viscosity)',
        latex_label=r'$\log_{10}(\eta_{\rm VI})$',
        description='Ice VI viscosity.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[10, 18],
        units='Pa·s',
        requires_structure_rebuild=False,
        rheology_constraint=None
    ),

    'log10_eta_sil': ParameterDef(
        id='log10_eta_sil',
        label='Log₁₀(Silicate Viscosity)',
        latex_label=r'$\log_{10}(\eta_{\rm sil})$',
        description='Silicate mantle viscosity. Negligible for tidal response '
                    'at satellite periods, but affects long-term evolution.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[12, 22],
        units='Pa·s',
        requires_structure_rebuild=False,
        rheology_constraint=None
    ),

    # -------------------------------------------------------------------------
    # Rheology Parameters (Maxwell only)
    # -------------------------------------------------------------------------
    'log10_mu_Ih': ParameterDef(
        id='log10_mu_Ih',
        label='Log₁₀(Ice Ih Shear Modulus)',
        latex_label=r'$\log_{10}(\mu_{\rm Ih})$',
        description='Ice Ih shear modulus for Maxwell rheology. Typically '
                    '3-4 GPa at satellite temperatures.',
        category='rheology',
        default_prior='uniform',
        default_bounds=[8, 10],
        units='Pa',
        requires_structure_rebuild=False,
        rheology_constraint='maxwell'
    ),

    # -------------------------------------------------------------------------
    # Structural Parameters (require structure rebuild)
    # -------------------------------------------------------------------------
    'Tb_K': ParameterDef(
        id='Tb_K',
        label='Bottom Temperature',
        latex_label=r'$T_b$ (K)',
        description='Ocean-ice interface temperature. Controls HP ice layer '
                    'thicknesses and phase transitions. Requires pre-computed '
                    'structure grid for inference.',
        category='structure',
        default_prior='uniform',
        default_bounds=[251.2, 269.0],
        units='K',
        requires_structure_rebuild=True,
        rheology_constraint=None
    ),

    'D_iceIh_km': ParameterDef(
        id='D_iceIh_km',
        label='Ice Shell Thickness',
        latex_label=r'$D_{\rm ice\,Ih}$ (km)',
        description='Ice-Ih shell thickness (v5 reparameterization). Sampled '
                    'in place of Tb_K so the informative prior (29 +/- 10 km, '
                    'truncated Gaussian) is placed on the salinity-INDEPENDENT '
                    'ice thickness rather than on Tb (a uniform Tb prior is an '
                    'implicit salinity-dependent thickness prior). Tb is derived '
                    'per draw by inverting the cached D_iceIh(Tb, w) field at the '
                    'drawn salinity. Requires a 2D (Tb x w) structure grid.',
        category='structure',
        default_prior='truncated_gaussian',
        default_bounds=[5.0, 60.0],
        units='km',
        requires_structure_rebuild=True,
        rheology_constraint=None
    ),

    'rhoSilInput_kgm3': ParameterDef(
        id='rhoSilInput_kgm3',
        label='Silicate Density',
        latex_label=r'$\rho_{\rm sil}$ (kg m$^{-3}$)',
        description=(
            'Self-consistent-grid control parameter for Titan silicate density. '
            'Maps to Planet.Sil.rhoSilWithCore_kgm3 and forces '
            'Planet.Do.CONSTANT_INNER_DENSITY=True during structure-cache generation. '
            'This gives direct density control but relaxes full EOS self-consistency.'
        ),
        category='structure',
        default_prior='uniform',
        default_bounds=[2400.0, 3500.0],
        units='kg/m^3',
        requires_structure_rebuild=True,
        rheology_constraint=None
    ),

    'D_hydro_km': ParameterDef(
        id='D_hydro_km',
        label='Hydrosphere Thickness',
        latex_label=r'$D_{\rm hydro}$ (km)',
        description=(
            'Hybrid PPTest45 structural control: total hydrosphere thickness '
            'from surface to silicate top. Hydrosphere PT remains controlled by '
            'Tb_K; Mtot_kg and CMR2mean are computed outputs.'
        ),
        category='structure',
        default_prior='uniform',
        default_bounds=[300.0, 700.0],
        units='km',
        requires_structure_rebuild=True,
        rheology_constraint=None
    ),

    'R_core_km': ParameterDef(
        id='R_core_km',
        label='Core Radius',
        latex_label=r'$R_{\rm core}$ (km)',
        description=(
            'Radius of a differentiated metallic/rock core. Sampled jointly with '
            'rho_core_kgm3; silicate mantle density (rho_sil_kgm3) is then derived '
            'from mass conservation against the cached (core-independent) '
            'hydrosphere and rejected if outside its configured bounds. Does not '
            'require regenerating the Tb-grid structure cache -- the core is '
            'assembled at MCMC runtime from the cached hydrosphere layers.'
        ),
        category='structure',
        default_prior='uniform',
        default_bounds=[0.0, 2000.0],
        units='km',
        requires_structure_rebuild=False,
        rheology_constraint=None
    ),

    'rho_core_kgm3': ParameterDef(
        id='rho_core_kgm3',
        label='Core Density',
        latex_label=r'$\rho_{\rm core}$ (kg m$^{-3}$)',
        description=(
            'Bulk density of a differentiated core, sampled jointly with '
            'R_core_km. Used with the derived silicate density (mass '
            'conservation, rho_sil_kgm3) to reconstruct CMR2 for the core + '
            'silicate + cached-hydrosphere shell stack. Bounds are '
            'body/scenario-specific (e.g. rock-differentiation core: '
            '[2591, 3600] kg/m^3 per Petricca et al. 2025 for Titan; metallic '
            'Fe/Fe-S core: much higher, e.g. Callisto [5000, 8000] kg/m^3).'
        ),
        category='structure',
        default_prior='uniform',
        default_bounds=[2591.0, 3600.0],
        units='kg/m^3',
        requires_structure_rebuild=False,
        rheology_constraint=None
    ),

    # Note: Future structural parameters to add when grid caching is implemented:
    # - 'ocean_salinity_ppt': Ocean salinity (affects ocean density, freezing point)
    # - 'ice_shell_thickness_km': Direct ice shell thickness constraint
}


# ============================================================================
# Category Display Configuration
# ============================================================================

CATEGORY_LABELS = {
    'rheology': '🔧 Rheology Parameters',
    'structure': '🌍 Structural Parameters',
    'ocean': '🌊 Ocean Properties',
    'magnetic': '🧲 Magnetic Properties',
    'gravity': '⚖️ Gravity (non-hydrostatic offsets)',
}

CATEGORY_ORDER = ['rheology', 'structure', 'ocean', 'magnetic', 'gravity']


# ============================================================================
# Presets
# ============================================================================

PARAMETER_PRESETS = {
    'andrade_titan': {
        'name': 'Andrade Rheology (Titan)',
        'description': 'Extended 7-param Andrade model for Titan no-ocean (PPTest41). '
                       'Uses per-phase zeta (Ih/HP/sil independent), unlike Test41 which uses '
                       'a single log10_zeta across all solid layers (5-param). '
                       'Petricca et al. 2025 observational constraints.',
        'parameters': ['alpha', 'log10_zeta_Ih', 'log10_zeta_HP', 'log10_zeta_sil',
                        'log10_eta_Ih', 'log10_eta_HP', 'log10_eta_sil'],
        'observables': {
            'Re_k2': (0.608, 0.048),
            'abs_Im_k2': (0.135, 0.035),
            'CMR2': (0.343, 0.001),  # Petricca et al. (2025)
        },
        'test_module': 'PlanetProfile.Test.PPTest41',
        'rheology': 'andrade'
    },

    'andrade_titan_noocean_8D': {
        'name': 'Andrade Rheology (Titan No-Ocean 8D)',
        'description': 'Test50 8D no-ocean Titan, Andrade rheology, Yao 2014 spherical convection, Tb sampled across triple-pt depression band.',
        'parameters': ['alpha', 'log10_zeta', 'log10_eta_Ih', 'log10_eta_III',
                       'log10_eta_V', 'log10_eta_VI', 'log10_eta_sil', 'Tb_K'],
        # NOTE: keys here byte-match test50_titan_noocean_andrade_8D.json's
        # "observables" block -- 'Im_k2' (NOT 'abs_Im_k2') and NO CMR2, matching
        # the committed Test50 production MCMC result used as the SBI crosscheck
        # reference. Test50 samples no core parameters, so model CMR2 varies only
        # through the 9-pt Tb cache lookup (peak-to-peak 2.5e-5 = 0.025 sigma_obs)
        # -- a CMR2 term adds no posterior information and would make GUI-launched
        # runs diverge from the CLI reference. mcmc_runner.py accepts 'Im_k2' /
        # 'abs_Im_k2' as aliases.
        'observables': {
            'Re_k2': (0.608, 0.048),
            'Im_k2': (0.135, 0.035),
        },
        'test_module': 'PlanetProfile.Test.Test50_mcmc_andrade_noocean_yao2014',
        'rheology': 'andrade'
    },

    'andrade_titan_noocean_diff_10D': {
        'name': 'Andrade Rheology (Titan No-Ocean, Differentiated Core, 10D)',
        'description': (
            'Test52 10D no-ocean Titan: Test50 8D rheology/thermal parameter '
            'space (alpha, log10_zeta, log10_eta_{Ih,III,V,VI,sil}, Tb_K) plus a '
            'sampled differentiated core (R_core_km, rho_core_kgm3) with rho_sil '
            'derived from mass conservation against the cached hydrosphere. '
            'Petricca et al. (2025) Re_k2/Im_k2/CMR2 constraints.'
        ),
        'parameters': ['alpha', 'log10_zeta', 'log10_eta_Ih', 'log10_eta_III',
                       'log10_eta_V', 'log10_eta_VI', 'log10_eta_sil', 'Tb_K',
                       'R_core_km', 'rho_core_kgm3'],
        # NOTE: keys here byte-match test52_titan_noocean_andrade_10D.json's
        # "observables" block -- 'Im_k2' (NOT 'abs_Im_k2', unlike the other
        # presets below). mcmc_runner.py accepts both spellings as aliases,
        # but this preset must match the config file verbatim per the Test52
        # Phase 1 plan (plans/test52-differentiated-titan-plan.md).
        'observables': {
            'Re_k2': (0.608, 0.048),
            'Im_k2': (0.135, 0.035),
            'CMR2': (0.343, 0.001),  # Petricca et al. (2025)
        },
        'test_module': 'PlanetProfile.Test.PPTest50',
        'rheology': 'andrade'
    },

    'maxwell_titan': {
        'name': 'Maxwell Rheology + Ocean (Titan)',
        'description': 'Maxwell viscoelastic ocean model for Titan (PPTest42) — 4-param space matching Test42_mcmc_maxwell_ocean.py',
        'parameters': ['log10_eta_Ih', 'log10_eta_HP', 'log10_eta_sil', 'Tb_K'],
        'observables': {
            'Re_k2': (0.608, 0.048),
            'abs_Im_k2': (0.135, 0.035),
            'CMR2': (0.343, 0.001),  # Petricca et al. (2025)
        },
        'test_module': 'PlanetProfile.Test.PPTest42',
        'rheology': 'maxwell'
    },

    'maxwell_titan_hybrid_hydrosphere': {
        'name': 'Maxwell Rheology + Hybrid Hydrosphere Thickness (Titan)',
        'description': (
            'PPTest45 prototype. Keeps PlanetProfile hydrosphere PT self-consistent '
            'in Tb_K, but samples D_hydro_km instead of using MoI closure to set '
            'the hydrosphere/silicate boundary.'
        ),
        'parameters': ['log10_eta_Ih', 'log10_eta_HP', 'log10_eta_sil', 'Tb_K', 'D_hydro_km'],
        'observables': {
            'Re_k2': (0.608, 0.048),
            'abs_Im_k2': (0.135, 0.035),
            'CMR2': (0.343, 0.001),
            'Mtot_kg': (1.3452e23, 1.0e20),
        },
        'test_module': 'PlanetProfile.Test.PPTest45',
        'rheology': 'maxwell'
    },

    'andrade_europa': {
        # LEGACY (2026-07-12): 5D PPTest46 reference model (Titan-derived
        # PureH2O structure). NOT the Europa campaign model — that is the v2
        # 7D seawater config configs/europa_seawater_andrade_7D.json (sampled
        # Tb/core, mass-conservation rho_sil, Mazarico-assumption k2/h2,
        # |Ae_synodic| > 0.7 induction bound). Retired from the GUI preset
        # radio; kept for provenance of the Variant A/B study.
        # Variant B reference structure: intermediate Fe-FeS core (xFeS=0.55, rhoCore≈5865 kg/m³).
        # Thin Ice Ih shell only; no HP ice phases for Europa's ~25 km shell.
        # k2 values are theoretical placeholders — Europa Clipper will constrain these.
        # CMR2 from Gomez Casajus et al. (2021) via Petricca et al. (2025).
        'name': 'Andrade Rheology (Europa, thin shell)',
        'description': (
            'Andrade rheology for Europa thin Ice Ih shell (PPTest46, Variant B). '
            'Fe-FeS core with xFeS=0.55. MOI from Gomez Casajus 2021. '
            'k2 values are theoretical placeholders pending Europa Clipper data.'
        ),
        'parameters': ['alpha', 'log10_zeta_Ih', 'log10_zeta_sil',
                       'log10_eta_Ih', 'log10_eta_sil'],
        'observables': {
            'CMR2': (0.3547, 0.0024),        # Gomez Casajus et al. (2021)
            'Re_k2': (0.328, 0.15),          # Theoretical estimate; large σ reflects ignorance
            'abs_Im_k2': (0.015, 0.010),     # Theoretical estimate; placeholder
        },
        'test_module': 'PlanetProfile.Test.PPTest46',
        'rheology': 'andrade'
    },

    'andrade_europa_nocore': {
        # Variant A: partially undifferentiated Europa (no metallic core).
        # Represents upper limit on partial differentiation per Petricca et al. (2025).
        'name': 'Andrade Rheology (Europa, no core)',
        'description': (
            'Andrade rheology for undifferentiated Europa (PPTest46 Variant A, Fe_CORE=False). '
            'Upper limit on partial differentiation per Petricca et al. (2025).'
        ),
        'parameters': ['alpha', 'log10_zeta_Ih', 'log10_zeta_sil',
                       'log10_eta_Ih', 'log10_eta_sil'],
        'observables': {
            'CMR2': (0.3547, 0.0024),
            'Re_k2': (0.328, 0.15),
            'abs_Im_k2': (0.015, 0.010),
        },
        'test_module': 'PlanetProfile.Test.PPTest46',
        'rheology': 'andrade'
    },
}


# ============================================================================
# Utility Functions
# ============================================================================

def get_parameters_by_category(category: str) -> List[ParameterDef]:
    """Get all parameters in a given category."""
    return [p for p in PARAMETER_REGISTRY.values() if p.category == category]


def get_parameters_by_rheology(rheology: Optional[str] = None) -> List[ParameterDef]:
    """
    Get parameters compatible with a rheology type.

    Args:
        rheology: 'andrade', 'maxwell', or None (all parameters)

    Returns:
        List of ParameterDef objects
    """
    if rheology is None:
        return list(PARAMETER_REGISTRY.values())

    return [
        p for p in PARAMETER_REGISTRY.values()
        if p.rheology_constraint is None or p.rheology_constraint == rheology
    ]


def validate_parameter_combination(param_ids: List[str]) -> Dict[str, Any]:
    """
    Validate that a set of parameters is compatible.

    Args:
        param_ids: List of parameter IDs to check

    Returns:
        Dict with:
            'valid': bool
            'rheology': 'andrade', 'maxwell', or None (inferred)
            'warnings': List[str] (validation messages)
            'requires_rebuild': bool (any param needs structure rebuild)
    """
    warnings = []
    requires_rebuild = False

    # Check for parameters requiring structure rebuild
    for param_id in param_ids:
        param = PARAMETER_REGISTRY.get(param_id)
        if param and param.requires_structure_rebuild:
            requires_rebuild = True
            warnings.append(
                f"Parameter '{param.label}' requires pre-computed structure grid. "
                f"Generate grid with prepare_structure_variants.py before inference."
            )

    # Infer rheology constraint
    andrade_only = [p for p in param_ids if PARAMETER_REGISTRY[p].rheology_constraint == 'andrade']
    maxwell_only = [p for p in param_ids if PARAMETER_REGISTRY[p].rheology_constraint == 'maxwell']

    if andrade_only and maxwell_only:
        warnings.append(
            f"Parameter set mixes Andrade-only ({andrade_only}) and Maxwell-only ({maxwell_only}) parameters. "
            f"Choose parameters compatible with one rheology."
        )
        return {'valid': False, 'rheology': None, 'warnings': warnings, 'requires_rebuild': requires_rebuild}

    # Determine rheology
    if andrade_only:
        rheology = 'andrade'
    elif maxwell_only:
        rheology = 'maxwell'
    else:
        rheology = None  # All parameters are rheology-agnostic
        warnings.append(
            "No rheology-specific parameters selected. Add 'alpha' and 'log10_zeta_Ih'/'log10_zeta_HP'/'log10_zeta_sil' for Andrade, "
            "or 'log10_mu_Ih' for Maxwell."
        )

    return {
        'valid': len(warnings) == 0 or not any('mixes' in w for w in warnings),
        'rheology': rheology,
        'warnings': warnings,
        'requires_rebuild': requires_rebuild
    }
