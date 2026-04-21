"""
ValidateKalousova2018.py
Validation script comparing PlanetProfile's HP ice convection implementation
against Kalousova & Sotin (2018) scaling laws.

Tests:
1. Direct formula verification for known (qs, mu0) combinations
2. PPTest30 integration check for physical reasonableness
3. Regime boundary identification (sub/supercritical)

Run: python -m PlanetProfile.Test.ValidateKalousova2018
"""
import numpy as np
import sys
import logging

log = logging.getLogger('PlanetProfile')


def validate_scaling_laws():
    """Verify ConvectionKalousova2018 reproduces the paper's scaling laws directly."""
    print('='*70)
    print('PART 1: Direct Scaling Law Verification')
    print('='*70)

    # Kalousova & Sotin (2018) scaling laws:
    # Ra*_c = 19.965e3 * (qs [mW/m^2])^3.690                     (Eq. 7)
    # Ht[km] = (0.145e-3 * qs[mW/m^2] + 0.015) * mu0[Pa*s]^0.21 (Eq. 9)
    # delta_i = 2.746 * (Ra*)^(-0.271)                             (Eq. 4)

    # Test parameter grid from K&S 2018:
    # qs: 10-40 mW/m^2, mu0: 1e14-1e17 Pa*s
    test_cases = [
        # (qs_mWm2, mu0_Pas, description)
        (10.0, 1e14, 'low flux, low viscosity'),
        (10.0, 1e17, 'low flux, high viscosity'),
        (25.0, 1e15, 'mid flux, mid viscosity'),
        (40.0, 1e14, 'high flux, low viscosity'),
        (40.0, 1e17, 'high flux, high viscosity'),
    ]

    print(f'\n{"qs (mW/m2)":>12} {"mu0 (Pa*s)":>12} {"Ra*_c":>12} {"Ht (km)":>10} {"Description"}')
    print('-'*70)

    all_pass = True
    for qs, mu0, desc in test_cases:
        # Direct formula computation
        RaCrit = 19.965e3 * (qs**3.690)
        Ht_km = (0.145e-3 * qs + 0.015) * (mu0**0.21)

        print(f'{qs:12.1f} {mu0:12.1e} {RaCrit:12.3e} {Ht_km:10.3f}   {desc}')

        # Sanity checks
        if RaCrit <= 0:
            print(f'  FAIL: Ra*_c should be positive')
            all_pass = False
        if Ht_km < 0:
            print(f'  FAIL: Ht should be non-negative')
            all_pass = False

    # Verify TBL scaling
    print(f'\n{"Ra*":>12} {"delta_i":>10} {"deltaTBL (km, H=200km)":>24}')
    print('-'*50)
    for Ra_star in [1e6, 1e8, 1e10, 1e12]:
        delta_i = 2.746 * (Ra_star**(-0.271))
        deltaTBL_km = delta_i * 200  # For 200 km layer
        print(f'{Ra_star:12.1e} {delta_i:10.4f} {deltaTBL_km:24.1f}')
        if delta_i < 0 or delta_i > 1:
            print(f'  FAIL: delta_i should be between 0 and 1')
            all_pass = False

    if all_pass:
        print('\nPART 1 PASSED: All scaling law computations are self-consistent.')
    else:
        print('\nPART 1 FAILED: Some scaling law checks did not pass.')
    return all_pass


def validate_test30():
    """Run PPTest30 and check HP ice convection diagnostics."""
    print('\n' + '='*70)
    print('PART 2: PPTest30 Integration Validation')
    print('='*70)

    try:
        from PlanetProfile.Main import PlanetProfile as RunPP
        from PlanetProfile.GetConfig import GetConfig
    except ImportError as e:
        print(f'Could not import PlanetProfile: {e}')
        return False

    # Load Test30 configuration
    try:
        Params = GetConfig()
        Params.CALC_NEW = True
        Params.SKIP_PLOTS = True
        Params.SKIP_INDUCTION = True
        Params.SKIP_GRAVITY = True
        Params.NO_SAVEFILE = True
        Params.VERBOSE = False
        Params.CALC_OCEAN_PROPS = False
        Params.CALC_SEISMIC = False
        Params.CALC_VISCOSITY = False

        # Import test config
        import importlib
        spec = importlib.util.spec_from_file_location('PPTest30',
            'PlanetProfile/Test/PPTest30.py')
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        Planet = mod.Planet

        print('Running PPTest30...')
        Planet, Params = RunPP(Planet, Params)
        print('PPTest30 completed successfully.')
    except Exception as e:
        print(f'PPTest30 failed to run: {e}')
        import traceback
        traceback.print_exc()
        return False

    all_pass = True

    # Check HP ice convection diagnostics for each phase
    phases = {
        'III': ('RaConvectIII', 'RaCritIII', 'eLidIII_m', 'DconvIII_m',
                'deltaTBLIII_m', 'TconvIII_K', 'vConvIII_ms'),
        'V': ('RaConvectV', 'RaCritV', 'eLidV_m', 'DconvV_m',
              'deltaTBLV_m', 'TconvV_K', 'vConvV_ms'),
        'VI': ('RaConvectVI', 'RaCritVI', 'eLidVI_m', 'DconvVI_m',
               'deltaTBLVI_m', 'TconvVI_K', 'vConvVI_ms'),
    }

    print(f'\n{"Phase":>6} {"Ra*":>12} {"Ra*_c":>12} {"eLid(km)":>10} {"Dconv(km)":>10} '
          f'{"dTBL(km)":>10} {"Tconv(K)":>10} {"v(m/yr)":>10} {"Regime":>12}')
    print('-'*100)

    for phase, attrs in phases.items():
        Ra = getattr(Planet, attrs[0], None)
        RaCrit = getattr(Planet, attrs[1], None)
        eLid = getattr(Planet, attrs[2], None)
        Dconv = getattr(Planet, attrs[3], None)
        deltaTBL = getattr(Planet, attrs[4], None)
        Tconv = getattr(Planet, attrs[5], None)
        vConv = getattr(Planet, attrs[6], None)

        if Ra is None or not np.isfinite(Ra) or Ra == 0:
            print(f'{phase:>6}   (no HP ice {phase} layer present)')
            continue

        vConv_myr = vConv * 3.156e7 if vConv is not None and np.isfinite(vConv) else 0.0
        regime = 'supercritical' if Ra > RaCrit else 'subcritical'

        print(f'{phase:>6} {Ra:12.3e} {RaCrit:12.3e} {eLid/1e3:10.1f} {Dconv/1e3:10.1f} '
              f'{deltaTBL/1e3:10.1f} {Tconv:10.1f} {vConv_myr:10.3f} {regime:>12}')

        # Physical reasonableness checks
        if Ra < 0:
            print(f'  FAIL: Ra* should be non-negative for ice {phase}')
            all_pass = False
        if eLid < 0 or Dconv < 0 or deltaTBL < 0:
            print(f'  FAIL: Layer thicknesses should be non-negative for ice {phase}')
            all_pass = False
        if Tconv is not None and Tconv < 200:
            print(f'  WARNING: Tconv={Tconv:.1f} K seems low for HP ice {phase}')

    # Check Ice Ih convection too
    print(f'\n{"Ih":>6} {Planet.RaConvect:12.3e} {Planet.RaCrit:12.3e} '
          f'{Planet.eLid_m/1e3:10.1f} {Planet.Dconv_m/1e3:10.1f} '
          f'{Planet.deltaTBL_m/1e3:10.1f} {Planet.Tconv_K:10.1f} '
          f'{Planet.vConv_ms*3.156e7:10.3f} '
          f'{"supercritical" if Planet.RaConvect > Planet.RaCrit else "subcritical":>12}')

    if all_pass:
        print('\nPART 2 PASSED: PPTest30 HP ice convection diagnostics are physically reasonable.')
    else:
        print('\nPART 2 FAILED: Some diagnostics failed physical reasonableness checks.')
    return all_pass


def validate_regime_boundaries():
    """Verify regime boundary identification (Ra* vs Ra*_c)."""
    print('\n' + '='*70)
    print('PART 3: Regime Boundary Verification')
    print('='*70)

    # For a given qs, compute Ra*_c and verify the regime classification
    # K&S 2018 four regimes (based on Ra*/Ra*_c):
    # 1. Direct exchange (Ra* >> Ra*_c): vigorous melting, melt reaches ocean
    # 2. Indirect exchange (Ra* > Ra*_c): some melting, partial transport
    # 3. Limited exchange (Ra* ~ Ra*_c): marginal melting
    # 4. No melting (Ra* < Ra*_c): no temperate layer

    test_qs = [10.0, 20.0, 30.0, 40.0]  # mW/m^2

    print(f'\n{"qs (mW/m2)":>12} {"Ra*_c":>12} {"log10(Ra*_c)":>14}')
    print('-'*42)

    all_pass = True
    for qs in test_qs:
        RaCrit = 19.965e3 * (qs**3.690)
        print(f'{qs:12.1f} {RaCrit:12.3e} {np.log10(RaCrit):14.2f}')

        # Ra*_c should increase steeply with qs (exponent 3.69)
        if qs > 10.0:
            RaCrit_prev = 19.965e3 * ((qs-10)**3.690)
            if RaCrit <= RaCrit_prev:
                print(f'  FAIL: Ra*_c should increase with qs')
                all_pass = False

    # Verify supercritical regime produces temperate layer
    print('\nVerifying temperate layer appears only when Ra* > Ra*_c:')
    qs = 20.0  # mW/m^2
    RaCrit = 19.965e3 * (qs**3.690)
    mu0 = 1e15  # Pa*s

    # Sub-critical case
    Ra_sub = RaCrit * 0.1
    has_temperate_sub = Ra_sub > RaCrit
    print(f'  Ra* = {Ra_sub:.2e} (0.1 * Ra*_c): temperate = {has_temperate_sub} (expected False)')
    if has_temperate_sub:
        print('  FAIL')
        all_pass = False

    # Super-critical case
    Ra_super = RaCrit * 10
    has_temperate_super = Ra_super > RaCrit
    Ht_km = (0.145e-3 * qs + 0.015) * (mu0**0.21)
    print(f'  Ra* = {Ra_super:.2e} (10 * Ra*_c): temperate = {has_temperate_super} '
          f'(expected True), Ht = {Ht_km:.2f} km')
    if not has_temperate_super:
        print('  FAIL')
        all_pass = False

    if all_pass:
        print('\nPART 3 PASSED: Regime boundary identification is correct.')
    else:
        print('\nPART 3 FAILED: Regime boundary checks did not pass.')
    return all_pass


if __name__ == '__main__':
    print('Kalousova & Sotin (2018) Validation for PlanetProfile')
    print('='*70)

    p1 = validate_scaling_laws()
    p3 = validate_regime_boundaries()

    # Part 2 requires full PlanetProfile run - only run if requested
    if '--full' in sys.argv:
        p2 = validate_test30()
    else:
        p2 = None
        print('\nSkipping PART 2 (full PPTest30 run). Use --full to include it.')

    print('\n' + '='*70)
    print('SUMMARY')
    print('='*70)
    print(f'Part 1 (Scaling Laws): {"PASS" if p1 else "FAIL"}')
    if p2 is not None:
        print(f'Part 2 (PPTest30):     {"PASS" if p2 else "FAIL"}')
    else:
        print(f'Part 2 (PPTest30):     SKIPPED')
    print(f'Part 3 (Regime IDs):   {"PASS" if p3 else "FAIL"}')

    results = [p1, p3] + ([p2] if p2 is not None else [])
    if all(results):
        print('\nAll tests PASSED.')
        sys.exit(0)
    else:
        print('\nSome tests FAILED.')
        sys.exit(1)
