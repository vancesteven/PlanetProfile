""" Regression test for the MgSO4 ocean electrical-conductivity models in
    PlanetProfile.Thermodynamics.MgSO4.MgSO4Props.MgSO4Conduct.

    Two models are covered:
      - 'Vance2018' (Larionov & Kryukov 1984 extrapolation, Vance et al. 2018) —
        the default; must remain unchanged.
      - 'Pan2020'  (Pan, Yong & Secco 2020, 10 wt% MgSO4) — ported from the
        MATLAB reference Thermodynamics/MgSO4/getSigmaMgSO4_Pan.m.

    The Pan2020 values below are locked to the MATLAB regression (reproduced to
    4 significant figures at the time of the port, scientific-reviewer verified
    2026-07-13). A change to these numbers means the conductivity model changed
    and must be re-reviewed against Pan et al. (2020) Fig. 5.

    Run standalone:  python -m PlanetProfile.Test.TestMgSO4Conductivity
    Or import TestMgSO4Conductivity() from the build-test harness.
"""
import logging
import numpy as np

from PlanetProfile.Thermodynamics.MgSO4.MgSO4Props import MgSO4Conduct

log = logging.getLogger('PlanetProfile')


def TestMgSO4Conductivity():
    """ Assert the MgSO4 conductivity models reproduce their reference values.
        Raises AssertionError on any regression; returns a dict of the sampled
        conductivities on success.
    """
    results = {}

    # --- Pan2020: locked regression values (S/m) from getSigmaMgSO4_Pan.m ---
    # Ported constants: 10 wt% (100 ppt), molality 0.9231. Sample points chosen
    # to bracket the shallow-warm and deep-cold ends of the ocean column.
    pan = MgSO4Conduct(100.0, 'Pan2020')
    pan_shallow = float(pan(10.0, 270.0))    # ~near-surface, warm
    pan_deep = float(pan(800.0, 255.0))      # deep, cold (table edge)
    results['Pan2020'] = {'sigma_10MPa_270K': pan_shallow,
                          'sigma_800MPa_255K': pan_deep}
    log.info(f'Pan2020 sigma(10 MPa, 270 K)  = {pan_shallow:.4f} S/m')
    log.info(f'Pan2020 sigma(800 MPa, 255 K) = {pan_deep:.4f} S/m')

    # Locked to the MATLAB port (4 sig figs). Tight tolerance: this is an exact
    # numerical reproduction, not a physical-range check.
    assert abs(pan_shallow - 2.4460) < 1e-3, \
        f'Pan2020 sigma(10 MPa,270 K) regressed: {pan_shallow:.4f} != 2.4460 S/m'
    assert abs(pan_deep - 0.4998) < 1e-3, \
        f'Pan2020 sigma(800 MPa,255 K) regressed: {pan_deep:.4f} != 0.4998 S/m'

    # Physical band over the table-backed subregion (P<=800 MPa, T>=253.15 K):
    # Pan et al. (2020) Fig. 5 contour levels span 0.17-3.3 S/m.
    P_grid = np.array([10.0, 200.0, 500.0, 800.0])
    T_grid = np.array([255.0, 265.0, 275.0, 290.0])
    sig_grid = pan(P_grid, T_grid, grid=True)
    assert np.all(sig_grid > 0.1) and np.all(sig_grid < 4.0), \
        (f'Pan2020 table-backed sigma outside physical band [0.1, 4.0] S/m: '
         f'min {sig_grid.min():.4f}, max {sig_grid.max():.4f}')

    # --- Vance2018 (default): must still construct and give a sane value ---
    vance = MgSO4Conduct(100.0, 'Vance2018', rhoType='Millero',
                         scalingType='Vance2018')
    vance_shallow = float(vance(10.0, 270.0))
    results['Vance2018'] = {'sigma_10MPa_270K': vance_shallow}
    log.info(f'Vance2018 sigma(10 MPa, 270 K) = {vance_shallow:.4f} S/m')
    assert 0.1 < vance_shallow < 10.0, \
        f'Vance2018 sigma(10 MPa,270 K) unphysical: {vance_shallow:.4f} S/m'

    # --- Unknown model must raise (guard against silent typos) ---
    try:
        MgSO4Conduct(100.0, 'NotAModel')
    except ValueError:
        pass
    else:
        raise AssertionError('MgSO4Conduct accepted an unknown elecType without error')

    log.info('TestMgSO4Conductivity: all assertions passed.')
    return results


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    out = TestMgSO4Conductivity()
    print('PASS — MgSO4 conductivity regression:')
    for model, vals in out.items():
        print(f'  {model}: ' + ', '.join(f'{k}={v:.4f}' for k, v in vals.items()))
