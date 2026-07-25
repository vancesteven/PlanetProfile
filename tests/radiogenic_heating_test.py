"""Validation of Thermodynamics.RadiogenicHeating against McDonough,
Šrámek & Wipperfurth (2020), doi:10.1029/2019GC008865.

- Per-element present-day specific powers must reproduce the paper's
  Equation-1 coefficients (h[uW/m3] = rho*(3.387e-3 K + 0.01139 Rb +
  0.04595 Sm + 26.18 Th + 98.29 U)) to <1%.
- The BSE inventory times the silicate-Earth mass (4.04e24 kg) must
  reproduce the paper's 19.9 +/- 3.0 TW present-day radiogenic power.
- Paleo-heating must increase monotonically back in time and be
  40K-dominated in the early solar system.
- SetupInit opt-in must override Qrad_Wkg; default None must not.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Thermodynamics.RadiogenicHeating import (
    INVENTORY_BSE_MCD20, INVENTORY_CI_MS95, LLR_SYSTEMS,
    scale_inventory, specific_heating_Wkg,
    specific_power_W_per_kg_element)

# Paper Equation 1 coefficients, uW per kg ELEMENT
EQ1_UW_PER_KG = {'K': 3.387e-3, 'Rb': 0.01139, 'Sm': 0.04595,
                 'Th': 26.18, 'U': 98.29}


def test_eq1_coefficients():
    for el, coeff in EQ1_UW_PER_KG.items():
        p = specific_power_W_per_kg_element(el) * 1e6  # uW/kg
        assert abs(p - coeff) / coeff < 0.01, (el, p, coeff)


def test_bse_total_power():
    M_BSE_KG = 4.04e24
    P_TW = specific_heating_Wkg(INVENTORY_BSE_MCD20) * M_BSE_KG / 1e12
    assert abs(P_TW - 19.9) < 1.0, P_TW


def test_preset_names():
    assert specific_heating_Wkg('BSE') == pytest.approx(
        specific_heating_Wkg(INVENTORY_BSE_MCD20))
    assert specific_heating_Wkg('CI') == pytest.approx(
        specific_heating_Wkg(INVENTORY_CI_MS95))
    with pytest.raises(ValueError):
        specific_heating_Wkg('nope')
    with pytest.raises(ValueError):
        specific_heating_Wkg({'Pu': 1e-9})


def test_paleo_monotonic_and_40K_dominance():
    ts = [0.0, 1.0, 2.0, 3.0, 4.0, 4.5]
    qs = [specific_heating_Wkg(INVENTORY_CI_MS95, t) for t in ts]
    assert all(b > a for a, b in zip(qs, qs[1:])), qs
    # At 4.5 Ga, 40K (t_half 1.25 Ga) dominates a CI inventory
    pK = (INVENTORY_CI_MS95['K']
          * specific_power_W_per_kg_element('K', 4.5))
    assert pK / qs[-1] > 0.5, (pK, qs[-1])
    # Early heating substantially exceeds present (paper Figure 4 sense)
    assert qs[-1] / qs[0] > 4.0, qs


def test_scale_inventory():
    half = scale_inventory(INVENTORY_CI_MS95, 0.5)
    assert specific_heating_Wkg(half) == pytest.approx(
        0.5 * specific_heating_Wkg(INVENTORY_CI_MS95))


def test_setupinit_optin_override():
    from PlanetProfile.Thermodynamics.RadiogenicHeating import (
        qrad_from_inventory)

    class _Sil:
        pass

    class _P:
        pass

    P = _P()
    P.Sil = _Sil()
    P.Sil.radionuclideInventory = None
    assert qrad_from_inventory(P) is None  # default: no override

    P.Sil.radionuclideInventory = 'CI'
    P.Sil.inventoryTime_Ga = 0.0
    q = qrad_from_inventory(P)
    assert q == pytest.approx(specific_heating_Wkg('CI'))
    # CI rock present-day: ~3.4e-12 W/kg ballpark (sanity window)
    assert 2e-12 < q < 6e-12, q
