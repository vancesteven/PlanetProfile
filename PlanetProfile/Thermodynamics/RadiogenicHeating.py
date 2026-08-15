"""Radiogenic heat production from a radionuclide inventory.

Implements the long-lived radionuclide (LLR) heat production of
McDonough, Šrámek & Wipperfurth (2020), "Radiogenic Power and
Geoneutrino Luminosity of the Earth and Other Terrestrial Bodies
Through Time", G-cubed 21, doi:10.1029/2019GC008865
(papers/mcdonough2020radiogenic.pdf; data supplement:
https://data.4tu.nl/articles/_/12681371/2). Per-isotope decay
constants, heat energies per decay (Qh — antineutrino losses already
excluded), and natural isotopic mole fractions are their Table 3.

Opt-in wiring: set ``Planet.Sil.radionuclideInventory`` (element mass
fractions in kg-element/kg-rock, or a named preset) and SetupInit
derives ``Sil.Qrad_Wkg`` from it; the default ``None`` preserves the
existing constant-Qrad_Wkg behavior exactly.

Validation (see tests/radiogenic_heating_test.py): the per-element
present-day specific powers reproduce the paper's Equation-1
coefficients (98.29 uW/kg-U, 26.18 uW/kg-Th, 3.387e-3 uW/kg-K, ...)
to <1%, and the BSE inventory times the silicate-Earth mass (4.04e24
kg) reproduces the paper's 19.9 +/- 3.0 TW.

Scaling to other bodies: icy-moon rock is commonly taken as CI-like;
``INVENTORY_CI_MS95`` (McDonough & Sun 1995 CI chondrite) is provided
as a starting point, and ``scale_inventory`` applies a leaching /
depletion factor. Choosing that factor for a given body is a science
decision — document it in the run config.
"""
from __future__ import annotations

import logging
import numpy as np

log = logging.getLogger('PlanetProfile')

N_AVOGADRO = 6.02214076e23  # 1/mol
S_PER_GA = 3.15576e16       # seconds per Ga (Julian)

# McDonough et al. (2020) Table 3. x_nat: natural mole fraction of the
# isotope in the element; Ar_std: standard atomic weight of the ELEMENT
# (g/mol) — note 235U/238U both use the U standard weight, as in the
# paper; lam_pers: decay constant (1/s); Qh_pJ: heat deposited per
# decay (pJ), antineutrino energy excluded.
LLR_SYSTEMS = {
    '238U':  dict(element='U',  x_nat=0.992740,  Ar_std=238.0289,
                  lam_pers=4.916e-18,  Qh_pJ=7.633),
    '235U':  dict(element='U',  x_nat=0.0072049, Ar_std=238.0289,
                  lam_pers=31.223e-18, Qh_pJ=7.110),
    '232Th': dict(element='Th', x_nat=1.0000,    Ar_std=232.038,
                  lam_pers=1.56e-18,   Qh_pJ=6.475),
    '40K':   dict(element='K',  x_nat=1.167e-4,  Ar_std=39.098,
                  lam_pers=17.40e-18,  Qh_pJ=0.108),
    '87Rb':  dict(element='Rb', x_nat=0.2783,    Ar_std=85.468,
                  lam_pers=0.4428e-18, Qh_pJ=0.013),
    '147Sm': dict(element='Sm', x_nat=0.14993,   Ar_std=150.362,
                  lam_pers=0.2066e-18, Qh_pJ=0.370),
}

# Bulk silicate Earth element mass fractions (kg/kg), McDonough et al.
# (2020) Table 3 — reproduces their 19.9 TW with M_BSE = 4.04e24 kg.
INVENTORY_BSE_MCD20 = {
    'U': 2.00e-8, 'Th': 7.54e-8, 'K': 2.80e-4,
    'Rb': 6.00e-7, 'Sm': 4.06e-7,
}

# CI chondrite element mass fractions (kg/kg), McDonough & Sun (1995),
# The composition of the Earth, Chem. Geol. 120 — the usual undepleted
# reference for icy-moon rock.
INVENTORY_CI_MS95 = {
    'U': 7.4e-9, 'Th': 2.9e-8, 'K': 5.50e-4,
    'Rb': 2.30e-6, 'Sm': 1.48e-7,
}

INVENTORY_PRESETS = {
    'BSE': INVENTORY_BSE_MCD20,
    'CI': INVENTORY_CI_MS95,
}


def scale_inventory(inventory, factor):
    """Uniformly scale an element inventory (e.g. leaching/depletion of
    rock in an undifferentiated or altered body)."""
    return {el: w * float(factor) for el, w in inventory.items()}


def specific_power_W_per_kg_element(element, t_beforePresent_Ga=0.0):
    """Present-day (or paleo) heat production per kg of ELEMENT (W/kg),
    summing that element's LLR systems. t_beforePresent_Ga > 0 rewinds
    each isotope by exp(+lambda t); element mass fractions are defined
    at PRESENT DAY (the standard geochemical convention)."""
    t_s = float(t_beforePresent_Ga) * S_PER_GA
    p = 0.0
    for iso in LLR_SYSTEMS.values():
        if iso['element'] != element:
            continue
        n_per_kg = iso['x_nat'] * N_AVOGADRO / (iso['Ar_std'] * 1e-3)
        p += (n_per_kg * iso['lam_pers'] * iso['Qh_pJ'] * 1e-12
              * np.exp(iso['lam_pers'] * t_s))
    return p


def specific_heating_Wkg(inventory, t_beforePresent_Ga=0.0):
    """Radiogenic heat production per kg of ROCK (W/kg) for an element
    mass-fraction inventory {element: kg-element/kg-rock}."""
    if isinstance(inventory, str):
        try:
            inventory = INVENTORY_PRESETS[inventory]
        except KeyError:
            raise ValueError(
                f"Unknown radionuclide inventory preset '{inventory}'. "
                f"Available: {sorted(INVENTORY_PRESETS)}") from None
    known = {el for iso in LLR_SYSTEMS.values() for el in [iso['element']]}
    extra = set(inventory) - known
    if extra:
        raise ValueError(f'Unknown elements in radionuclide inventory: '
                         f'{sorted(extra)} (known: {sorted(known)})')
    return float(sum(
        w * specific_power_W_per_kg_element(el, t_beforePresent_Ga)
        for el, w in inventory.items()))


def qrad_from_inventory(Planet):
    """Derive Sil.Qrad_Wkg from Planet.Sil.radionuclideInventory (dict
    of element mass fractions, or preset name 'BSE'/'CI'), evaluated at
    Planet.Sil.inventoryTime_Ga before present (default 0 = today).
    Returns the W/kg value; the caller assigns it."""
    inv = getattr(Planet.Sil, 'radionuclideInventory', None)
    if inv is None:
        return None
    t_Ga = float(getattr(Planet.Sil, 'inventoryTime_Ga', 0.0) or 0.0)
    q = specific_heating_Wkg(inv, t_beforePresent_Ga=t_Ga)
    log.info(f'Radiogenic heating from radionuclide inventory '
             f'(McDonough et al. 2020 LLR): Qrad_Wkg = {q:.3e} W/kg'
             + (f' at {t_Ga:g} Ga before present' if t_Ga else '')
             + f' [inventory: {inv if isinstance(inv, str) else "custom"}]')
    return q
