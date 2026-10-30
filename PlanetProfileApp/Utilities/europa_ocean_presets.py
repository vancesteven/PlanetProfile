"""Published Europa ocean composition models (mol/kg), for seeding the
CustomSolution ion table in Ocean Settings.

Source: user-provided literature table (papers/, screenshot 2026-07-17)
summarizing published predictions of Europa's ocean composition. The
'Seawater w/ NH3' entry requires an ammonia-capable database
(frezchemNH3). Concentrations are as published; some models are not
charge-balanced as transcribed (the source table's minor species are
truncated) — Reaktoro equilibrates what it is given.
"""

EUROPA_OCEAN_PRESETS = {
    'Reduced Ocean (Zolotov 2009)': {
        'Na+': 0.301, 'Cl-': 0.301, 'K+': 0.02, 'Ca+2': 0.0001},
    'Sulfate Rich (Zolotov 2009)': {
        'Na+': 0.301, 'Cl-': 0.301, 'K+': 0.01, 'Ca+2': 0.1,
        'SO4-2': 0.301, 'Mg+2': 0.1},
    'Bulk Silicate (Zolotov 2001)': {
        'Na+': 0.049, 'Cl-': 0.021, 'K+': 0.002, 'Ca+2': 0.01,
        'SO4-2': 0.087, 'Mg+2': 0.062},
    'Melwani Daswani 2021 (5 km crust)': {
        'Na+': 0.1, 'Cl-': 0.02, 'K+': 0.008, 'Ca+2': 0.01,
        'SO4-2': 0.07, 'Mg+2': 0.07},
    'Melwani Daswani 2021 (30 km crust)': {
        'Na+': 0.1, 'Cl-': 0.02, 'K+': 0.008, 'Ca+2': 0.02,
        'SO4-2': 0.09, 'Mg+2': 0.08},
    'Seawater w/ NH3': {
        'Na+': 0.486, 'Cl-': 0.566, 'K+': 0.011, 'Ca+2': 0.011,
        'SO4-2': 0.029, 'Mg+2': 0.055, 'NH3': 0.5},
    'Standard Seawater (Millero 2008)': {
        'Na+': 0.486, 'Cl-': 0.566, 'K+': 0.011, 'Ca+2': 0.011,
        'SO4-2': 0.029, 'Mg+2': 0.055},
}
