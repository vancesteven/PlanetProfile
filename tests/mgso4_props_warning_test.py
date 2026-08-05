"""Regression coverage for MgSO4 property-table pressure clamping."""

import importlib
import logging

import numpy as np


mgso4_props = importlib.import_module(
    'PlanetProfile.Thermodynamics.MgSO4.MgSO4Props')


def _pressure_clamp_messages(records):
    return [
        record.getMessage() for record in records
        if record.getMessage().startswith(
            'MgSO4 properties lookup table pressure ceiling')
    ]


def test_pressure_clamp_warns_once_without_changing_properties(caplog, monkeypatch):
    P_cross_MPa = np.array([0.0, 400.0, 900.0, 1800.0])
    P_clamped_MPa = np.array([0.0, 400.0, 800.0, 800.0])
    T_K = np.array([260.0, 290.0])

    monkeypatch.setattr(mgso4_props, '_MGSO4_PRESSURE_CLAMP_WARNED', False)
    with caplog.at_level(logging.WARNING, logger='PlanetProfile'):
        crossed = mgso4_props.MgSO4Props(
            P_cross_MPa, T_K, wOcean_ppt=100.0, EXTRAP=False)
        mgso4_props.MgSO4Props(
            np.array([0.0, 800.0, 1200.0]), T_K,
            wOcean_ppt=100.0, EXTRAP=False)

    messages = _pressure_clamp_messages(caplog.records)
    assert len(messages) == 1
    assert '800 MPa' in messages[0]
    assert '1800 MPa' in messages[0]
    assert 'clamped to the nearest table value' in messages[0]

    explicitly_clamped = mgso4_props.MgSO4Props(
        P_clamped_MPa, T_K, wOcean_ppt=100.0, EXTRAP=False)
    for actual, expected in zip(crossed, explicitly_clamped):
        assert np.array_equal(actual, expected, equal_nan=True)


def test_in_range_pressure_does_not_warn(caplog, monkeypatch):
    monkeypatch.setattr(mgso4_props, '_MGSO4_PRESSURE_CLAMP_WARNED', False)
    with caplog.at_level(logging.WARNING, logger='PlanetProfile'):
        mgso4_props.MgSO4Props(
            np.array([0.0, 400.0, 800.0]), np.array([260.0, 290.0]),
            wOcean_ppt=100.0, EXTRAP=False)

    assert _pressure_clamp_messages(caplog.records) == []
