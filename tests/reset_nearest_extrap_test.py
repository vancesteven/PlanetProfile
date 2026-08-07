"""Regression: ResetNearestExtrap must not truncate float domain bounds
when handed int-typed scalars (the Tb=252 K cache-defect root cause,
plans/MACHINE-B-HANDOFF.md section 0.15)."""
import numpy as np

from PlanetProfile.Utilities.DataManip import ResetNearestExtrap


def test_int_scalar_matches_float_scalar():
    # int 230 beyond max1=229.96 must clamp to 229.96, not truncate to 229
    o1_int, o2_int = ResetNearestExtrap(230, 260, 0.0, 229.96, 240.0, 280.0)
    o1_flt, o2_flt = ResetNearestExtrap(230.0, 260.0, 0.0, 229.96, 240.0,
                                        280.0)
    assert float(o1_int) == float(o1_flt) == 229.96
    assert float(o2_int) == float(o2_flt) == 260.0
    assert np.asarray(o1_int).dtype.kind == 'f'


def test_int_scalar_temperature_axis():
    _, o2 = ResetNearestExtrap(100.0, 240, 0.0, 500.0, 245.5, 280.0)
    assert float(o2) == 245.5
    assert np.asarray(o2).dtype.kind == 'f'


def test_in_domain_values_unchanged():
    o1, o2 = ResetNearestExtrap(100, 260, 0.0, 229.96, 240.0, 280.0)
    assert float(o1) == 100.0 and float(o2) == 260.0


def test_array_inputs_unaffected():
    p = np.array([100, 230, 300])
    t = np.array([250.0, 260.0, 300.0])
    o1, o2 = ResetNearestExtrap(p, t, 0.0, 229.96, 240.0, 280.0)
    assert np.allclose(o1, [100.0, 229.96, 229.96])
    assert np.allclose(o2, [250.0, 260.0, 280.0])
