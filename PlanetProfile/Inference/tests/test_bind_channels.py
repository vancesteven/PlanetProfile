"""
Unit tests for the Bind_ induction channel family (Europa Clipper v2).

Bind_<label>_<comp>_<part> is the induced dipole coefficient expressed as
equivalent surface field:

    Bind_comp(f) = Ae(f) * Be_comp(f)   [nT, complex]

conditioned on the SIGNED real/imag parts (never abs-folded, unlike
Im_k2/Im_h2). These tests cover:

  1. _parse_bind_channel grammar (multi-word labels, comp/part suffix).
  2. Likelihood chi² == analytic complex product on a synthetic Ae grid.
  3. Signed-Im: a sign flip in Be or Ae imag moves the residual (not folded).
  4. NaN / -1e30 rejection when the excitation label or Ae cache is missing.
  5. _load_be_excitation closest-period mapping against the real Europa file.
  6. SBIRunner._channel_conventions marks Bind_ 'signed', Im_k2 per-global.
  7. SBIRunner does not abs-fold a Bind_ input in _x_obs_vector.
"""
import os
import sys
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from PlanetProfile.Inference.mcmc_runner import (
    MCMCRunner,
    _parse_bind_channel,
    BE_PERIOD_MATCH_TOL_HR,
)


# ---------------------------------------------------------------------------
# 1. Channel-name grammar
# ---------------------------------------------------------------------------

def test_parse_bind_channel_simple():
    assert _parse_bind_channel('Bind_synodic_x_real') == ('synodic', 'x', 'real')
    assert _parse_bind_channel('Bind_synodic_z_imag') == ('synodic', 'z', 'imag')


def test_parse_bind_channel_multiword_label():
    # Labels may contain spaces and digits; the _<comp>_<part> suffix is
    # stripped from the right, so the whole remainder is the label.
    assert _parse_bind_channel('Bind_synodic 2nd_y_imag') == ('synodic 2nd', 'y', 'imag')
    assert _parse_bind_channel('Bind_adjusted orbital_x_real') == ('adjusted orbital', 'x', 'real')


def test_parse_bind_channel_rejects_non_bind():
    assert _parse_bind_channel('Ae_synodic_real') is None
    assert _parse_bind_channel('Re_k2') is None
    assert _parse_bind_channel('Bind_synodic_x') is None          # no part
    assert _parse_bind_channel('Bind_synodic_q_real') is None     # bad comp
    assert _parse_bind_channel('Bind__x_real') is None            # empty label


# ---------------------------------------------------------------------------
# Stub plumbing for the likelihood closure
# ---------------------------------------------------------------------------

class _StubConfig:
    def __init__(self, **kw):
        self.sampler_settings = kw.get('sampler_settings') or {}
        self.derived_params = kw.get('derived_params') or {}
        self.induction_bounds = kw.get('induction_bounds') or {}
        self.bodyname = kw.get('bodyname', 'Europa')


def _make_bind_ll(observables, ae_grid_cache, be_excitation, structure_data,
                  param_names=('Tb_K',)):
    """Build the log-likelihood with a stub carrying the induction caches.

    Mirrors MCMCRunner's attribute surface for the induction fast-path:
    self._ae_grid_cache (Tb-index -> {label: complex Ae}), self._be_excitation
    ({label: {'x'/'y'/'z': complex Be}}), and a grid structure_data.
    """
    class _Stub:
        pass

    stub = _Stub()
    stub.param_names = list(param_names)
    stub.param_groups = {}
    stub.fixed_params = {}
    stub.config = _StubConfig()
    # v5 reparam: the theta-expansion calls _inject_derived_Tb; these stubs
    # sample Tb_K directly, so bind the real method with the D_iceIh flag off
    # (it early-returns, a no-op for Tb-sampled configs).
    stub._samples_D_iceIh = False
    stub._inject_derived_Tb = MCMCRunner._inject_derived_Tb.__get__(stub, MCMCRunner)
    stub._gravity_clairaut_active = MCMCRunner._gravity_clairaut_active.__get__(stub, MCMCRunner)
    stub._ae_grid_cache = ae_grid_cache
    stub._be_excitation = be_excitation
    # The induction Ae lookup lives in MCMCRunner._blended_ae_dict (shared by
    # the likelihood and the SBI support guard; nearest-Tb for 1D caches,
    # bilinear (Tb, log w) for 2D v3.0 caches). Bind the real method + the
    # structure_data it reads so the stub exercises the production lookup.
    stub.structure_data = structure_data
    stub._blended_ae_dict = MCMCRunner._blended_ae_dict.__get__(stub, MCMCRunner)
    return MCMCRunner._make_flexible_log_likelihood(stub, observables, structure_data)


@pytest.fixture
def mock_k2(monkeypatch):
    """Constant k2/h2 + identity apply_parameters so induction is isolated."""
    import PlanetProfile.Inference.forward_models as fm

    def _mock_fwd(theta_dict, structure_data, **kw):
        # (Re_k2, Im_k2, Re_h2, Im_h2, heating)
        return 0.25, -0.05, 1.2, -0.1, None

    def _mock_apply(theta_dict, structure_data):
        # Return a dict with no phases/P/T so the no-ocean guard is inert.
        return dict(structure_data)

    monkeypatch.setattr(fm, 'forward_model_k2_flexible', _mock_fwd)
    monkeypatch.setattr(fm, 'apply_parameters', _mock_apply)
    yield


def _grid_sd(tb_grid):
    return {'Tb_K_grid': np.asarray(tb_grid, dtype=float),
            'structures': [{} for _ in tb_grid]}


# ---------------------------------------------------------------------------
# 2 & 3. Likelihood chi² == analytic complex product; signed Im
# ---------------------------------------------------------------------------

def test_bind_likelihood_matches_analytic_product(mock_k2):
    tb_grid = [260.0, 266.0, 271.0]
    # Synthetic Ae per grid point (complex, dimensionless).
    Ae_by_idx = {0: {'synodic': complex(0.75, -0.20)},
                 1: {'synodic': complex(0.90, -0.10)},
                 2: {'synodic': complex(0.94, -0.05)}}
    ae_grid = {i: Ae_by_idx[i] for i in range(3)}
    # Synthetic Be_x for synodic (complex, nT).
    Be = {'synodic': {'x': complex(128.5, -170.0),
                      'y': complex(-65.0, -38.0),
                      'z': complex(-4.8, -15.2)}}

    # Condition Bind_synodic_x_real / _imag exactly at the model value for
    # grid idx 1 (Tb=266) -> zero residual -> chi²=0 -> ll = 0.
    Ae = Ae_by_idx[1]['synodic']
    Bind_x = Ae * Be['synodic']['x']
    observables = {
        'Bind_synodic_x_real': (Bind_x.real, 1.5),
        'Bind_synodic_x_imag': (Bind_x.imag, 1.5),
    }
    ll = _make_bind_ll(observables, ae_grid, Be, _grid_sd(tb_grid))
    val = ll(np.array([266.0]))
    assert np.isfinite(val)
    assert abs(val) < 1e-9, f"expected ~0 ll at the exact model value, got {val}"

    # Offset the real channel by exactly 3 sigma -> chi² = 9 -> ll = -4.5.
    observables_off = dict(observables)
    observables_off['Bind_synodic_x_real'] = (Bind_x.real + 3 * 1.5, 1.5)
    ll_off = _make_bind_ll(observables_off, ae_grid, Be, _grid_sd(tb_grid))
    val_off = ll_off(np.array([266.0]))
    assert abs(val_off - (-4.5)) < 1e-9, f"expected -4.5, got {val_off}"


def test_bind_imag_is_signed_not_folded(mock_k2):
    """A signed Im residual must respond to the SIGN of the imaginary part;
    folding would collapse +Im and -Im onto the same |Im|."""
    tb_grid = [266.0]
    Ae = complex(0.9, -0.10)
    ae_grid = {0: {'synodic': Ae}}
    Be = {'synodic': {'x': complex(100.0, 0.0),
                      'y': complex(0.0, 0.0),
                      'z': complex(0.0, 0.0)}}
    Bind_x = Ae * Be['synodic']['x']   # imag = 0.9*0 + (-0.10)*100 = -10 nT

    # Condition on +|Im| (folded value). Signed model (-10) vs +10 -> residual
    # of 20/1.5, chi² large. If the code folded, residual would be 0.
    observables = {'Bind_synodic_x_imag': (abs(Bind_x.imag), 1.5)}
    ll = _make_bind_ll(observables, ae_grid, Be, _grid_sd(tb_grid))
    val = ll(np.array([266.0]))
    expected_chi2 = ((Bind_x.imag - abs(Bind_x.imag)) / 1.5) ** 2
    assert abs(val - (-0.5 * expected_chi2)) < 1e-9
    assert val < -1.0, "signed-Im residual collapsed to ~0 -> channel was folded"


# ---------------------------------------------------------------------------
# 4. Missing label / missing cache -> rejection
# ---------------------------------------------------------------------------

def test_bind_missing_be_label_rejected(mock_k2):
    tb_grid = [266.0]
    ae_grid = {0: {'synodic': complex(0.9, -0.1)}}
    Be = {}  # excitation for 'synodic' absent
    observables = {'Bind_synodic_x_real': (100.0, 1.5)}
    ll = _make_bind_ll(observables, ae_grid, Be, _grid_sd(tb_grid))
    assert ll(np.array([266.0])) == -1e30


def test_bind_missing_ae_label_rejected(mock_k2):
    tb_grid = [266.0]
    ae_grid = {0: {}}  # Ae for 'synodic' absent at this grid point
    Be = {'synodic': {'x': complex(100.0, 0.0), 'y': 0j, 'z': 0j}}
    observables = {'Bind_synodic_x_real': (100.0, 1.5)}
    ll = _make_bind_ll(observables, ae_grid, Be, _grid_sd(tb_grid))
    assert ll(np.array([266.0])) == -1e30


# ---------------------------------------------------------------------------
# 5. Be loader: closest-period mapping against the real Europa file
# ---------------------------------------------------------------------------

def test_load_be_excitation_europa_closest_period():
    """_load_be_excitation maps canonical labels to Be1xyz_Europa.txt rows by
    closest period, tolerating the 'orbital' (cache) vs 'adjusted orbital'
    (file) label drift, and returns the documented complex components."""
    class _Stub:
        pass
    stub = _Stub()
    stub.config = _StubConfig(bodyname='Europa')

    observables = {
        'Bind_synodic_x_real': (0.0, 1.5),
        'Bind_synodic 2nd_x_real': (0.0, 1.5),
        'Bind_orbital_x_real': (0.0, 1.5),
    }
    be = MCMCRunner._load_be_excitation(stub, observables)
    assert be is not None
    assert set(be) == {'synodic', 'synodic 2nd', 'orbital'}

    # synodic row: Bex = 128.4663 - 170.0841j (file row period 11.233 hr).
    assert abs(be['synodic']['x'].real - 128.466302323619) < 1e-6
    assert abs(be['synodic']['x'].imag - (-170.084146795136)) < 1e-6
    # |Bex| synodic ~ 213 nT (matches the plan's excitation-constants table).
    assert abs(abs(be['synodic']['x']) - 213.2) < 0.5

    # 'orbital' (cache label, 85.238 hr) must resolve to the 'adjusted
    # orbital' file row (85.213 hr), NOT to a true-anomaly row (~84.6 hr).
    assert abs(be['orbital']['x'].real - (-7.21158212931032)) < 1e-6
    assert abs(be['orbital']['x'].imag - 7.51709570390079) < 1e-6

    # synodic 2nd row (5.617 hr).
    assert abs(be['synodic 2nd']['x'].real - 16.8101001632006) < 1e-6


def test_load_be_excitation_none_without_bind():
    class _Stub:
        pass
    stub = _Stub()
    stub.config = _StubConfig(bodyname='Europa')
    assert MCMCRunner._load_be_excitation(stub, {'Re_k2': (0.6, 0.03)}) is None


def test_load_be_excitation_europa_all_labels_within_margin():
    """Every real Europa Bind_ label must map to a file row well inside the
    period-match tolerance (scientific review margin guard). If a future
    excitation-table edit collapses the orbital/true-anomaly separation this
    fails loudly rather than binding to the wrong excitation."""
    class _Stub:
        pass
    stub = _Stub()
    stub.config = _StubConfig(bodyname='Europa')
    observables = {f'Bind_{lbl}_x_real': (0.0, 1.5)
                   for lbl in ('synodic', 'synodic 2nd', 'orbital')}
    be = MCMCRunner._load_be_excitation(stub, observables)
    # All three must have survived the margin guard.
    assert set(be) == {'synodic', 'synodic 2nd', 'orbital'}
    assert BE_PERIOD_MATCH_TOL_HR <= 0.1


def test_load_be_excitation_margin_guard_rejects_far_label(monkeypatch):
    """A label whose period is far from any file row is refused (NaN channel),
    not silently bound to the nearest unrelated excitation."""
    from PlanetProfile.MagneticInduction.Moments import Excitations

    class _Stub:
        pass
    stub = _Stub()
    stub.config = _StubConfig(bodyname='Europa')

    # Inject a bogus label 1000 hr from every Europa file row.
    orig = dict(Excitations.Texc_hr['Europa'])
    patched = dict(orig)
    patched['synodic'] = 1000.0
    monkeypatch.setitem(Excitations.Texc_hr, 'Europa', patched)
    try:
        be = MCMCRunner._load_be_excitation(
            stub, {'Bind_synodic_x_real': (0.0, 1.5)})
    finally:
        Excitations.Texc_hr['Europa'] = orig
    # Only label requested was out-of-range -> no valid excitation -> None.
    assert be is None


# ---------------------------------------------------------------------------
# 6 & 7. SBIRunner channel-convention metadata + no-fold guarantee
# ---------------------------------------------------------------------------

def _sbi_stub(obs_names, imag_convention='abs'):
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    obj = SBIRunner.__new__(SBIRunner)
    obj.obs_names = list(obs_names)
    obj.imag_convention = imag_convention
    return obj


def test_channel_conventions_marks_bind_signed():
    obj = _sbi_stub(['Re_k2', 'abs_Im_k2', 'CMR2',
                     'Bind_synodic_x_real', 'Bind_synodic_x_imag'])
    conv = obj._channel_conventions()
    assert conv['Bind_synodic_x_real'] == 'signed'
    assert conv['Bind_synodic_x_imag'] == 'signed'
    assert conv['abs_Im_k2'] == 'abs'
    # Purely-real / magnitude channels carry no Im ambiguity -> omitted.
    assert 'Re_k2' not in conv
    assert 'CMR2' not in conv


def test_x_obs_vector_does_not_fold_bind():
    obj = _sbi_stub(['abs_Im_k2', 'Bind_synodic_x_imag'])
    x_obs = {'abs_Im_k2': -0.135, 'Bind_synodic_x_imag': -42.0}
    vec = obj._x_obs_vector(x_obs)
    # Im_k2 folded to +0.135; Bind imag preserved signed at -42.
    assert abs(vec[0] - 0.135) < 1e-12
    assert abs(vec[1] - (-42.0)) < 1e-12


# ---------------------------------------------------------------------------
# 8. Per-frequency, degree-based induction support cuts (v5 open-interp pivot)
#
# The v5 "open interpretation" redesign replaces the single global Ae guard
# (amp_min 0.7 / im_abs_max 0.4 on the synodic band) with per-label bounds and
# a *degree-based* phase-delay cap (`phase_deg_max`). The literature basis
# (Vance et al. 2021, Fig 6 / Table 2) is that the tight synodic constraint
# (|Ae|>0.70, phase<30 deg) does NOT hold at the longer-period orbital band:
# low-salinity thick-ice models reach |Ae_orbital| ~= 0.37 at a ~65 deg phase
# delay. Those samples MUST be admitted by the orbital bound and rejected by
# the old global bound. The |Im| proxy (im_abs_max) mismaps to phase at low
# |Ae|; the degree cap is amplitude-independent and is the correct knob.
#
# Both the MCMC likelihood (mcmc_runner ~line 1003) and the SBI support guard
# (_check_induction_bounds ~line 2289) apply the SAME predicate, so the
# training support == reference-MCMC support. The guard is a closure inside
# generate_sbi_dataset (not directly callable), so these tests exercise the
# likelihood site (the code path the reference MCMC actually runs) plus a
# standalone replica of the shared predicate to assert the two agree.
# ---------------------------------------------------------------------------

def _make_bound_ll(induction_bounds, ae_by_label, tb_grid=(266.0,),
                   observables=None):
    """Log-likelihood stub with `induction_bounds` configured on the config.

    No Ae_/Bind_ observables by default -> chi2 stays 0, so a bound-passing
    draw yields ll == 0 and a bound-violating draw yields -1e30. This isolates
    the support cut from any residual term.
    """
    class _Stub:
        pass
    stub = _Stub()
    stub.param_names = ['Tb_K']
    stub.param_groups = {}
    stub.fixed_params = {}
    stub.config = _StubConfig(induction_bounds=induction_bounds)
    stub._samples_D_iceIh = False
    stub._inject_derived_Tb = MCMCRunner._inject_derived_Tb.__get__(stub, MCMCRunner)
    # Clairaut gravity is inactive for this stub config -> real predicate
    # returns False, so the k2/gravity block is skipped and the closure reaches
    # the induction support cut we are exercising.
    stub._gravity_clairaut_active = MCMCRunner._gravity_clairaut_active.__get__(stub, MCMCRunner)
    # One Ae dict, reused at every (single) grid node.
    ae_grid = {i: dict(ae_by_label) for i in range(len(tb_grid))}
    stub._ae_grid_cache = ae_grid
    stub._be_excitation = {}
    stub.structure_data = _grid_sd(tb_grid)
    stub._blended_ae_dict = MCMCRunner._blended_ae_dict.__get__(stub, MCMCRunner)
    ll = MCMCRunner._make_flexible_log_likelihood(
        stub, observables or {}, stub.structure_data)
    return ll


def _predicate_violates(spec, Ae):
    """Standalone replica of the shared bound predicate (both call sites).

    Kept byte-faithful to mcmc_runner's amp_min / im_abs_max / phase_deg_max
    checks so the test fails loudly if the two production sites ever drift from
    this contract. Returns True when the sample must be rejected.
    """
    import numpy as _np
    Ae = complex(Ae)
    amp_min = spec.get('amp_min')
    if amp_min is not None and abs(Ae) < float(amp_min):
        return True
    im_abs_max = spec.get('im_abs_max')
    if im_abs_max is not None and abs(Ae.imag) > float(im_abs_max):
        return True
    phase_deg_max = spec.get('phase_deg_max')
    if (phase_deg_max is not None
            and abs(float(_np.degrees(_np.angle(Ae)))) > float(phase_deg_max)):
        return True
    return False


# Vance et al. 2021 low-salinity thick-ice orbital point: |Ae|~=0.374 at a
# ~65.1 deg phase delay (Ae^{-i phi} convention -> Im < 0).
_AE_ORBITAL_THICK = 0.374 * np.exp(-1j * np.deg2rad(65.1))
# A representative synodic point (well inside the tight synodic bound).
_AE_SYNODIC = 0.853 * np.exp(-1j * np.deg2rad(16.0))

_V5_ORBITAL_BOUND = {'amp_min': 0.3, 'phase_deg_max': 70.0}
_V5_SYNODIC_BOUND = {'amp_min': 0.7, 'phase_deg_max': 30.0}
# The retired pre-pivot global guard (synodic-tuned, |Im| proxy).
_OLD_GLOBAL_BOUND = {'amp_min': 0.7, 'im_abs_max': 0.4}


def test_orbital_thick_ice_accepted_under_v5_bound(mock_k2):
    """|Ae_orbital|~=0.37 @ 65 deg is ACCEPTED by the per-freq orbital bound."""
    assert not _predicate_violates(_V5_ORBITAL_BOUND, _AE_ORBITAL_THICK)
    ll = _make_bound_ll({'orbital': _V5_ORBITAL_BOUND},
                        {'orbital': _AE_ORBITAL_THICK})
    val = ll(np.array([266.0]))
    assert val == 0.0, f"expected ll=0 (accepted, no residual terms), got {val}"


def test_orbital_thick_ice_rejected_under_old_global_bound(mock_k2):
    """The SAME point is REJECTED by the retired global 0.7/0.4 guard
    (|Ae|=0.37 < 0.7). This is exactly the low-salinity support the pivot
    restores."""
    assert _predicate_violates(_OLD_GLOBAL_BOUND, _AE_ORBITAL_THICK)
    ll = _make_bound_ll({'orbital': _OLD_GLOBAL_BOUND},
                        {'orbital': _AE_ORBITAL_THICK})
    assert ll(np.array([266.0])) == -1e30


def test_synodic_tight_point_accepted_under_v5_bound(mock_k2):
    """A canonical synodic point (|Ae|=0.85 @ 16 deg) passes the tight
    synodic bound the pivot keeps."""
    assert not _predicate_violates(_V5_SYNODIC_BOUND, _AE_SYNODIC)
    ll = _make_bound_ll({'synodic': _V5_SYNODIC_BOUND},
                        {'synodic': _AE_SYNODIC})
    assert ll(np.array([266.0])) == 0.0


def test_orbital_thick_ice_rejected_by_synodic_bound(mock_k2):
    """Applying the tight SYNODIC bound to the orbital point rejects it
    (|Ae|=0.37 < 0.7 and phase 65 > 30) — i.e. a single global bound cannot
    serve both bands, which is the whole motivation for per-freq cuts."""
    assert _predicate_violates(_V5_SYNODIC_BOUND, _AE_ORBITAL_THICK)
    ll = _make_bound_ll({'orbital': _V5_SYNODIC_BOUND},
                        {'orbital': _AE_ORBITAL_THICK})
    assert ll(np.array([266.0])) == -1e30


def test_phase_cap_operates_in_degrees_not_radians(mock_k2):
    """The phase cap is applied to |angle(Ae)| in DEGREES. A 65 deg point sits
    just under a 70 deg cap (accept) and just over a 60 deg cap (reject). If
    the code compared radians (|angle|=1.14 rad), a 70 'deg' cap would never
    bite and the reject case below would spuriously pass."""
    ae = 0.374 * np.exp(-1j * np.deg2rad(65.1))
    # 70 deg cap -> accept
    assert not _predicate_violates({'amp_min': 0.3, 'phase_deg_max': 70.0}, ae)
    assert _make_bound_ll({'orbital': {'amp_min': 0.3, 'phase_deg_max': 70.0}},
                          {'orbital': ae})(np.array([266.0])) == 0.0
    # 60 deg cap -> reject (65.1 > 60)
    assert _predicate_violates({'amp_min': 0.3, 'phase_deg_max': 60.0}, ae)
    assert _make_bound_ll({'orbital': {'amp_min': 0.3, 'phase_deg_max': 60.0}},
                          {'orbital': ae})(np.array([266.0])) == -1e30


def test_likelihood_and_guard_predicate_agree(mock_k2):
    """The likelihood-site accept/reject must match the standalone predicate
    (which is byte-faithful to the SBI support-guard site) across a grid of
    (Ae, bound) cases — the reviewer-binding 'training support == likelihood
    support' contract."""
    ae_cases = [
        _AE_ORBITAL_THICK,
        _AE_SYNODIC,
        0.374 * np.exp(-1j * np.deg2rad(65.1)),
        0.20 * np.exp(-1j * np.deg2rad(80.0)),   # below every amp floor
        0.95 * np.exp(-1j * np.deg2rad(5.0)),    # strong, near-in-phase
    ]
    bound_cases = [_V5_ORBITAL_BOUND, _V5_SYNODIC_BOUND, _OLD_GLOBAL_BOUND]
    for ae in ae_cases:
        for spec in bound_cases:
            pred_reject = _predicate_violates(spec, ae)
            ll = _make_bound_ll({'orbital': spec}, {'orbital': ae})
            ll_reject = (ll(np.array([266.0])) == -1e30)
            assert pred_reject == ll_reject, (
                f"disagreement: Ae={ae:.3f} spec={spec} "
                f"predicate_reject={pred_reject} likelihood_reject={ll_reject}")
