"""Externally prescribed temperature profiles for solid (silicate) layers.

Provides a body- and dataset-agnostic way to supply T(depth) from outside
PlanetProfile, for cases where PlanetProfile's own conductive/adiabatic
silicate propagation is not trusted for the body in question.

Only the temperature assignment is replaced. The structural recursion --
pressure, density, gravity, layer mass, moment of inertia and the core-radius
search -- is untouched, and the EOS is still evaluated at the prescribed
temperature. A model built this way is hydrostatically self-consistent
*conditional on* the prescribed thermal field; it is NOT thermally
self-consistent.

Silent extrapolation is deliberately impossible: a query outside the supplied
depth range raises.
"""
import logging

import numpy as np

log = logging.getLogger('PlanetProfile')


class PrescribedTemperatureRangeError(ValueError):
    """Raised when a prescribed-T profile is queried outside its data range."""
    pass


class PrescribedTemperatureProfile:
    """Callable T_K = f(depth_m) built from a tabulated profile.

    Parameters
    ----------
    depth_m, T_K : array_like
        Tabulated profile, strictly increasing in depth.
    source : str
        Provenance string, carried into model output.
    method : {'linear'}
        Interpolation scheme. Linear is the default and is deliberately the
        only option: it is shape-preserving, introduces no overshoot, and
        cannot manufacture structure between tabulated levels. A spline is
        NOT offered, because splines ring across sharp thermal gradients and
        would interpolate across physical boundaries.
    """

    def __init__(self, depth_m, T_K, source='unspecified', method='linear',
                 deepContinuation=None):
        depth_m = np.asarray(depth_m, dtype=float)
        T_K = np.asarray(T_K, dtype=float)
        if depth_m.ndim != 1 or T_K.ndim != 1:
            raise ValueError('Prescribed T profile must be 1-D.')
        if depth_m.size != T_K.size:
            raise ValueError('depth_m and T_K must have the same length.')
        if depth_m.size < 2:
            raise ValueError('Prescribed T profile needs at least two levels.')
        if not np.all(np.isfinite(depth_m)) or not np.all(np.isfinite(T_K)):
            raise ValueError('Prescribed T profile contains non-finite values.')
        if not np.all(np.diff(depth_m) > 0):
            raise ValueError('Prescribed T profile depths must be strictly increasing.')
        if method != 'linear':
            raise ValueError(f'Unsupported interpolation method {method!r}; only "linear" is supported.')

        self.depth_m = depth_m
        self.T_K = T_K
        self.source = source
        self.method = method
        self.zMin_m = float(depth_m[0])
        self.zMax_m = float(depth_m[-1])
        self.TMin_K = float(np.min(T_K))
        self.TMax_K = float(np.max(T_K))

        # Deep continuation policy.
        #
        # PlanetProfile propagates EVERY candidate silicate profile from its
        # start radius down to r = 0 (Silicates.py: rSilEnd_m = 0 when
        # Fe_CORE is set); the core only replaces the lower part afterwards.
        # No finite thermal profile can cover r -> 0, so strict mode is
        # architecturally impossible for a cored body.
        #
        #   None (default)  -- raise below the data. Correct for any use that
        #                      does not hit PlanetProfile's candidate scratch
        #                      space, and the safe default for other bodies.
        #   'isothermal'    -- hold the deepest tabulated temperature constant
        #                      below the data, and RECORD every such use.
        #                      Legitimate only if the final selected model is
        #                      afterwards shown not to depend on it; the
        #                      counters below exist to make that checkable.
        if deepContinuation not in (None, 'isothermal', 'linear'):
            raise ValueError(f'Unsupported deepContinuation {deepContinuation!r}; '
                             'expected None, "isothermal" or "linear".')
        self.deepContinuation = deepContinuation
        self.nDeepContinuationCalls = 0
        self.maxContinuationDepth_m = None
        # Deep gradient for 'linear' continuation.
        #
        # Defined by an unweighted least-squares fit to the deepest
        # DEEP_GRADIENT_NODES source nodes, rather than the final pair, so a
        # single noisy interval cannot set the slope. The node count is fixed
        # and documented so the rule is reproducible across compositions.
        self.DEEP_GRADIENT_NODES = 5
        self.deepGradient_KperM = None
        if deepContinuation == 'linear':
            n = min(self.DEEP_GRADIENT_NODES, depth_m.size)
            zz, TT = depth_m[-n:], T_K[-n:]
            self.deepGradient_KperM = float(np.polyfit(zz, TT, 1)[0])
            self.deepGradientNodesUsed = int(n)

    def __repr__(self):
        return (f'PrescribedTemperatureProfile({self.depth_m.size} levels, '
                f'z {self.zMin_m/1e3:.1f}-{self.zMax_m/1e3:.1f} km, '
                f'T {self.TMin_K:.1f}-{self.TMax_K:.1f} K, source={self.source!r})')

    def CheckRange(self, depth_m, label=''):
        """Raise if any requested depth lies outside the tabulated range.

        Depths deeper than the data are tolerated only when an explicit
        deepContinuation policy is set, and every such use is recorded.
        """
        z = np.atleast_1d(np.asarray(depth_m, dtype=float))
        below = z < self.zMin_m
        above = z > self.zMax_m
        if np.any(above) and self.deepContinuation is not None:
            self.nDeepContinuationCalls += int(np.count_nonzero(above))
            deepest = float(np.max(z[above]))
            if self.maxContinuationDepth_m is None or deepest > self.maxContinuationDepth_m:
                self.maxContinuationDepth_m = deepest
            above = np.zeros_like(above)
        if not (np.any(below) or np.any(above)):
            return
        parts = [f'Prescribed temperature profile queried outside its data range'
                 f'{" for " + label if label else ""}.',
                 f'  profile:   z {self.zMin_m/1e3:.3f} to {self.zMax_m/1e3:.3f} km '
                 f'({self.depth_m.size} levels, source: {self.source})',
                 f'  requested: z {np.min(z)/1e3:.3f} to {np.max(z)/1e3:.3f} km']
        if np.any(below):
            parts.append(f'  {int(np.count_nonzero(below))} node(s) shallower than the profile, '
                         f'by up to {(self.zMin_m - np.min(z[below]))/1e3:.3f} km')
        if np.any(above):
            parts.append(f'  {int(np.count_nonzero(above))} node(s) deeper than the profile, '
                         f'by up to {(np.max(z[above]) - self.zMax_m)/1e3:.3f} km')
        parts.append('  Extrapolation is not permitted. Either restrict the model domain to '
                     'the profile, or supply a profile covering it.')
        raise PrescribedTemperatureRangeError('\n'.join(parts))

    def __call__(self, depth_m, label=''):
        z = np.asarray(depth_m, dtype=float)
        self.CheckRange(z, label=label)
        # np.interp clamps to the end values outside the data range, which is
        # exactly the 'isothermal' continuation. CheckRange has already either
        # rejected out-of-range depths or recorded the continuation, so no
        # extrapolation can occur silently here.
        T = np.interp(z, self.depth_m, self.T_K)
        if self.deepContinuation == 'linear' and self.deepGradient_KperM is not None:
            deep = np.asarray(z) > self.zMax_m
            if np.any(deep):
                T = np.array(T, dtype=float, copy=True)
                T[deep] = self.T_K[-1] + self.deepGradient_KperM * (np.asarray(z)[deep] - self.zMax_m)
        return T


def LoadPrescribedTProfileCSV(fpath, source=None, deepContinuation=None):
    """Load a two-column ``depth_km,T_K`` CSV with ``#`` comment lines."""
    rows = []
    with open(fpath) as fh:
        for line in fh:
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            try:
                rows.append([float(v) for v in line.split(',')])
            except ValueError:
                # A single non-numeric column-header row is tolerated; anything
                # else is a malformed file and should surface as an error below.
                continue
    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or data.size == 0:
        raise ValueError(f'No numeric data found in {fpath}.')
    if data.shape[1] != 2:
        raise ValueError(f'Expected two columns (depth_km, T_K) in {fpath}, got {data.shape[1]}.')
    return PrescribedTemperatureProfile(data[:, 0] * 1e3, data[:, 1],
                                        source=source if source is not None else str(fpath),
                                        deepContinuation=deepContinuation)
