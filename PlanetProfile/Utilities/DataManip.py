"""
DataManip: Miscellaneous functions for adjusting data that don't belong in
a more specific location.
"""

import numpy as np
from PlanetProfile.Utilities.defineStructs import EOSlist
from PlanetProfile.Utilities.defineStructs import Constants
import logging
from scipy.interpolate import griddata

# Assign logger
log = logging.getLogger('PlanetProfile')

def ResetNearestExtrap(var1, var2, min1, max1, min2, max2):
    """ Choose the nearest value when EOS function inputs are outside
        the generation domain.

        var1 (float, shape N): First variable, typically P_MPa.
        var2 (float, shape M): First variable, typically T_K.
        min1, max1, min2, max2 (float): Domain boundaries.
    """
    outVar1 = var1 + 0.0
    outVar2 = var2 + 0.0
    if np.size(var1) == 1:
        outVar1 = np.array(var1)
    if np.size(var2) == 1:
        outVar2 = np.array(var2)

    outVar1[var1 < min1] = min1
    outVar1[var1 > max1] = max1
    outVar2[var2 < min2] = min2
    outVar2[var2 > max2] = max2

    return outVar1, outVar2

class EOSDomainError(ValueError):
    """Raised when an EOS is queried outside its tabulated P-T domain under a
    strict domain policy."""
    pass


# Recognized EOS-domain policies. See EnforceEOSDomain.
EOS_DOMAIN_POLICIES = ('legacy_clamp', 'warn', 'error')


def EOSDomainExceedance(var1, var2, min1, max1, min2, max2):
    """Quantify how far a set of (P, T) queries falls outside an EOS domain.

    Returns a dict of diagnostics; all counts are ints and all exceedances are
    in the same units as the inputs. An in-domain query set yields
    nOutside == 0 and zero exceedances.
    """
    v1 = np.atleast_1d(np.asarray(var1, dtype=float))
    v2 = np.atleast_1d(np.asarray(var2, dtype=float))
    # Broadcast so paired (P, T) queries are counted per node, not per axis.
    try:
        b1, b2 = np.broadcast_arrays(v1, v2)
    except ValueError:
        b1, b2 = v1, v2

    low1, high1 = b1 < min1, b1 > max1
    low2, high2 = b2 < min2, b2 > max2
    outside = low1 | high1 | low2 | high2

    def _worst(arr, mask, bound, sign):
        if not np.any(mask):
            return 0.0
        return float(np.max(sign * (arr[mask] - bound)))

    return {
        'nTotal': int(np.size(outside)),
        'nOutside': int(np.count_nonzero(outside)),
        'n1Low': int(np.count_nonzero(low1)), 'n1High': int(np.count_nonzero(high1)),
        'n2Low': int(np.count_nonzero(low2)), 'n2High': int(np.count_nonzero(high2)),
        'max1Below': _worst(b1, low1, min1, -1.0), 'max1Above': _worst(b1, high1, max1, 1.0),
        'max2Below': _worst(b2, low2, min2, -1.0), 'max2Above': _worst(b2, high2, max2, 1.0),
        'min1Requested': float(np.min(b1)) if np.size(b1) else np.nan,
        'max1Requested': float(np.max(b1)) if np.size(b1) else np.nan,
        'min2Requested': float(np.min(b2)) if np.size(b2) else np.nan,
        'max2Requested': float(np.max(b2)) if np.size(b2) else np.nan,
    }


def FormatEOSDomainMessage(report, min1, max1, min2, max2, EOSname='EOS',
                           var1name='P_MPa', var2name='T_K'):
    """Human-readable description of an EOS-domain violation."""
    parts = [
        f'{EOSname}: {report["nOutside"]} of {report["nTotal"]} node(s) requested '
        f'outside the tabulated EOS domain.',
        f'  domain:    {var1name} [{min1:.6g}, {max1:.6g}], {var2name} [{min2:.6g}, {max2:.6g}]',
        f'  requested: {var1name} [{report["min1Requested"]:.6g}, {report["max1Requested"]:.6g}], '
        f'{var2name} [{report["min2Requested"]:.6g}, {report["max2Requested"]:.6g}]',
    ]
    if report['n1Low']:
        parts.append(f'  {var1name} below min: {report["n1Low"]} node(s), max exceedance {report["max1Below"]:.6g}')
    if report['n1High']:
        parts.append(f'  {var1name} above max: {report["n1High"]} node(s), max exceedance {report["max1Above"]:.6g}')
    if report['n2Low']:
        parts.append(f'  {var2name} below min: {report["n2Low"]} node(s), max exceedance {report["max2Below"]:.6g}')
    if report['n2High']:
        parts.append(f'  {var2name} above max: {report["n2High"]} node(s), max exceedance {report["max2Above"]:.6g}')
    return '\n'.join(parts)


def EnforceEOSDomain(var1, var2, min1, max1, min2, max2, policy='legacy_clamp',
                     CLAMP=True, EOSname='EOS', var1name='P_MPa', var2name='T_K'):
    """Apply an EOS-domain policy to a set of (P, T) queries.

    Policies:
      'legacy_clamp' -- historical behavior: silently clamp to the nearest
                        domain boundary (only when CLAMP is True). No warning.
      'warn'         -- identical numerics to 'legacy_clamp', but log an
                        explicit warning naming the EOS, the requested and
                        tabulated ranges, the number of affected nodes and the
                        maximum exceedance.
      'error'        -- raise EOSDomainError before any property is evaluated.

    CLAMP mirrors ``not EOS.EXTRAP``: when extrapolation is permitted the
    values are passed through unchanged, but 'warn'/'error' still apply, so a
    deliberate extrapolation is never silent.

    Returns (var1, var2), clamped or unchanged per the policy.
    """
    if policy not in EOS_DOMAIN_POLICIES:
        raise ValueError(f'Unrecognized EOS domain policy {policy!r}. '
                         f'Expected one of {EOS_DOMAIN_POLICIES}.')

    if policy == 'legacy_clamp':
        # Fast path: preserve historical behavior exactly, with no extra work.
        if CLAMP:
            return ResetNearestExtrap(var1, var2, min1, max1, min2, max2)
        return var1, var2

    report = EOSDomainExceedance(var1, var2, min1, max1, min2, max2)
    if report['nOutside'] > 0:
        msg = FormatEOSDomainMessage(report, min1, max1, min2, max2,
                                     EOSname=EOSname, var1name=var1name, var2name=var2name)
        if policy == 'error':
            raise EOSDomainError(msg)
        log.warning(msg)

    if CLAMP:
        return ResetNearestExtrap(var1, var2, min1, max1, min2, max2)
    return var1, var2


def ReAssignPT(P_MPa, T_K, Pmin, Pmax, Tmin, Tmax, MELT=False, propsStepReductionFactor=1):
    # Find indices where P_MPa is within bounds [Pmin, Pmax]
    deltaP = np.mean(np.diff(P_MPa))
    deltaT = np.mean(np.diff(T_K))
    P_MPa = np.arange(Pmin, Pmax+deltaP, deltaP)
    T_K = np.arange(Tmin, Tmax+deltaT, deltaT)
    if MELT:
        PropsP_MPa = np.linspace(P_MPa[0], P_MPa[-1], 10)
        PropsT_K = np.linspace(T_K[0], T_K[-1], 11)
        Pphase_MPa = P_MPa
        Tphase_K = T_K
    else:
        PropsP_MPa = np.linspace(P_MPa[0], P_MPa[-1], int(np.ceil(np.size(P_MPa)/propsStepReductionFactor)))
        PropsT_K = np.linspace(T_K[0], T_K[-1], int(np.ceil(np.size(T_K)/propsStepReductionFactor)))
        Pphase_MPa = P_MPa
        Tphase_K = T_K
    return PropsP_MPa, PropsT_K, Pphase_MPa, Tphase_K


def smoothGrid(x, y, zs, factor):
    """ Smooth a grid by a factor """
    # Create finer resolution grid for smoothing
    x_min, x_max = np.nanmin(x), np.nanmax(x)
    y_min, y_max = np.nanmin(y), np.nanmax(y)
    
    # Ensure x and y are above machine error to prevent interpolation errors
    if (abs(x_min) < 1e-10 and abs(x_max) < 1e-10) or (abs(y_min) < 1e-10 and abs(y_max) < 1e-10):
        log.debug('Attempting to smooth 2d grid for plotting, but x or y is below machine error. Smoothing will not be performed.')
        x_fine = x
        y_fine = y
    else: 
        # Get original grid dimensions
        original_x_size = x.shape[0]
        original_y_size = x.shape[1]
        
        # Create finer grid with higher resolution
        n_x_fine = original_x_size * factor
        n_y_fine = original_y_size * factor
        
        # Create evenly spaced, strictly ascending fine grids
        x_fine_1d = np.linspace(x_min, x_max, n_x_fine)
        y_fine_1d = np.linspace(y_min, y_max, n_y_fine)
        x_fine, y_fine = np.meshgrid(x_fine_1d, y_fine_1d, indexing='ij')
        
        # Flatten original data for interpolation
        x_flat = x.flatten()
        y_flat = y.flatten()
        for i, z in enumerate(zs):
            z_flat = z.flatten()
            
            # Remove invalid points for interpolation
            valid_mask = ~np.isnan(z_flat)
            x_valid = x_flat[valid_mask]
            y_valid = y_flat[valid_mask]
            z_valid = z_flat[valid_mask]
            if len(z_valid) > 0:
                try:
                    z_fine = griddata((x_valid, y_valid), z_valid, 
                                    (x_fine, y_fine), method='linear', fill_value=np.nan)
                    # Update x, y, z to use the finer resolution data
                    zs[i] = z_fine
                except Exception as e:
                    log.debug(f'Error smoothing 2d grid for plotting: {e}')
                    zs[i] = z
            else:
                zs[i] = np.full((x_fine.shape[0], y_fine.shape[1]), np.nan)
    return x_fine, y_fine, zs
class ReturnZeros:
    """ Returns an array or tuple of arrays of zeros, for functions of properties
        not modeled that still work with querying routines. We have to run things
        this way and not with a lambda because anonymous functions can't be used
        in parallel processing.
    """
    def __init__(self, nVar):
        self.nVar = nVar

    def __call__(self, P, T, grid=False):
        if grid:
            out = np.zeros((np.size(P), np.size(T)))
        else:
            nPs = np.size(P)
            nTs = np.size(T)
            if self.nVar > 1:
                out = tuple(np.zeros(np.maximum(nPs, nTs)) for _ in range(self.nVar))
            else:
                out = np.zeros(np.maximum(nPs, nTs))
        return out


class ReturnConstantPTw:
    """ Returns an array or tuple of arrays of a constant value, for functions of
        properties not modeled that still work with querying routines. We have to
        run things this way and not with a lambda because anonymous functions
        can't be used in parallel processing.
    """
    def __init__(self, const=None, nVar=1):
        if const is None:
            self.const = 0
        else:
            self.const = const
        self.nVar = nVar

    def __call__(self, P, T, w, grid=False, KEEP_wDIM=False):
        nPs = np.size(P)
        nTs = np.size(T)
        nws = np.size(w)
        if grid or nPs != nTs:
            if KEEP_wDIM:
                out = np.zeros((nPs, nTs, nws)) + self.const
            else:
                out = np.squeeze(np.zeros((nPs, nTs, nws)) + self.const)
        else:
            if self.nVar > 1:
                out = (np.zeros(np.max([nPs, nTs, nws])) for _ in range(self.nVar))
            else:
                out = np.zeros(np.max([nPs, nTs, nws]))
        return out

class ReturnConstantSpecies:
    """ Returns 3 arrays of constant values, for functions of
        properties not modeled that still work with querying routines. We have to
        run things this way and not with a lambda because anonymous functions
        can't be used in parallel processing.
        reference_ppt: Only applicable to Seawater. Gives the reference ppt of the speciation in the associated
        Standard Seawater dictionary and its species
    """
    def __init__(self, ppt, reference_ppt, pH, species):
        if reference_ppt is None:
            speciation_factor = ppt
        else:
            speciation_factor = ppt/reference_ppt # Relative ratio between ppt of HydroEOS and reference ppt in
            # Constant dictionary
        if pH is None:
            self.pH = np.nan
        else:
            self.pH = pH
        if species is None:
            self.species_names = np.array([])
            self.speciation = np.array([])
        else:
            self.species_names = np.array(list(species.keys()))
            self.speciation = np.array(list(species.values())) * speciation_factor
        # Add 55.5mol/kg of water to speciation and species_name
        self.species_names = np.append(self.species_names,'H2O(aq)')
        self.speciation = np.append(self.speciation, 1/Constants.m_gmol['H2O']*1000)

    def __call__(self, P, T, grid = False, reactionEquation = None):
        nPs = np.size(P)
        nTs = np.size(T)
        pH = (np.zeros(nPs)) + self.pH

        speciation = np.zeros((len(self.species_names), nPs))

        # Populate the speciation array
        for row_index, value in enumerate(self.speciation):
            speciation[row_index] = value  # Fill the entire row with the corresponding speciation value
        species_names = np.array(self.species_names)
        return pH, speciation, species_names

class Nearest2DInterpolator:
    """Simplified version that may be faster for many use cases.
    We use our own 2D interpolator to reduce memory overhead by allowing x, y, and z to keep their data types, rather than upcasting to float64 like scipy does.
    This is particularly relevant for the phase lookup tables, where we can store the phase lookup table as a uint8 array."""
    def __init__(self, x, y, z):
        self.x = np.asarray(x)
        self.y = np.asarray(y) 
        self.z = np.asarray(z)
        self.nx = len(x)
        self.ny = len(y)
    
    def __call__(self, xi, yi, grid=False):
        if grid:
            return self._interpolate_grid(xi, yi)
        else:
            return self._interpolate(xi, yi)
    
    def _interpolate_grid(self, xi, yi):
        """Fast grid interpolation without creating intermediate meshgrids"""
        xi = np.asarray(xi)
        yi = np.asarray(yi)
        
        # Find nearest indices for each dimension independently
        ix = np.searchsorted(self.x, xi)
        iy = np.searchsorted(self.y, yi)
        
        # Handle boundaries
        ix_left = ix == 0
        iy_left = iy == 0
        
        ix = np.minimum(ix, self.nx - 1)
        iy = np.minimum(iy, self.ny - 1)
        
        # Distance comparisons
        ix_adj = ~ix_left & ((xi - self.x[ix-1]) < (self.x[ix] - xi))
        iy_adj = ~iy_left & ((yi - self.y[iy-1]) < (self.y[iy] - yi))
        
        ix -= ix_adj.astype(np.intp)
        iy -= iy_adj.astype(np.intp)
        
        # Use advanced indexing to get grid result directly
        # This creates the meshgrid effect without storing intermediate arrays
        return self.z[np.ix_(ix, iy)]
    
    def _interpolate(self, xi, yi):
        xi = np.asarray(xi)
        yi = np.asarray(yi)
        
        # Binary search
        ix = np.searchsorted(self.x, xi)
        iy = np.searchsorted(self.y, yi)
        
        # Boundary handling with single operations
        ix_left = ix == 0
        iy_left = iy == 0
        
        # Clamp and adjust in one step
        ix = np.minimum(ix, self.nx - 1)
        iy = np.minimum(iy, self.ny - 1)
        
        # Distance comparison (avoid where ix/iy is 0)
        ix_adj = ~ix_left & ((xi - self.x[ix-1]) < (self.x[ix] - xi))
        iy_adj = ~iy_left & ((yi - self.y[iy-1]) < (self.y[iy] - yi))
        
        ix -= ix_adj.astype(np.intp)
        iy -= iy_adj.astype(np.intp)
        
        return self.z[ix, iy]


class EOSwrapper:
    """ Lightweight wrapper for accessing EOS functions stored in the EOSlist dict. """

    def __init__(self, key):
        self.key = key

        # Assign only those attributes we reference in functions
        if EOSlist.loaded[self.key].EOStype == 'ice':
            self.phaseID = EOSlist.loaded[self.key].phaseID
            self.POROUS = EOSlist.loaded[self.key].POROUS
            self.Tmin = EOSlist.loaded[self.key].Tmin
        elif EOSlist.loaded[self.key].EOStype == 'ocean':
            self.Pmin = EOSlist.loaded[self.key].Pmin
            self.Pmax = EOSlist.loaded[self.key].Pmax
            self.propsPmax = EOSlist.loaded[self.key].propsPmax
            self.Tmin = EOSlist.loaded[self.key].Tmin
            self.Tmax = EOSlist.loaded[self.key].Tmax
            self.deltaP = EOSlist.loaded[self.key].deltaP
            self.deltaT = EOSlist.loaded[self.key].deltaT
            self.EOSdeltaP = EOSlist.loaded[self.key].EOSdeltaP
            self.EOSdeltaT = EOSlist.loaded[self.key].EOSdeltaT
            self.comp = EOSlist.loaded[self.key].comp
            self.w_ppt = EOSlist.loaded[self.key].w_ppt
        elif EOSlist.loaded[self.key].EOStype == 'inner':
            self.comp = EOSlist.loaded[self.key].comp

    def fn_phase(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_phase(P_MPa, T_K, grid=grid)
    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_rho_kgm3(P_MPa, T_K, grid=grid)
    def fn_Cp_JkgK(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_Cp_JkgK(P_MPa, T_K, grid=grid)
    def fn_alpha_pK(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_alpha_pK(P_MPa, T_K, grid=grid)
    def fn_kTherm_WmK(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_kTherm_WmK(P_MPa, T_K, grid=grid)
    def fn_phi_frac(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_phi_frac(P_MPa, T_K, grid=grid)
    def fn_porosCorrect(self, propBulk, propPore, phi, J):
        return EOSlist.loaded[self.key].fn_porosCorrect(propBulk, propPore, phi, J)
    def fn_sigma_Sm(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_sigma_Sm(P_MPa, T_K, grid=grid)
    def fn_eta_Pas(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_eta_Pas(P_MPa, T_K, grid=grid)
    def updateConvectionViscosity(self, etaConv_Pas, Tconv_K):
        return EOSlist.loaded[self.key].updateConvectionViscosity(etaConv_Pas, Tconv_K)
    def fn_Seismic(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_Seismic(P_MPa, T_K, grid=grid)
    def fn_species(self, P_MPa, T_K, grid = False, reactionEquation=None):
        return EOSlist.loaded[self.key].fn_species(P_MPa, T_K, grid=grid, reactionEquation=reactionEquation)
    def fn_averageValuesAccordingtoRule(self, prop1, prop2, rule):
        return EOSlist.loaded[self.key].fn_averageValuesAccordingtoRule(prop1, prop2, rule)
