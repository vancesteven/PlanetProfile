"""Build a traceable Fe-S core EOS table from the werami tables in this repository.

Why: the distributed PlanetProfile/Thermodynamics/EOStables/Perple_X/Fe-S_3D_EOS.mat was
built by stackFeSPerplexTables.py, which converts bar with the atmosphere factor
(bar2MPa = 1.01325e-1, bar2GPa = 1.01325e-4). Its pressure axis and its K_S and G_S are
therefore scaled by +1.325 %, its Fe95-S5 slice has moduli 1e5-1e6 times too small, and the
werami sources of its eight compositions are not in the repository.

This script stacks the five werami tables that ARE in the repository, without any P-T
resampling, using only exact unit conversions:
    P:      bar * 0.1    -> MPa
    KS, GS: bar * 1e-4   -> GPa
    Cp:     J/K/m^3 / rho -> J/kg/K
Every output cell reproduces its source cell after conversion, and the result is checked
before it is written. The product is identified by a semantic SHA-256 over its arrays, which
does not depend on MAT-file timestamps.

Output (default): Fe-S_3D_EOS_5slice_traceable_v1.mat in the PlanetProfile Perple_X cache,
the location PerplexEOSStruct reads. Select it with
    Planet.Core.coreEOS = 'Fe-S_3D_EOS_5slice_traceable_v1.mat'
    Planet.Core.wFe_ppt = 950   # one of 800, 850, 900, 950, 1000
The five slices are used as they are; compositions between slices are interpolated by the
loader, which has not been validated scientifically for this table.

Usage:
    python Thermodynamics/Perple_X/stackFeSTraceable.py [output.mat]

Same construction as the product used by the lunar interior study
(lunar-bsm-interior-models: scripts/build_fe_s_5slice_core_eos.py, docs/FE_S_CORE_EOS_PRODUCT.md),
whose semantic SHA-256 is SEMANTIC_SHA256 below.
"""
import hashlib
import os
import sys

import numpy as np
from scipy.io import savemat

HERE = os.path.dirname(os.path.abspath(__file__))
WERAMI_DIR = os.path.join(HERE, 'werami_tables')
SOURCE_FILES = ['80Fe20S_core.tab', '85Fe15S_core.tab', '90Fe10S_core.tab', '95Fe5S_core.tab',
                '100Fe0S_core.tab']
WFE_PPT = np.array([800.0, 850.0, 900.0, 950.0, 1000.0])
OUTPUT_NAME = 'Fe-S_3D_EOS_5slice_traceable_v1.mat'
SEMANTIC_SHA256 = '452147a82062916219dc8b83e53dcb53114ccfe8f60030b8afb4a9db72834445'

BAR2MPA = 1e-1  # exact
BAR2GPA = 1e-4  # exact
FIELDS = ['wFe_ppt', 'wS_ppt', 'P_MPa', 'T_K', 'rho_kgm3', 'VP_kms', 'VS_kms', 'Cp_JkgK',
          'alpha_pK', 'KS_GPa', 'GS_GPa']
PROPS = FIELDS[4:]


def ParseWeramiTable(fpath, nHeaders=13):
    """ Parse a native werami table with T as the first (fastest) axis and P as the second.
        Returns the axes and the property arrays with shape (nP, nT), in the table's units:
        rows follow P and columns follow T, the layout of the 3D Fe-S tables read by
        PerplexEOSStruct.
    """
    with open(fpath) as f:
        lines = [f.readline() for _ in range(nHeaders)]
    if not (lines[3].strip().startswith('T') and lines[7].strip().startswith('P')):
        raise ValueError(f'{fpath}: expected T as the first axis and P as the second.')
    Tmin, dT, nT = float(lines[4]), float(lines[5]), int(lines[6])
    Pmin, dP, nP = float(lines[8]), float(lines[9]), int(lines[10])
    data = np.loadtxt(fpath, skiprows=nHeaders)
    if data.shape[0] != nT * nP:
        raise ValueError(f'{fpath}: expected {nT * nP} rows, found {data.shape[0]}.')
    cols = [data[:, i].reshape(nP, nT) for i in range(9)]
    return dict(T_K=Tmin + dT * np.arange(nT), P_bar=Pmin + dP * np.arange(nP), T_col=cols[0],
                P_col=cols[1], rho_kgm3=cols[2], VP_kms=cols[3], VS_kms=cols[4], Cp_Jm3K=cols[5],
                alpha_pK=cols[6], KS_bar=cols[7], GS_bar=cols[8])


def StackTraceableFeS(weramiDir=WERAMI_DIR):
    """ Stack the five werami tables into the arrays of a 3D Fe-S EOS (slice, P, T). """
    srcs = [ParseWeramiTable(os.path.join(weramiDir, f)) for f in SOURCE_FILES]
    T_K, P_bar = srcs[0]['T_K'], srcs[0]['P_bar']
    for f, s in zip(SOURCE_FILES, srcs):
        if not (np.array_equal(s['T_K'], T_K) and np.array_equal(s['P_bar'], P_bar)):
            raise ValueError(f'{f}: P-T axes differ from {SOURCE_FILES[0]}.')
        if not (np.allclose(s['T_col'], T_K[None, :]) and np.allclose(s['P_col'], P_bar[:, None])):
            raise ValueError(f'{f}: data columns do not follow the header axes.')
    out = dict(wFe_ppt=WFE_PPT.copy(), wS_ppt=1e3 - WFE_PPT, P_MPa=P_bar * BAR2MPA, T_K=T_K.copy())
    out['rho_kgm3'] = np.stack([s['rho_kgm3'] for s in srcs])
    out['VP_kms'] = np.stack([s['VP_kms'] for s in srcs])
    out['VS_kms'] = np.stack([s['VS_kms'] for s in srcs])
    out['Cp_JkgK'] = np.stack([s['Cp_Jm3K'] / s['rho_kgm3'] for s in srcs])
    out['alpha_pK'] = np.stack([s['alpha_pK'] for s in srcs])
    out['KS_GPa'] = np.stack([s['KS_bar'] * BAR2GPA for s in srcs])
    out['GS_GPa'] = np.stack([s['GS_bar'] * BAR2GPA for s in srcs])
    CheckTraceableFeS(out, srcs)
    return out


def CheckTraceableFeS(out, srcs):
    """ Raise if the stacked arrays are not a faithful, physical copy of their sources. """
    for k in PROPS:
        if not np.all(np.isfinite(out[k])):
            raise ValueError(f'{k}: non-finite values.')
    for k in ['rho_kgm3', 'VP_kms', 'VS_kms', 'KS_GPa', 'GS_GPa']:
        if not np.all(out[k] > 0):
            raise ValueError(f'{k}: {np.sum(out[k] <= 0)} non-positive values.')
    for i, s in enumerate(srcs):
        for k, src in [('rho_kgm3', s['rho_kgm3']), ('VP_kms', s['VP_kms']), ('VS_kms', s['VS_kms']),
                       ('KS_GPa', s['KS_bar'] * BAR2GPA), ('GS_GPa', s['GS_bar'] * BAR2GPA)]:
            if not np.array_equal(out[k][i], src):
                raise ValueError(f'{k}, slice {SOURCE_FILES[i]}: not an exact copy of its source.')
    # Elastic identity VS = sqrt(GS / rho), at the precision of the tabulated values
    VS_from_G = np.sqrt(out['GS_GPa'] * 1e9 / out['rho_kgm3']) / 1e3
    if np.median(np.abs(VS_from_G - out['VS_kms'])) * 1e3 > 0.1:
        raise ValueError('VS does not match sqrt(GS/rho) at table precision.')


def SemanticSHA256(out):
    """ SHA-256 over field name, dtype, shape and C-contiguous little-endian float64 bytes. """
    h = hashlib.sha256()
    for k in FIELDS:
        a = np.ascontiguousarray(out[k], dtype='<f8')
        h.update(k.encode()); h.update(a.dtype.str.encode()); h.update(str(a.shape).encode())
        h.update(a.tobytes())
    return h.hexdigest()


def WriteTraceableFeS(out, fpath):
    """ Write the MAT file in the schema PerplexEOSStruct reads for 3D Fe-S tables. """
    d = {k: out[k].reshape(1, -1) for k in FIELDS[:4]}
    d.update({k: out[k] for k in PROPS})
    d['Perplex version'] = np.array(['6.6.6'], dtype='U5')
    savemat(fpath, d, format='5', do_compression=False)


def main(argv):
    if len(argv) > 1:
        fpath = argv[1]
    else:
        from PlanetProfile import _PERPLEXCACHE
        fpath = os.path.join(_PERPLEXCACHE, OUTPUT_NAME)
    out = StackTraceableFeS()
    sha = SemanticSHA256(out)
    if sha != SEMANTIC_SHA256:
        raise ValueError(f'Semantic SHA-256 {sha} differs from the reference product {SEMANTIC_SHA256}.')
    WriteTraceableFeS(out, fpath)
    print(f'Wrote {fpath}\nsemantic SHA-256 {sha} (matches the reference product)')


if __name__ == '__main__':
    main(sys.argv)
