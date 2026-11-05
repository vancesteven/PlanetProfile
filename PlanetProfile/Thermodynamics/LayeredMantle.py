""" Layered silicate mantle from two Perple_X tables.

A mantle made of two compositions, one above and one below an interface (for example the
cumulates left by a partial magma ocean over a pristine lower mantle), is represented by a
single Perple_X table: the rows of the upper table for P < Pinterface_MPa and the rows of the
lower table for P >= Pinterface_MPa, on the identical P-T grid of the two tables. The result is
an ordinary table, so the solver, the EOS-domain checks and the seismic evaluation treat it
like any other mantle EOS:

    from PlanetProfile.Thermodynamics.LayeredMantle import BuildLayeredPerplexTable
    BuildLayeredPerplexTable('upper.tab', 'lower.tab', 2333.8, 'layered.tab')
    Planet.Sil.mantleEOS = 'layered.tab'

The interface is fixed in pressure, not in depth: its depth follows each structure. Its
resolution is one pressure step of the tables. Choosing Pinterface_MPa (for instance the
pressure at a given depth in a reference structure) is left to the caller.
"""
import logging
import os

import numpy as np

from PlanetProfile import _PERPLEXCACHE

log = logging.getLogger('PlanetProfile')


def _ReadTable(fpath, nHeaders):
    with open(fpath) as f:
        lines = f.read().splitlines()
    return lines[:nHeaders], lines[nHeaders:]


def BuildLayeredPerplexTable(upperEOSfname, lowerEOSfname, Pinterface_MPa, outEOSfname,
                             EOSdir=_PERPLEXCACHE, nHeaders=13, overwrite=False):
    """ Write a Perple_X table made of the upper table's rows below Pinterface_MPa and the
        lower table's rows at or above it.

        Parameters
        ----------
        upperEOSfname, lowerEOSfname : str
            Table file names in EOSdir, with identical P-T grids.
        Pinterface_MPa : float
            Interface pressure (MPa). Tables store P in bar, converted exactly (1 bar = 0.1 MPa).
        outEOSfname : str
            Output file name in EOSdir; its header is the lower table's, renamed.
        Returns
        -------
        dict with the output path and the number of rows taken from each table.
    """
    upHdr, upRows = _ReadTable(os.path.join(EOSdir, upperEOSfname), nHeaders)
    loHdr, loRows = _ReadTable(os.path.join(EOSdir, lowerEOSfname), nHeaders)
    if len(upRows) != len(loRows) or upHdr[2:] != loHdr[2:]:
        raise ValueError(f'{upperEOSfname} and {lowerEOSfname} do not share the same P-T grid header.')
    cols = loHdr[-1].split()
    iP = next((i for i, c in enumerate(cols) if c.startswith('P(')), None)
    iT = next((i for i, c in enumerate(cols) if c.startswith('T(')), None)
    if iP is None or iT is None or cols[iP] != 'P(bar)':
        raise ValueError(f'{lowerEOSfname}: expected T(K) and P(bar) columns, found {cols[:2]}.')
    gUp = np.array([[float(r.split()[iT]), float(r.split()[iP])] for r in upRows])
    gLo = np.array([[float(r.split()[iT]), float(r.split()[iP])] for r in loRows])
    if not np.allclose(gUp, gLo, rtol=1e-6, atol=1e-6):
        raise ValueError(f'{upperEOSfname} and {lowerEOSfname}: P-T nodes differ.')
    Pint_bar = Pinterface_MPa * 10.0
    useUpper = gLo[:, 1] < Pint_bar
    if not useUpper.any() or useUpper.all():
        log.warning(f'Interface at {Pinterface_MPa} MPa lies outside the table pressure range; '
                    f'the layered table is a copy of the {"upper" if useUpper.all() else "lower"} table.')
    outPath = os.path.join(EOSdir, outEOSfname)
    if os.path.exists(outPath) and not overwrite:
        raise FileExistsError(f'{outPath} exists; pass overwrite=True to replace it.')
    hdr = list(loHdr)
    hdr[1] = f'{outEOSfname[:-4]}_1.tab'.ljust(len(loHdr[1]))
    body = [u if up else l for u, l, up in zip(upRows, loRows, useUpper)]
    with open(outPath, 'w') as f:
        f.write('\n'.join(hdr + body) + '\n')
    return dict(path=outPath, nUpperRows=int(useUpper.sum()), nLowerRows=int((~useUpper).sum()),
                Pinterface_MPa=float(Pinterface_MPa))
