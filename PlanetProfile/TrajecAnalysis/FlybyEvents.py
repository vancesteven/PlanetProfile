""" Lists of significant events and related information needed to load
    trajectories from SPICE kernels and measurements from MAG data files.
"""

import logging
import numpy as np
import spiceypy as spice
from PlanetProfile.GetConfig import Params
from PlanetProfile.TrajecAnalysis.SpiceFuncs import LoadKernels, BodyDist_km, BiTrajec, spiceSCname

# Assign logger
log = logging.getLogger('PlanetProfile')

scPlanets = {
    'Galileo': 'Jupiter',
    'Cassini': 'Saturn',
    'Juno': 'Jupiter'
}
scNames = list(scPlanets.keys())
[LoadKernels(Params, parent, scName) for scName, parent in scPlanets.items()]

class FlybyCAStruct:
    def __init__(self, scName):
        self.scName = scName
        self.spiceName = spiceSCname[self.scName]

        self.rCA_km = None  # Dict of body: flyby ID: closest approach radii in km
        # Dict of body: flyby ID: closest approach UTC timestamp in SPICE-readable format
        # NOTE: These times are evaluated from the archived SPICE kernels for trajectory
        # reconstruction as of 2022-07-02.
        if self.scName == 'Galileo':
            self.tCA_UTC = {
                'Io': {
                    'I0': '1995-12-07T17:45:58.461',
                    'I24': '1999-10-11T04:33:02.580',
                    'I27': '2000-02-22T13:46:41.423',
                    'I31': '2001-08-06T04:59:20.516',
                    'I32': '2001-10-16T01:23:20.598'
                },
                'Europa': {
                    'E4': '1996-12-19T06:52:57.758',
                    'E11': '1997-11-06T20:31:44.224',
                    'E12': '1997-12-16T12:03:19.869',
                    'E14': '1998-03-29T13:21:05.161',
                    'E15': '1998-05-31T21:12:56.608',
                    'E19': '1999-02-01T02:19:50.002',
                    'E26': '2000-01-03T17:59:42.595'
                },
                'Ganymede': {
                    'G1': '1996-06-27T06:29:06.687',
                    'G2': '1996-09-06T18:59:33.836',
                    'G7': '1997-04-05T07:09:58.113',
                    'G8': '1997-05-07T15:56:09.556',
                    'G28': '2000-05-20T10:10:09.662',
                    'G29': '2000-12-28T08:25:26.659'
                },
                'Callisto': {
                    'C3': '1996-11-04T13:34:27.726',
                    'C9': '1997-06-25T13:47:49.949',
                    'C10': '1997-09-17T00:18:54.790',
                    'C30': '2001-05-25T11:23:57.770'
                }
            }
        elif self.scName == 'Cassini':
            self.tCA_UTC = {
                'Enceladus': {
                    'E0': '2005-02-17T03:30:28.745',
                    'E1': '2005-03-09T09:08:02.313',
                    'E2': '2005-07-14T19:55:20.924',
                    'E3': '2008-03-12T19:06:11.827',
                    'E4': '2008-08-11T21:06:18.646',
                    'E5': '2008-10-09T19:06:39.789',
                    'E6': '2008-10-31T17:14:51.514',
                    'E7': '2009-11-02T07:41:57.693',
                    'E8': '2009-11-21T02:09:56.440',
                    'E9': '2010-04-28T00:10:16.852',
                    'E10': '2010-05-18T06:04:39.708',
                    'E11': '2010-08-13T22:30:51.619',
                    'E12': '2010-11-30T11:53:59.407',
                    'E13': '2010-12-21T01:08:26.768',
                    'E14': '2011-10-01T13:52:25.736',
                    'E15': '2011-10-19T09:22:11.790',
                    'E16': '2011-11-06T04:58:52.984',
                    'E17': '2012-03-27T18:30:09.015',
                    'E18': '2012-04-14T14:01:37.801',
                    'E19': '2012-05-02T09:31:28.806',
                    'E20': '2015-10-14T10:41:29.053',
                    'E21': '2015-10-28T15:22:41.596',
                    'E22': '2015-12-19T17:49:16.159'
                }
            }
        elif self.scName == 'Juno':
            self.tCA_UTC = {
                'Ganymede': {
                    'G34': '2021-06-07T16:56:08.733'
                }
            }

        self.etCA = {
            body: {flybyID: spice.str2et(tCA) for flybyID, tCA in self.tCA_UTC[body].items()}
        for body in self.tCA_UTC.keys()}

    def GetrCA(self):
        self.rCA_km = {body: {flybyID: rCA for flybyID, rCA in zip(self.tCA_UTC[body].keys(),
                                BodyDist_km(self.spiceName, body, list(self.etCA[body].values()))[-1])}
                                for body in self.tCA_UTC.keys()}

    def PrinttCA(self, flybyID=None):
        # Get the actual time of closest approach according to trajectory
        # reconstruction and print.
        print(f' - {self.scName.upper()} - ')
        if flybyID is None:
            for body in self.tCA_UTC.keys():
                print(f'{body}:')
                for ID, tCA_UTC in self.tCA_UTC[body].items():
                    tCA = GetActualCA(self.spiceName, tCA_UTC, body)
                    print(f'{ID}: {tCA}')
        else:
            body = next(bname for bname in self.tCA_UTC.keys() if bname[0] == flybyID[0])
            tCA = GetActualCA(self.spiceName, self.tCA_UTC[body][flybyID], body)
            print(f'{flybyID}: {tCA}')
                    

def GetActualCA(spiceSCname, t_UTC, bodyname, range_min=5, res_s=0.001):
    # Look within +/- range_min of t_UTC to within res_s precision
    # to get the UTC time string of the closest approach.

    spiceBody = bodyname.upper()
    etApprox = spice.str2et(t_UTC)
    ets = np.arange(etApprox - range_min * 60, etApprox + range_min * 60, res_s)
    pos, _ = spice.spkpos(spiceSCname, ets, f'IAU_{spiceBody}', 'NONE', spiceBody)
    r_km = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
    etCA = ets[np.argmin(r_km)]
    tCA_UTC = spice.et2utc(etCA, 'ISOC', 3)

    return tCA_UTC


FlybyCA = {scName: FlybyCAStruct(scName) for scName in scNames}
[FlybyCA[sc].GetrCA() for sc in scNames]
