""" Lists of significant events and related information needed to load
    trajectories from SPICE kernels and measurements from MAG data files.
"""

import logging
import numpy as np
import spiceypy as spice
from PlanetProfile.GetConfig import Params
from PlanetProfile.TrajecAnalysis.SpiceFuncs import LoadKernels, BodyDist_km, BodyVel_kms, \
    spiceSCname

# Assign logger
log = logging.getLogger('PlanetProfile')

scPlanets = {
    'Cassini': 'Saturn',
    'Clipper': 'Jupiter',
    'Galileo': 'Jupiter',
    'JUICE': 'Jupiter',
    'Juno': 'Jupiter',
    'Voyager 1': ['Jupiter', 'Saturn'],
    'Voyager 2': ['Jupiter', 'Saturn', 'Uranus', 'Neptune']
}
scNames = list(scPlanets.keys())

scTargets = {
    'Cassini': ['Mimas', 'Enceladus', 'Tethys', 'Rhea', 'Titan', 'Iapetus'],
    'Clipper': ['Europa', 'Ganymede', 'Callisto'],
    'Galileo': ['Io', 'Europa', 'Ganymede', 'Callisto'],
    'JUICE': ['Europa', 'Ganymede', 'Callisto'],
    'Juno': ['Io', 'Europa', 'Ganymede'],
    'Voyager 1': ['Jupiter', 'Saturn'],
    'Voyager 2': ['Jupiter', 'Saturn', 'Uranus', 'Neptune']
}


class FlybyCAStruct:
    def __init__(self, scName):
        [LoadKernels(Params, parent, scName) for scName, parent in scPlanets.items()]
        self.scName = scName
        self.spiceName = spiceSCname[self.scName]

        self.rCA_km = {}  # Dict of body: flyby ID: closest approach radii in km
        self.vCA_kms = {}  # Dict of body: flyby ID: closest approach relative speed in km/s
        self.latCA_deg = {}  # Dict of body: flyby ID: closest approach IAU latitude in degrees
        self.lonCA_deg = {}  # Dict of body: flyby ID: closest approach IAU longitude in degrees
        self.latSunCA_deg = {}  # Dict of body: flyby ID: sub-solar point IAU latitude in degrees
        self.lonSunCA_deg = {}  # Dict of body: flyby ID: sub-solar point IAU longitude in degrees
        # Planetary Longitude System (ULS and NLS for Uranus and Neptune respectively) coordinates
        self.rPLS_km = {}  # Dict of body: flyby ID: PLS radii in km
        self.latPLS_deg = {}  # Dict of body: flyby ID: PLS latitude in degrees
        self.lonPLS_deg = {}  # Dict of body: flyby ID: PLS longitude in degrees
        # Dict of body: flyby ID: closest approach UTC timestamp in SPICE-readable format
        # NOTE: These times are evaluated from the archived SPICE kernels for trajectory
        # reconstruction as of 2022-07-02.
        if self.scName == 'Cassini':
            self.tCA_UTC = {
                'Enceladus': {
                    '00': '2005-02-17T03:30:28.745',
                    '01': '2005-03-09T09:08:02.313',
                    '02': '2005-07-14T19:55:20.924',
                    '03': '2008-03-12T19:06:11.827',
                    '04': '2008-08-11T21:06:18.646',
                    '05': '2008-10-09T19:06:39.789',
                    '06': '2008-10-31T17:14:51.514',
                    '07': '2009-11-02T07:41:57.693',
                    '08': '2009-11-21T02:09:56.440',
                    '09': '2010-04-28T00:10:16.852',
                    '10': '2010-05-18T06:04:39.708',
                    '11': '2010-08-13T22:30:51.619',
                    '12': '2010-11-30T11:53:59.407',
                    '13': '2010-12-21T01:08:26.768',
                    '14': '2011-10-01T13:52:25.736',
                    '15': '2011-10-19T09:22:11.790',
                    '16': '2011-11-06T04:58:52.984',
                    '17': '2012-03-27T18:30:09.015',
                    '18': '2012-04-14T14:01:37.801',
                    '19': '2012-05-02T09:31:28.806',
                    '20': '2015-10-14T10:41:29.053',
                    '21': '2015-10-28T15:22:41.596',
                    '22': '2015-12-19T17:49:16.159'
                }
            }
        elif self.scName == 'Clipper':
            self.tCA_UTC = {
                'Europa': {
                    '01': '2030-06-07T16:56:08.733'
                }
            }
        elif self.scName == 'Galileo':
            self.tCA_UTC = {
                'Io': {
                    '00': '1995-12-07T17:45:58.461',
                    '24': '1999-10-11T04:33:02.580',
                    '27': '2000-02-22T13:46:41.423',
                    '31': '2001-08-06T04:59:20.516',
                    '32': '2001-10-16T01:23:20.598'
                },
                'Europa': {
                    '04': '1996-12-19T06:52:57.758',
                    '11': '1997-11-06T20:31:44.224',
                    '12': '1997-12-16T12:03:19.869',
                    '14': '1998-03-29T13:21:05.161',
                    '15': '1998-05-31T21:12:56.608',
                    '19': '1999-02-01T02:19:50.002',
                    '26': '2000-01-03T17:59:42.595'
                },
                'Ganymede': {
                    '01': '1996-06-27T06:29:06.687',
                    '02': '1996-09-06T18:59:33.836',
                    '07': '1997-04-05T07:09:58.113',
                    '08': '1997-05-07T15:56:09.556',
                    '28': '2000-05-20T10:10:09.662',
                    '29': '2000-12-28T08:25:26.659'
                },
                'Callisto': {
                    '03': '1996-11-04T13:34:27.726',
                    '09': '1997-06-25T13:47:49.949',
                    '10': '1997-09-17T00:18:54.790',
                    '20': '1999-05-05T13:56:18.118',
                    '21': '1999-06-30T07:46:49.675',
                    '22': '1999-08-14T08:30:51.767',
                    '23': '1999-09-16T17:27:01.828',
                    '30': '2001-05-25T11:23:57.770'
                }
            }
        elif self.scName == 'JUICE':
            self.tCA_UTC = {
            }
        elif self.scName == 'Juno':
            self.tCA_UTC = {
                'Io':       {
                  '56': '2023-10-15T06:47:15.907',
                  '57': '2023-12-30T08:37:44.609',
                  '58': '2024-02-03T17:48:35.293'
                },
                'Europa':   {'45': '2022-09-29T09:36:28.827'},
                'Ganymede': {'34': '2021-06-07T16:56:08.733'}
            }
        elif self.scName == 'Voyager 1':
            self.tCA_UTC = {
                'Jupiter': {'J': '1979-03-05T12:04:35.405'},
                'Saturn':  {'S': '1980-11-12T23:45:42.725'}
            }
        elif self.scName == 'Voyager 2':
            self.tCA_UTC = {
                'Jupiter': {'J': '1979-07-09T22:29:01.964'},
                'Saturn':  {'S': '1981-08-26T03:24:04.760'},
                'Uranus':  {'U': '1986-01-24T17:58:51.346'},
                'Neptune': {'N': '1989-08-25T03:55:40.076'}
            }
            self.tPLS_UTC = {
                'Uranus':  {'U': '1986-01-24T18:00:00.000'},
                'Neptune': {'N': '1989-08-25T03:56:00.000'}
            }

        self.etCA = {
            body: {flybyID: spice.str2et(tCA) for flybyID, tCA in self.tCA_UTC[body].items() if tCA is not None}
        for body in self.tCA_UTC.keys()}

        if self.scName == 'Voyager 2':
            # Dict of body: flyby ID: specific times for definitions of planetary longitude system (ULS, NLS) frames
            self.etPLS = {
                body: {flybyID: spice.str2et(tPLS) for flybyID, tPLS in self.tPLS_UTC[body].items() if tPLS is not None}
            for body in self.tPLS_UTC.keys()}
            self.GetrPLS()
        else:
            self.etPLS = None

        self.GetCAvals()
        self.GetSunDir()

    def GetCAvals(self):
        for body in self.tCA_UTC.keys():
            x_km, y_km, z_km, r_km = BodyDist_km(self.spiceName, body, list(self.etCA[body].values()))
            _, _, _, v_kms = BodyVel_kms(self.spiceName, body, list(self.etCA[body].values()))
            lat_deg = np.degrees(np.arctan2(z_km, np.sqrt(x_km**2 + y_km**2)))
            lon_deg = np.degrees(np.arctan2(y_km, x_km))
            self.rCA_km[body] = {flybyID: rCA for flybyID, rCA in zip(self.tCA_UTC[body].keys(), r_km)}
            self.vCA_kms[body] = {flybyID: vCA for flybyID, vCA in zip(self.tCA_UTC[body].keys(), v_kms)}
            self.latCA_deg[body] = {flybyID: latCA for flybyID, latCA in zip(self.tCA_UTC[body].keys(), lat_deg)}
            self.lonCA_deg[body] = {flybyID: lonCA for flybyID, lonCA in zip(self.tCA_UTC[body].keys(), lon_deg)}

    def GetrPLS(self):
        if self.etPLS is not None:
            for body in self.tPLS_UTC.keys():
                x_km, y_km, z_km, r_km = BodyDist_km(self.spiceName, body, list(self.etPLS[body].values()))
                lat_deg = np.degrees(np.arctan2(z_km, np.sqrt(x_km**2 + y_km**2)))
                lon_deg = np.degrees(np.arctan2(y_km, x_km))
                self.rPLS_km[body] = {flybyID: rPLS for flybyID, rPLS in zip(self.tPLS_UTC[body].keys(), r_km)}
                self.latPLS_deg[body] = {flybyID: latPLS for flybyID, latPLS in zip(self.tPLS_UTC[body].keys(), lat_deg)}
                self.lonPLS_deg[body] = {flybyID: lonPLS for flybyID, lonPLS in zip(self.tPLS_UTC[body].keys(), lon_deg)}

    def GetSunDir(self):
        for body in self.tCA_UTC.keys():
            x_km, y_km, z_km, _ = BodyDist_km('SUN', body, list(self.etCA[body].values()))
            lat_deg = np.degrees(np.arctan2(z_km, np.sqrt(x_km**2 + y_km**2)))
            lon_deg = np.degrees(np.arctan2(y_km, x_km))
            self.latSunCA_deg[body] = {flybyID: latCA for flybyID, latCA in zip(self.tCA_UTC[body].keys(), lat_deg)}
            self.lonSunCA_deg[body] = {flybyID: lonCA for flybyID, lonCA in zip(self.tCA_UTC[body].keys(), lon_deg)}

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


def GetFlybyCA():
    return {scName: FlybyCAStruct(scName) for scName in scNames}

if __name__ == '__main__':
    FlybyCA = GetFlybyCA()
    a=0
