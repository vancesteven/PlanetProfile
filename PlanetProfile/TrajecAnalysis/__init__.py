"""
Trajectory analysis functions for manipulating and analyzing spacecraft data in comparing to
interior structure models.
"""

import os

import numpy as np

dirName = 'SpacecraftMAGdata'
_MAGdir = os.path.join(os.getcwd(), dirName)
if not os.path.exists(_MAGdir):
    yn = input(f'{dirName} directory not found in pwd. Create it? [y]/n ')
    if yn in ['', 'y', 'Y', 'yes', 'Yes']:
        os.mkdir(_MAGdir)
        _scList = next(os.walk(_MAGdir))[1]
    else:
        print('Aborting.')
        _scList = []

else:
    _scList = next(os.walk(_MAGdir))[1]

_MAGdataList = {
    'Cassini': {
        'Saturn': {revID: f'YYDDD_YYDDD_{revID}_FGM_KRTP_1S.TAB' for revID in
                   ['00', '0A', '0B', '0C'] + [f'{rev:02d}' for rev in range(3, 292)]}
    },
    'Galileo': {
        'Io': {f'{fbID:02d}': f'ORB{fbID:02d}_IO_SYS3.TAB' for fbID in [0, 24, 27, 31, 32]},
        'Europa': {f'{fbID:02d}': f'ORB{fbID:02d}_EUR_SYS3.TAB' for fbID in [4, 11, 12, 14, 15, 19, 26]},
        'Ganymede': {f'{fbID:02d}': f'ORB{fbID:02d}_GAN_SYS3.TAB' for fbID in [1, 2, 7, 8, 28, 29]},
        'Callisto': {f'{fbID:02d}': f'ORB{fbID:02d}_CALL_SYS3.TAB' for fbID in [3, 9, 10, 30]},
        'Jupiter': {f'{fbID:02d}': f'ORB{fbID:02d}_SYS3.TAB' for fbID in list(range(5)) + list(range(6, 35))}
    },
    'Juno': {
        'Europa': {'45': {'272': 'fgm_jno_l3_2022272pc_r1s_v02.sts'}},
        'Ganymede': {'34': {'158': 'fgm_jno_l3_2021158pc_r1s_v01.sts'}},
        'Jupiter': {
            '00': {f'{DDD:03d}': f'fgm_jno_l3_2016{DDD:03d}pc_r1s_v01.sts' for DDD in range(189, 213)},
            '01': {f'{DDD:03d}': f'fgm_jno_l3_2016{DDD:03d}pc_r1s_v01.sts' for DDD in range(214, 267)},
            '02': {**{f'{DDD:03d}': f'fgm_jno_l3_2016{DDD:03d}pc_r1s_v01.sts' for DDD in range(268, 293)},
                   **{f'{DDD:03d}': f'fgm_jno_l3_2016{DDD:03d}pc_r1s_v01.sts' for DDD in range(301, 320)}},
            '03': {**{f'{DDD:03d}': f'fgm_jno_l3_2016{DDD:03d}pc_r1s_v01.sts' for DDD in range(321, 366)},
                   **{f'{DDD:03d}': f'fgm_jno_l3_2017{DDD:03d}pc_r1s_v01.sts' for DDD in range(1, 7)}},
            '04': {f'{DDD:03d}': f'fgm_jno_l3_2017{DDD:03d}pc_r1s_v01.sts' for DDD in range(8, 59)},
            '05': {f'{DDD:03d}': f'fgm_jno_l3_2017{DDD:03d}pc_r1s_v01.sts' for DDD in range(60, 112)},
            '06': {f'{DDD:03d}': f'fgm_jno_l3_2017{DDD:03d}pc_r1s_v01.sts' for DDD in range(113, 165)},
            '07': {f'{DDD:03d}': f'fgm_jno_l3_2017{DDD:03d}pc_r1s_v01.sts' for DDD in range(166, 218)},
            '08': {f'{DDD:03d}': f'fgm_jno_l3_2017{DDD:03d}pc_r1s_v01.sts' for DDD in range(219, 271)},
            '09': {f'{DDD:03d}': f'fgm_jno_l3_2017{DDD:03d}pc_r1s_v01.sts' for DDD in range(272, 324)},
            '10': {**{f'{DDD:03d}': f'fgm_jno_l3_2017{DDD:03d}pc_r1s_v01.sts' for DDD in range(325, 365)},
                   **{f'{DDD:03d}': f'fgm_jno_l3_2018{DDD:03d}pc_r1s_v01.sts' for DDD in range(1, 12)}},
            '11': {f'{DDD:03d}': f'fgm_jno_l3_2018{DDD:03d}pc_r1s_v01.sts' for DDD in range(13, 64)},
            '12': {f'{DDD:03d}': f'fgm_jno_l3_2018{DDD:03d}pc_r1s_v01.sts' for DDD in range(65, 117)},
            '13': {f'{DDD:03d}': f'fgm_jno_l3_2018{DDD:03d}pc_r1s_v01.sts' for DDD in range(118, 170)},
            '14': {f'{DDD:03d}': f'fgm_jno_l3_2018{DDD:03d}pc_r1s_v01.sts' for DDD in range(171, 223)},
            '15': {f'{DDD:03d}': f'fgm_jno_l3_2018{DDD:03d}pc_r1s_v01.sts' for DDD in range(224, 276)},
            '16': {f'{DDD:03d}': f'fgm_jno_l3_2018{DDD:03d}pc_r1s_v01.sts' for DDD in range(277, 329)},
            '17': {**{f'{DDD:03d}': f'fgm_jno_l3_2018{DDD:03d}pc_r1s_v01.sts' for DDD in range(330, 365)},
                   **{f'{DDD:03d}': f'fgm_jno_l3_2019{DDD:03d}pc_r1s_v01.sts' for DDD in range(1, 17)}},
            '18': {f'{DDD:03d}': f'fgm_jno_l3_2019{DDD:03d}pc_r1s_v01.sts' for DDD in range(18, 70)},
            '19': {f'{DDD:03d}': f'fgm_jno_l3_2019{DDD:03d}pc_r1s_v01.sts' for DDD in range(71, 122)},
            '20': {f'{DDD:03d}': f'fgm_jno_l3_2019{DDD:03d}pc_r1s_v01.sts' for DDD in range(123, 175)},
            '21': {f'{DDD:03d}': f'fgm_jno_l3_2019{DDD:03d}pc_r1s_v01.sts' for DDD in range(176, 228)},
            '22': {f'{DDD:03d}': f'fgm_jno_l3_2019{DDD:03d}pc_r1s_v01.sts' for DDD in range(229, 281)},
            '23': {f'{DDD:03d}': f'fgm_jno_l3_2019{DDD:03d}pc_r1s_v01.sts' for DDD in range(282, 334)},
            '24': {**{f'{DDD:03d}': f'fgm_jno_l3_2019{DDD:03d}pc_r1s_v01.sts' for DDD in range(335, 365)},
                   **{f'{DDD:03d}': f'fgm_jno_l3_2020{DDD:03d}pc_r1s_v01.sts' for DDD in range(1, 22)}},
            '25': {f'{DDD:03d}': f'fgm_jno_l3_2020{DDD:03d}pc_r1s_v01.sts' for DDD in range(23, 75)},
            '26': {f'{DDD:03d}': f'fgm_jno_l3_2020{DDD:03d}pc_r1s_v01.sts' for DDD in range(76, 128)},
            '27': {f'{DDD:03d}': f'fgm_jno_l3_2020{DDD:03d}pc_r1s_v01.sts' for DDD in range(130, 180)},
            '28': {f'{DDD:03d}': f'fgm_jno_l3_2020{DDD:03d}pc_r1s_v01.sts' for DDD in range(181, 233)},
            '29': {f'{DDD:03d}': f'fgm_jno_l3_2020{DDD:03d}pc_r1s_v01.sts' for DDD in range(234, 286)},
            '30': {f'{DDD:03d}': f'fgm_jno_l3_2020{DDD:03d}pc_r1s_v01.sts' for DDD in range(287, 339)},
            '31': {**{f'{DDD:03d}': f'fgm_jno_l3_2020{DDD:03d}pc_r1s_v01.sts' for DDD in range(340, 366)},
                   **{f'{DDD:03d}': f'fgm_jno_l3_2021{DDD:03d}pc_r1s_v01.sts' for DDD in range(1, 26)}},
            '32': {f'{DDD:03d}': f'fgm_jno_l3_2021{DDD:03d}pc_r1s_v01.sts' for DDD in range(27, 79)},
            '33': {f'{DDD:03d}': f'fgm_jno_l3_2021{DDD:03d}pc_r1s_v01.sts' for DDD in range(80, 132)},
            '34': {f'{DDD:03d}': f'fgm_jno_l3_2021{DDD:03d}pc_r1s_v01.sts' for DDD in range(133, 180)},
            '35': {f'{DDD:03d}': f'fgm_jno_l3_2021{DDD:03d}pc_r1s_v01.sts' for DDD in range(181, 224)},
            '36': {f'{DDD:03d}': f'fgm_jno_l3_2021{DDD:03d}pc_r1s_v01.sts' for DDD in range(225, 267)},
            '37': {f'{DDD:03d}': f'fgm_jno_l3_2021{DDD:03d}pc_r1s_v01.sts' for DDD in range(268, 311)},
            '38': {f'{DDD:03d}': f'fgm_jno_l3_2021{DDD:03d}pc_r1s_v01.sts' for DDD in range(312, 355)},
            '39': {**{f'{DDD:03d}': f'fgm_jno_l3_2021{DDD:03d}pc_r1s_v01.sts' for DDD in range(356, 365)},
                   **{f'{DDD:03d}': f'fgm_jno_l3_2022{DDD:03d}pc_r1s_v01.sts' for DDD in range(1, 34)}},
            '40': {f'{DDD:03d}': f'fgm_jno_l3_2022{DDD:03d}pc_r1s_v01.sts' for DDD in range(35, 77)},
            '41': {f'{DDD:03d}': f'fgm_jno_l3_2022{DDD:03d}pc_r1s_v01.sts' for DDD in range(78, 121)},
            '42': {f'{DDD:03d}': f'fgm_jno_l3_2022{DDD:03d}pc_r1s_v01.sts' for DDD in range(122, 164)},
            '43': {f'{DDD:03d}': f'fgm_jno_l3_2022{DDD:03d}pc_r1s_v01.sts' for DDD in range(165, 208)},
            '44': {f'{DDD:03d}': f'fgm_jno_l3_2022{DDD:03d}pc_r1s_v01.sts' for DDD in range(209, 251)},
            '45': {f'{DDD:03d}': f'fgm_jno_l3_2022{DDD:03d}pc_r1s_v01.sts' for DDD in range(252, 291)},
            '46': {f'{DDD:03d}': f'fgm_jno_l3_2022{DDD:03d}pc_r1s_v01.sts' for DDD in range(292, 330)},
            '47': {**{f'{DDD:03d}': f'fgm_jno_l3_2022{DDD:03d}pc_r1s_v01.sts' for DDD in range(331, 348)},
                      f'003': f'fgm_jno_l3_2023003pc_r1s_v01.sts'},
            '48': {f'{DDD:03d}': f'fgm_jno_l3_2023{DDD:03d}pc_r1s_v01.sts' for DDD in range(4, 41)},
            '49': {f'{DDD:03d}': f'fgm_jno_l3_2023{DDD:03d}pc_r1s_v01.sts' for DDD in range(42, 79)},
            '50': {f'{DDD:03d}': f'fgm_jno_l3_2023{DDD:03d}pc_r1s_v01.sts' for DDD in range(80, 117)}
        }
    }
}

# Correct v02 items in Juno data
_MAGdataList['Juno']['Jupiter']['01']['240'] = 'fgm_jno_l3_2016240pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['03']['346'] = 'fgm_jno_l3_2016346pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['04']['033'] = 'fgm_jno_l3_2017033pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['05']['086'] = 'fgm_jno_l3_2017086pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['06']['139'] = 'fgm_jno_l3_2017139pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['07']['191'] = 'fgm_jno_l3_2017191pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['07']['192'] = 'fgm_jno_l3_2017192pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['08']['244'] = 'fgm_jno_l3_2017244pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['08']['245'] = 'fgm_jno_l3_2017245pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['09']['297'] = 'fgm_jno_l3_2017297pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['10']['350'] = 'fgm_jno_l3_2017350pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['10']['351'] = 'fgm_jno_l3_2017351pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['11']['038'] = 'fgm_jno_l3_2018038pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['12']['091'] = 'fgm_jno_l3_2018091pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['13']['144'] = 'fgm_jno_l3_2018144pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['14']['197'] = 'fgm_jno_l3_2018197pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['15']['249'] = 'fgm_jno_l3_2018249pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['15']['250'] = 'fgm_jno_l3_2018250pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['19']['096'] = 'fgm_jno_l3_2019096pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['20']['149'] = 'fgm_jno_l3_2019149pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['21']['201'] = 'fgm_jno_l3_2019201pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['21']['202'] = 'fgm_jno_l3_2019202pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['22']['254'] = 'fgm_jno_l3_2019254pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['22']['255'] = 'fgm_jno_l3_2019255pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['23']['307'] = 'fgm_jno_l3_2019307pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['23']['308'] = 'fgm_jno_l3_2019308pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['24']['360'] = 'fgm_jno_l3_2019360pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['25']['048'] = 'fgm_jno_l3_2020048pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['26']['101'] = 'fgm_jno_l3_2020101pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['27']['154'] = 'fgm_jno_l3_2020154pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['28']['207'] = 'fgm_jno_l3_2020207pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['29']['259'] = 'fgm_jno_l3_2020259pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['29']['260'] = 'fgm_jno_l3_2020260pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['30']['312'] = 'fgm_jno_l3_2020312pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['30']['313'] = 'fgm_jno_l3_2020313pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['31']['365'] = 'fgm_jno_l3_2020365pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['31']['366'] = 'fgm_jno_l3_2020366pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['32']['052'] = 'fgm_jno_l3_2021052pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['33']['105'] = 'fgm_jno_l3_2021105pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['33']['106'] = 'fgm_jno_l3_2021106pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['34']['159'] = 'fgm_jno_l3_2021159pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['35']['202'] = 'fgm_jno_l3_2021202pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['36']['245'] = 'fgm_jno_l3_2021245pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['36']['246'] = 'fgm_jno_l3_2021246pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['37']['289'] = 'fgm_jno_l3_2021289pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['38']['333'] = 'fgm_jno_l3_2021333pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['39']['012'] = 'fgm_jno_l3_2022012pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['40']['055'] = 'fgm_jno_l3_2022055pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['40']['056'] = 'fgm_jno_l3_2022056pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['41']['099'] = 'fgm_jno_l3_2022099pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['42']['142'] = 'fgm_jno_l3_2022142pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['42']['143'] = 'fgm_jno_l3_2022143pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['43']['186'] = 'fgm_jno_l3_2022186pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['44']['229'] = 'fgm_jno_l3_2022229pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['45']['272'] = 'fgm_jno_l3_2022272pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['46']['310'] = 'fgm_jno_l3_2022310pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['46']['311'] = 'fgm_jno_l3_2022311pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['47']['348'] = 'fgm_jno_l3_2022348pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['48']['022'] = 'fgm_jno_l3_2023022pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['49']['060'] = 'fgm_jno_l3_2023060pc_r1s_v02.sts'
_MAGdataList['Juno']['Jupiter']['50']['098'] = 'fgm_jno_l3_2023098pc_r1s_v02.sts'
