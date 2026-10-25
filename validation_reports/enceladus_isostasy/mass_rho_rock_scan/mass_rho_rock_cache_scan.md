# mass_rho_rock cache scan (no rebuilds)

cache: /Users/svance/scratch/enceladus_zbw_production_cache.pkl

## rho_ice = none, rho_s = 2200

ocean nodes 3399 (None 121, frozen/no-ocean skipped 0), feasible 3399, rejects 0 {}
- R_o - R_s km: {'min': 2.5904041182057407, 'median': 20.454767458628776, 'max': 49.63608377763088}
- derived C/MR2: {'min': 0.3424351245259057, 'median': 0.3427957398162167, 'max': 0.3469320456982632}
- R_s - R_sil(PP) km: {'min': 5.909233930252114, 'median': 7.519465964866715, 'max': 10.9552793412678}
- d_rho_o (mean[R_s,R_o]-mean[R_sil,R_o]) kg/m3: {'min': -0.26635464575269907, 'max': -0.16623371218349803}
- double-count shifts (sigma_obs; negligible < 0.05): {'C20': {'max_abs_shift_sigma_obs': 0.006946018789911278, 'at_node': [0, 39], 'negligible_lt_0.05': True}, 'C22': {'max_abs_shift_sigma_obs': 0.0051817059421708864, 'at_node': [0, 39], 'negligible_lt_0.05': True}, 'libration_deg': {'max_abs_shift_sigma_obs': 0.515984375069378, 'at_node': [0, 39], 'negligible_lt_0.05': False}}

binding nodes, thinnest ocean:
- (87,9) zb=45.00 w=0.492 R_o-R_s=2.590 km C=0.34244
- (87,10) zb=45.00 w=0.588 R_o-R_s=2.591 km C=0.34244
- (87,11) zb=45.00 w=0.702 R_o-R_s=2.591 km C=0.34244
- (87,12) zb=45.00 w=0.838 R_o-R_s=2.592 km C=0.34244
- (87,13) zb=45.00 w=1.000 R_o-R_s=2.592 km C=0.34244
- (87,0) zb=45.00 w=0.100 R_o-R_s=2.592 km C=0.34244
- (87,1) zb=45.00 w=0.119 R_o-R_s=2.592 km C=0.34244
- (87,2) zb=45.00 w=0.143 R_o-R_s=2.592 km C=0.34244
- (87,3) zb=45.00 w=0.170 R_o-R_s=2.592 km C=0.34244
- (87,4) zb=45.00 w=0.203 R_o-R_s=2.593 km C=0.34244

binding nodes, max |dR_s|:
- (14,39) zb=19.00 w=100.000 dR_s=+10.955 km
- (9,39) zb=14.00 w=100.000 dR_s=+10.878 km
- (17,39) zb=21.00 w=100.000 dR_s=+10.868 km
- (2,39) zb=7.00 w=100.000 dR_s=+10.854 km
- (4,39) zb=9.00 w=100.000 dR_s=+10.746 km
- (11,39) zb=16.00 w=100.000 dR_s=+10.731 km
- (20,39) zb=22.25 w=100.000 dR_s=+10.726 km
- (13,38) zb=18.00 w=83.768 dR_s=+10.689 km
- (18,39) zb=21.50 w=100.000 dR_s=+10.665 km
- (8,38) zb=13.00 w=83.768 dR_s=+10.663 km

Rejects:

## rho_ice = none, rho_s = 2600

ocean nodes 3399 (None 121, frozen/no-ocean skipped 0), feasible 3399, rejects 0 {}
- R_o - R_s km: {'min': 21.327339082038787, 'median': 39.14823933086972, 'max': 68.75092060967157}
- derived C/MR2: {'min': 0.324158090459218, 'median': 0.3250880224953872, 'max': 0.3318172653358165}
- R_s - R_sil(PP) km: {'min': -12.824926500385597, 'median': -11.249156771318376, 'max': -8.302158964541276}
- d_rho_o (mean[R_s,R_o]-mean[R_sil,R_o]) kg/m3: {'min': 0.1415104951174726, 'max': 0.2309800567967386}
- double-count shifts (sigma_obs; negligible < 0.05): {'C20': {'max_abs_shift_sigma_obs': 0.008621890351784569, 'at_node': [0, 16], 'negligible_lt_0.05': True}, 'C22': {'max_abs_shift_sigma_obs': 0.006021073686525515, 'at_node': [0, 39], 'negligible_lt_0.05': True}, 'libration_deg': {'max_abs_shift_sigma_obs': 0.2550419705196509, 'at_node': [0, 39], 'negligible_lt_0.05': False}}

binding nodes, thinnest ocean:
- (87,0) zb=45.00 w=0.100 R_o-R_s=21.327 km C=0.32416
- (87,1) zb=45.00 w=0.119 R_o-R_s=21.328 km C=0.32416
- (87,2) zb=45.00 w=0.143 R_o-R_s=21.328 km C=0.32416
- (87,3) zb=45.00 w=0.170 R_o-R_s=21.328 km C=0.32416
- (87,4) zb=45.00 w=0.203 R_o-R_s=21.329 km C=0.32416
- (87,5) zb=45.00 w=0.242 R_o-R_s=21.329 km C=0.32416
- (87,9) zb=45.00 w=0.492 R_o-R_s=21.330 km C=0.32416
- (87,6) zb=45.00 w=0.289 R_o-R_s=21.330 km C=0.32416
- (87,10) zb=45.00 w=0.588 R_o-R_s=21.331 km C=0.32416
- (87,7) zb=45.00 w=0.346 R_o-R_s=21.331 km C=0.32416

binding nodes, max |dR_s|:
- (85,7) zb=44.50 w=0.346 dR_s=-12.825 km
- (85,6) zb=44.50 w=0.289 dR_s=-12.824 km
- (85,5) zb=44.50 w=0.242 dR_s=-12.823 km
- (85,4) zb=44.50 w=0.203 dR_s=-12.822 km
- (85,3) zb=44.50 w=0.170 dR_s=-12.821 km
- (85,2) zb=44.50 w=0.143 dR_s=-12.821 km
- (85,1) zb=44.50 w=0.119 dR_s=-12.820 km
- (85,0) zb=44.50 w=0.100 dR_s=-12.820 km
- (82,10) zb=43.75 w=0.588 dR_s=-12.816 km
- (82,8) zb=43.75 w=0.412 dR_s=-12.815 km

Rejects:

## rho_ice = 915, rho_s = 2200

ocean nodes 3399 (None 121, frozen/no-ocean skipped 0), feasible 3399, rejects 0 {}
- R_o - R_s km: {'min': 1.8761061303655733, 'median': 19.97864568241223, 'max': 49.532191947694464}
- derived C/MR2: {'min': 0.34192773936877635, 'median': 0.34234710991691636, 'max': 0.34680301823528076}
- R_s - R_sil(PP) km: {'min': 6.534733267763514, 'median': 8.139902702302033, 'max': 11.323522355622845}
- d_rho_o (mean[R_s,R_o]-mean[R_sil,R_o]) kg/m3: {'min': -0.27893490552105504, 'max': -0.1777836878694643}
- double-count shifts (sigma_obs; negligible < 0.05): {'C20': {'max_abs_shift_sigma_obs': 0.006989826734304268, 'at_node': [0, 39], 'negligible_lt_0.05': True}, 'C22': {'max_abs_shift_sigma_obs': 0.0052372975636484805, 'at_node': [0, 39], 'negligible_lt_0.05': True}, 'libration_deg': {'max_abs_shift_sigma_obs': 0.5989019652192903, 'at_node': [0, 39], 'negligible_lt_0.05': False}}

binding nodes, thinnest ocean:
- (87,9) zb=45.00 w=0.492 R_o-R_s=1.876 km C=0.34193
- (87,10) zb=45.00 w=0.588 R_o-R_s=1.876 km C=0.34193
- (87,11) zb=45.00 w=0.702 R_o-R_s=1.877 km C=0.34193
- (87,12) zb=45.00 w=0.838 R_o-R_s=1.877 km C=0.34193
- (87,13) zb=45.00 w=1.000 R_o-R_s=1.878 km C=0.34193
- (87,14) zb=45.00 w=1.194 R_o-R_s=1.878 km C=0.34193
- (87,0) zb=45.00 w=0.100 R_o-R_s=1.878 km C=0.34193
- (87,17) zb=45.00 w=2.031 R_o-R_s=1.878 km C=0.34193
- (87,1) zb=45.00 w=0.119 R_o-R_s=1.878 km C=0.34193
- (87,2) zb=45.00 w=0.143 R_o-R_s=1.878 km C=0.34193

binding nodes, max |dR_s|:
- (14,39) zb=19.00 w=100.000 dR_s=+11.324 km
- (17,39) zb=21.00 w=100.000 dR_s=+11.271 km
- (9,39) zb=14.00 w=100.000 dR_s=+11.157 km
- (20,39) zb=22.25 w=100.000 dR_s=+11.151 km
- (18,39) zb=21.50 w=100.000 dR_s=+11.077 km
- (11,39) zb=16.00 w=100.000 dR_s=+11.045 km
- (13,38) zb=18.00 w=83.768 dR_s=+11.034 km
- (2,39) zb=7.00 w=100.000 dR_s=+10.998 km
- (13,39) zb=18.00 w=100.000 dR_s=+10.986 km
- (15,38) zb=20.00 w=83.768 dR_s=+10.949 km

Rejects:

## rho_ice = 915, rho_s = 2600

ocean nodes 3399 (None 121, frozen/no-ocean skipped 0), feasible 3399, rejects 0 {}
- R_o - R_s km: {'min': 20.678707431859163, 'median': 38.69401664403282, 'max': 68.65708575832547}
- derived C/MR2: {'min': 0.3233294029594757, 'median': 0.32443061918381433, 'max': 0.3316484339596599}
- R_s - R_sil(PP) km: {'min': -12.183502552679856, 'median': -10.746124807591347, 'max': -8.095945335930475}
- d_rho_o (mean[R_s,R_o]-mean[R_sil,R_o]) kg/m3: {'min': 0.13845064832469234, 'max': 0.22596637006881792}
- double-count shifts (sigma_obs; negligible < 0.05): {'C20': {'max_abs_shift_sigma_obs': 0.008391396908539907, 'at_node': [0, 27], 'negligible_lt_0.05': True}, 'C22': {'max_abs_shift_sigma_obs': 0.005979461284145832, 'at_node': [0, 39], 'negligible_lt_0.05': True}, 'libration_deg': {'max_abs_shift_sigma_obs': 0.2815070787720748, 'at_node': [0, 39], 'negligible_lt_0.05': False}}

binding nodes, thinnest ocean:
- (87,0) zb=45.00 w=0.100 R_o-R_s=20.679 km C=0.32333
- (87,1) zb=45.00 w=0.119 R_o-R_s=20.679 km C=0.32333
- (87,2) zb=45.00 w=0.143 R_o-R_s=20.679 km C=0.32333
- (87,3) zb=45.00 w=0.170 R_o-R_s=20.680 km C=0.32333
- (87,4) zb=45.00 w=0.203 R_o-R_s=20.680 km C=0.32333
- (87,5) zb=45.00 w=0.242 R_o-R_s=20.681 km C=0.32333
- (87,9) zb=45.00 w=0.492 R_o-R_s=20.681 km C=0.32333
- (87,6) zb=45.00 w=0.289 R_o-R_s=20.681 km C=0.32333
- (87,10) zb=45.00 w=0.588 R_o-R_s=20.682 km C=0.32333
- (87,7) zb=45.00 w=0.346 R_o-R_s=20.682 km C=0.32333

binding nodes, max |dR_s|:
- (79,12) zb=43.00 w=0.838 dR_s=-12.184 km
- (85,7) zb=44.50 w=0.346 dR_s=-12.182 km
- (82,10) zb=43.75 w=0.588 dR_s=-12.182 km
- (82,8) zb=43.75 w=0.412 dR_s=-12.181 km
- (85,6) zb=44.50 w=0.289 dR_s=-12.181 km
- (79,11) zb=43.00 w=0.702 dR_s=-12.180 km
- (57,11) zb=33.00 w=0.702 dR_s=-12.180 km
- (85,5) zb=44.50 w=0.242 dR_s=-12.180 km
- (82,9) zb=43.75 w=0.492 dR_s=-12.180 km
- (82,7) zb=43.75 w=0.346 dR_s=-12.180 km

Rejects:

## rho_ice = 935, rho_s = 2200

ocean nodes 3399 (None 121, frozen/no-ocean skipped 0), feasible 3399, rejects 0 {}
- R_o - R_s km: {'min': 2.824528380213364, 'median': 20.62462914997156, 'max': 49.67488398092764}
- derived C/MR2: {'min': 0.34260346600799785, 'median': 0.34295573390930184, 'max': 0.34698025667927995}
- R_s - R_sil(PP) km: {'min': 5.676726075559272, 'median': 7.326927299726842, 'max': 10.82535089782605}
- d_rho_o (mean[R_s,R_o]-mean[R_sil,R_o]) kg/m3: {'min': -0.2623733882335273, 'max': -0.1618558432262489}
- double-count shifts (sigma_obs; negligible < 0.05): {'C20': {'max_abs_shift_sigma_obs': 0.006931607105685298, 'at_node': [0, 39], 'negligible_lt_0.05': True}, 'C22': {'max_abs_shift_sigma_obs': 0.005160663215266293, 'at_node': [0, 39], 'negligible_lt_0.05': True}, 'libration_deg': {'max_abs_shift_sigma_obs': 0.48933458501246746, 'at_node': [0, 39], 'negligible_lt_0.05': False}}

binding nodes, thinnest ocean:
- (87,9) zb=45.00 w=0.492 R_o-R_s=2.825 km C=0.34260
- (87,10) zb=45.00 w=0.588 R_o-R_s=2.825 km C=0.34260
- (87,11) zb=45.00 w=0.702 R_o-R_s=2.825 km C=0.34260
- (87,12) zb=45.00 w=0.838 R_o-R_s=2.826 km C=0.34260
- (87,13) zb=45.00 w=1.000 R_o-R_s=2.826 km C=0.34260
- (87,0) zb=45.00 w=0.100 R_o-R_s=2.826 km C=0.34260
- (87,1) zb=45.00 w=0.119 R_o-R_s=2.826 km C=0.34260
- (87,2) zb=45.00 w=0.143 R_o-R_s=2.826 km C=0.34260
- (87,3) zb=45.00 w=0.170 R_o-R_s=2.827 km C=0.34260
- (87,4) zb=45.00 w=0.203 R_o-R_s=2.827 km C=0.34260

binding nodes, max |dR_s|:
- (14,39) zb=19.00 w=100.000 dR_s=+10.825 km
- (2,39) zb=7.00 w=100.000 dR_s=+10.800 km
- (9,39) zb=14.00 w=100.000 dR_s=+10.777 km
- (17,39) zb=21.00 w=100.000 dR_s=+10.727 km
- (4,39) zb=9.00 w=100.000 dR_s=+10.679 km
- (11,39) zb=16.00 w=100.000 dR_s=+10.618 km
- (20,39) zb=22.25 w=100.000 dR_s=+10.578 km
- (8,38) zb=13.00 w=83.768 dR_s=+10.569 km
- (13,38) zb=18.00 w=83.768 dR_s=+10.566 km
- (6,39) zb=11.00 w=100.000 dR_s=+10.566 km

Rejects:

## rho_ice = 935, rho_s = 2600

ocean nodes 3399 (None 121, frozen/no-ocean skipped 0), feasible 3399, rejects 0 {}
- R_o - R_s km: {'min': 21.539979383346218, 'median': 39.30233198138571, 'max': 68.78596486012684}
- derived C/MR2: {'min': 0.32443079205002545, 'median': 0.32531931140487846, 'max': 0.33188032030871084}
- R_s - R_sil(PP) km: {'min': -13.036132457706291, 'median': -11.425360155618925, 'max': -8.350371495879546}
- d_rho_o (mean[R_s,R_o]-mean[R_sil,R_o]) kg/m3: {'min': 0.14246231376978358, 'max': 0.23277162226406745}
- double-count shifts (sigma_obs; negligible < 0.05): {'C20': {'max_abs_shift_sigma_obs': 0.008769754282689734, 'at_node': [0, 16], 'negligible_lt_0.05': True}, 'C22': {'max_abs_shift_sigma_obs': 0.006035749346038974, 'at_node': [0, 39], 'negligible_lt_0.05': True}, 'libration_deg': {'max_abs_shift_sigma_obs': 0.24607044810949608, 'at_node': [0, 39], 'negligible_lt_0.05': False}}

binding nodes, thinnest ocean:
- (87,0) zb=45.00 w=0.100 R_o-R_s=21.540 km C=0.32443
- (87,1) zb=45.00 w=0.119 R_o-R_s=21.540 km C=0.32443
- (87,2) zb=45.00 w=0.143 R_o-R_s=21.541 km C=0.32443
- (87,3) zb=45.00 w=0.170 R_o-R_s=21.541 km C=0.32443
- (87,4) zb=45.00 w=0.203 R_o-R_s=21.541 km C=0.32443
- (87,5) zb=45.00 w=0.242 R_o-R_s=21.542 km C=0.32443
- (87,9) zb=45.00 w=0.492 R_o-R_s=21.542 km C=0.32443
- (87,6) zb=45.00 w=0.289 R_o-R_s=21.543 km C=0.32443
- (87,7) zb=45.00 w=0.346 R_o-R_s=21.544 km C=0.32443
- (87,10) zb=45.00 w=0.588 R_o-R_s=21.544 km C=0.32443

binding nodes, max |dR_s|:
- (85,7) zb=44.50 w=0.346 dR_s=-13.036 km
- (85,6) zb=44.50 w=0.289 dR_s=-13.035 km
- (85,5) zb=44.50 w=0.242 dR_s=-13.034 km
- (85,4) zb=44.50 w=0.203 dR_s=-13.033 km
- (85,3) zb=44.50 w=0.170 dR_s=-13.032 km
- (85,2) zb=44.50 w=0.143 dR_s=-13.032 km
- (85,1) zb=44.50 w=0.119 dR_s=-13.031 km
- (85,0) zb=44.50 w=0.100 dR_s=-13.031 km
- (82,10) zb=43.75 w=0.588 dR_s=-13.025 km
- (82,8) zb=43.75 w=0.412 dR_s=-13.024 km

Rejects:
