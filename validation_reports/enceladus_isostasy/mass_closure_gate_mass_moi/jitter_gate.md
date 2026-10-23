# Enceladus mass-closure jitter gate

Overall: **FAIL**

builds=90 deltaP/[1, 2, 4] nPoros=[8, 12, 16] mean build 76.53 s, total 1216.1 s; cache=/Users/svance/scratch/enceladus_zbw_production_cache.pkl

## Criteria

- C20_spread<=0.1sigma_obs_every_node: True
- C22_spread<=0.1sigma_obs_every_node: True
- every_node_build_feasible: True
- closed_stack_residual<=1e-12: True
- |offset - predicted|<=0.1sigma_every_node: True
- libration_spread<=0.05sigma_all_nodes: True
- local_quadratic<=3x_jitter_std: False
- posterior_median_shift<0.2_sigma: TODO: later step (sigma_model on vs off posterior runs)

## sigma_model = max(std, |bias|)

{
 "C20": {
  "max_std_over_sigma_obs": 0.00011893148387580932,
  "max_abs_bias_over_sigma_obs": 0.00028932456189413916,
  "sigma_model_over_sigma_obs": 0.00028932456189413916,
  "sigma_model_abs": 1.0702115544464208e-08,
  "recommended_sigma_model_add_abs": 5.300010805203809e-06,
  "existing_sigma_model_add": 5.3e-06
 },
 "C22": {
  "max_std_over_sigma_obs": 3.4525907189968694e-05,
  "max_abs_bias_over_sigma_obs": 7.915507539567409e-05,
  "sigma_model_over_sigma_obs": 7.915507539567409e-05,
  "sigma_model_abs": 1.163579608316409e-09,
  "recommended_sigma_model_add_abs": 1.7000003982109843e-06,
  "existing_sigma_model_add": 1.7e-06
 },
 "C30": {
  "max_std_over_sigma_obs": 9.501738942987449e-05,
  "max_abs_bias_over_sigma_obs": 0.0002295386586437231,
  "sigma_model_over_sigma_obs": 0.0002295386586437231,
  "sigma_model_abs": 7.671181971873227e-09,
  "recommended_sigma_model_add_abs": 4.400006687157742e-06,
  "existing_sigma_model_add": 4.4e-06
 },
 "libration_deg": {
  "max_std_over_sigma_obs": 0.00045316017549863833,
  "max_abs_bias_over_sigma_obs": 0.0006398551730957269,
  "sigma_model_over_sigma_obs": 0.0006398551730957269,
  "sigma_model_abs": 1.9195655192871808e-06,
  "recommended_sigma_model_add_abs": 0.0002500073693549509,
  "existing_sigma_model_add": 0.00025
 }
}

## C20-C22 jitter correlation

{"rho_obs_in_config": 0.47, "pooled": 0.9982027749338784, "per_node": {"(70,18)": 0.9999999999034058, "(77,10)": 0.9999999999439297, "(12,6)": 0.9999999998471523, "(25,30)": 0.9999999990479397, "(38,18)": 0.9999999996181993, "(50,34)": 0.9999999995824819, "(58,12)": 0.9999999997710215, "(64,26)": 0.9999999997800783, "(84,4)": 0.9999999999558693, "(86,22)": 0.9999999998630666}, "flag_differs_strongly": true, "note": "runner applies rho_obs to the INFLATED sigmas; flag = |pooled - rho_obs| > 0.3 (reporting only, not gating)"}

## Failed builds

{}

## Bias (production - finest), sigma_obs

- (70,18): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 9.197555402499849e-05, "C22": 2.4579748608265333e-05, "C30": -6.8645103890984e-05, "libration_deg": -0.0001310214366769443}
- (77,10): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.00021894430147766583, "C22": 5.7950685922530525e-05, "C30": -0.00016125534602569653, "libration_deg": -0.000277642574144273}
- (12,6): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -9.252962331976124e-05, "C22": -2.5411310846079164e-05, "C30": 7.66838520555273e-05, "libration_deg": 0.0003545699624740717}
- (25,30): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 7.314839925652366e-05, "C22": 2.1233240004447097e-05, "C30": -5.844068442781119e-05, "libration_deg": -0.000278545856963186}
- (38,18): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -0.00028932456189413916, "C22": -7.915507539567409e-05, "C30": 0.0002295386586437231, "libration_deg": 0.0006398551730957269}
- (50,34): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.000106068696192222, "C22": 3.188872155929922e-05, "C30": -8.175648375994281e-05, "libration_deg": -0.00044035679849322806}
- (58,12): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -1.707398197604901e-05, "C22": -4.588585865203293e-06, "C30": 1.313394268389828e-05, "libration_deg": 2.797186352913409e-05}
- (64,26): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.00013415179386087726, "C22": 3.6812884577066936e-05, "C30": -0.00010110126870474425, "libration_deg": -0.0002458268505746893}
- (84,4): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -7.392584043774203e-05, "C22": -1.950578230564144e-05, "C30": 5.400125901688994e-05, "libration_deg": 8.93718912625685e-05}
- (86,22): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.00012900950299839228, "C22": 3.448262576719815e-05, "C30": -9.367067892567164e-05, "libration_deg": -0.00017208284302146679}

## mass_moi closure per node (production build = first, dP/1 nPoros 8 if present; range over builds in brackets)

| node | feasible | R_s-R_sil(PP) km | rho_s kg/m3 | rho_s-rho_sil(PP) | MoI stack before | MoI after (target) | PP CMR2 | max abs M resid |
|---|---|---|---|---|---|---|---|---|
| (70,18) | True | -1.297 [-1.344, -0.513] | 2342.5 [2342.471, 2342.473] | 43.0 | 0.33673 [0.336, 0.337] | 0.33500 | 0.33528202108895777 | 3.33e-16 |
| (77,10) | True | -1.247 [-1.247, -0.464] | 2341.0 [2340.992, 2340.993] | 41.9 | 0.33668 [0.336, 0.337] | 0.33500 | 0.3352430354360986 | 3.33e-16 |
| (12,6) | True | -0.428 [-1.315, -0.428] | 2368.3 [2368.325, 2368.335] | 18.7 | 0.33567 [0.336, 0.336] | 0.33500 | 0.33505250310960405 | 1.11e-16 |
| (25,30) | True | -1.130 [-1.457, 1.185] | 2362.4 [2362.403, 2362.424] | 33.3 | 0.33634 [0.334, 0.336] | 0.33500 | 0.3356601931498706 | 1.11e-16 |
| (38,18) | True | 0.289 [-1.353, 0.289] | 2353.6 [2353.603, 2353.612] | 3.2 | 0.33498 [0.335, 0.337] | 0.33500 | 0.33440416732335076 | 0.00e+00 |
| (50,34) | True | -1.689 [-1.689, -0.171] | 2355.9 [2355.923, 2355.931] | 48.5 | 0.33697 [0.335, 0.337] | 0.33500 | 0.33583032998234635 | 3.33e-16 |
| (58,12) | True | -0.421 [-1.275, 0.361] | 2346.4 [2346.448, 2346.452] | 21.3 | 0.33575 [0.335, 0.336] | 0.33500 | 0.33479730866205515 | 3.33e-16 |
| (64,26) | True | -1.575 [-1.575, -0.017] | 2344.8 [2344.835, 2344.839] | 43.8 | 0.33686 [0.335, 0.337] | 0.33500 | 0.3359028522454197 | 1.11e-16 |
| (84,4) | True | -0.195 [-1.210, -0.195] | 2340.3 [2340.346, 2340.346] | 15.0 | 0.33549 [0.335, 0.336] | 0.33500 | 0.33469053336405263 | 3.33e-16 |
| (86,22) | True | -1.366 [-1.796, 0.195] | 2340.3 [2340.275, 2340.276] | 40.2 | 0.33669 [0.335, 0.337] | 0.33500 | 0.33566102083536614 | 1.11e-16 |

closure infeasible reasons: {"(70,18)": [], "(77,10)": [], "(12,6)": [], "(25,30)": [], "(38,18)": [], "(50,34)": [], "(58,12)": [], "(64,26)": [], "(84,4)": [], "(86,22)": []}

## Offset mass_moi - rock_density (finest deltaP), measured vs predicted = (moi_target - MoI_rd) dX/dC, sigma_obs

- C20: largest |offset| 0.9791002472599699, max |meas-pred| 0.003673943013029257
    - (70,18): signed offset +0.7477, predicted +0.7455, residual +0.0022
    - (77,10): signed offset +0.6069, predicted +0.6055, residual +0.0014
    - (12,6): signed offset +0.5881, predicted +0.5869, residual +0.0012
    - (25,30): signed offset +0.8307, predicted +0.8282, residual +0.0025
    - (38,18): signed offset +0.7948, predicted +0.7926, residual +0.0023
    - (50,34): signed offset +0.7660, predicted +0.7639, residual +0.0021
    - (58,12): signed offset +0.7501, predicted +0.7480, residual +0.0021
    - (64,26): signed offset +0.7385, predicted +0.7365, residual +0.0020
    - (84,4): signed offset +0.6940, predicted +0.6922, residual +0.0019
    - (86,22): signed offset +0.9791, predicted +0.9754, residual +0.0037
- C22: largest |offset| 0.7457443842049789, max |meas-pred| 0.0028964443103092208
    - (70,18): signed offset -0.5689, predicted -0.5672, residual -0.0017
    - (77,10): signed offset -0.4619, predicted -0.4608, residual -0.0011
    - (12,6): signed offset -0.4456, predicted -0.4447, residual -0.0009
    - (25,30): signed offset -0.6311, predicted -0.6291, residual -0.0020
    - (38,18): signed offset -0.6035, predicted -0.6017, residual -0.0018
    - (50,34): signed offset -0.5831, predicted -0.5814, residual -0.0017
    - (58,12): signed offset -0.5701, predicted -0.5685, residual -0.0016
    - (64,26): signed offset -0.5619, predicted -0.5603, residual -0.0016
    - (84,4): signed offset -0.5284, predicted -0.5269, residual -0.0015
    - (86,22): signed offset -0.7457, predicted -0.7428, residual -0.0029


## Frozen two-layer stack (M_stack-M)/M
{"n": 39, "min": -3.03345614782174e-16, "max": 3.03345614782174e-16, "median": 0.0, "min_abs": 0.0, "max_abs": 3.03345614782174e-16}


## Per node: spread / std in sigma_obs (closure ON | OFF)

| node (i,j) | zb | w | n | C20 sp | C22 sp | C30 sp | lib sp | C20 sp OFF | C22 sp OFF | lib sp OFF |
|---|---|---|---|---|---|---|---|---|---|---|
| (70,18) | 39.50 | 2.424 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.229 | 0.891 | 0.068 |
| (77,10) | 42.50 | 0.588 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.256 | 0.902 | 0.065 |
| (12,6) | 17.00 | 0.289 | 9 | 0.000 | 0.000 | 0.000 | 0.001 | 0.817 | 0.596 | 0.157 |
| (25,30) | 23.50 | 20.309 | 9 | 0.000 | 0.000 | 0.000 | 0.002 | 1.547 | 1.164 | 0.154 |
| (38,18) | 26.75 | 2.424 | 9 | 0.000 | 0.000 | 0.000 | 0.001 | 1.415 | 1.052 | 0.119 |
| (50,34) | 29.75 | 41.246 | 9 | 0.000 | 0.000 | 0.000 | 0.001 | 1.701 | 1.262 | 0.149 |
| (58,12) | 33.50 | 0.838 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.638 | 1.210 | 0.103 |
| (64,26) | 36.50 | 10.000 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.704 | 1.258 | 0.100 |
| (84,4) | 44.25 | 0.203 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 0.892 | 0.663 | 0.036 |
| (86,22) | 44.75 | 4.924 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.660 | 1.225 | 0.070 |

## Local quadratic check (cache neighbourhood, closure ON)

- (70,18): {"C20": {"dep_sigma_obs": 9.146386455469101e-05, "jitter_std_sigma_obs": 9.118523540174588e-05, "ratio_dep_over_jitter_std": 1.0030556389059866, "pass": true}, "C22": {"dep_sigma_obs": 4.2843159398576996e-05, "jitter_std_sigma_obs": 2.436950843130319e-05, "ratio_dep_over_jitter_std": 1.7580641611770873, "pass": true}, "C30": {"dep_sigma_obs": 8.25000945500398e-05, "jitter_std_sigma_obs": 6.805458818220461e-05, "ratio_dep_over_jitter_std": 1.2122635189439364, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.00024541504189494007, "jitter_std_sigma_obs": 0.00012993618096006617, "ratio_dep_over_jitter_std": 1.8887352243357414, "pass": true}}
- (77,10): {"C20": {"dep_sigma_obs": 0.00015518108126992938, "jitter_std_sigma_obs": 9.077141785864286e-05, "ratio_dep_over_jitter_std": 1.7095808893454858, "pass": true}, "C22": {"dep_sigma_obs": 5.8274841425612745e-06, "jitter_std_sigma_obs": 2.4025611889685633e-05, "ratio_dep_over_jitter_std": 0.24255299591612295, "pass": true}, "C30": {"dep_sigma_obs": 8.715920823721974e-05, "jitter_std_sigma_obs": 6.685432215747888e-05, "ratio_dep_over_jitter_std": 1.3037183748854955, "pass": true}, "libration_deg": {"dep_sigma_obs": 2.7523497636109763e-05, "jitter_std_sigma_obs": 0.00011510779254376789, "ratio_dep_over_jitter_std": 0.23911063732409252, "pass": true}}
- (12,6): {"C20": {"dep_sigma_obs": 0.00018193425867268692, "jitter_std_sigma_obs": 7.714171663696156e-05, "ratio_dep_over_jitter_std": 2.358441924865271, "pass": true}, "C22": {"dep_sigma_obs": 4.24936866389364e-05, "jitter_std_sigma_obs": 2.1184765198137397e-05, "ratio_dep_over_jitter_std": 2.0058606381284094, "pass": true}, "C30": {"dep_sigma_obs": 7.329912308547531e-05, "jitter_std_sigma_obs": 6.39314365577331e-05, "ratio_dep_over_jitter_std": 1.1465270770082376, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.005633721588775102, "jitter_std_sigma_obs": 0.00029554778442724165, "ratio_dep_over_jitter_std": 19.061965223975548, "pass": false}}
- (25,30): {"C20": {"dep_sigma_obs": 0.0001853410025696878, "jitter_std_sigma_obs": 0.00011893148387580932, "ratio_dep_over_jitter_std": 1.5583846810759094, "pass": true}, "C22": {"dep_sigma_obs": 4.362773575122995e-05, "jitter_std_sigma_obs": 3.4525907189968694e-05, "ratio_dep_over_jitter_std": 1.2636231543803065, "pass": true}, "C30": {"dep_sigma_obs": 0.00014033028244445753, "jitter_std_sigma_obs": 9.501738942987449e-05, "ratio_dep_over_jitter_std": 1.4768905280019848, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.0004768411999006572, "jitter_std_sigma_obs": 0.00045316017549863833, "ratio_dep_over_jitter_std": 1.0522575144118993, "pass": true}}
- (38,18): {"C20": {"dep_sigma_obs": 0.00015400001573086266, "jitter_std_sigma_obs": 9.58416669646102e-05, "ratio_dep_over_jitter_std": 1.606816957677996, "pass": true}, "C22": {"dep_sigma_obs": 3.880476975982088e-05, "jitter_std_sigma_obs": 2.622069750808415e-05, "ratio_dep_over_jitter_std": 1.4799289663387831, "pass": true}, "C30": {"dep_sigma_obs": 0.00011937094115215753, "jitter_std_sigma_obs": 7.603710046951075e-05, "ratio_dep_over_jitter_std": 1.569903907632863, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.00030175633793557327, "jitter_std_sigma_obs": 0.0002119460570014672, "ratio_dep_over_jitter_std": 1.4237412207837599, "pass": true}}
- (50,34): {"C20": {"dep_sigma_obs": 0.00011227552511189082, "jitter_std_sigma_obs": 5.0658062898189724e-05, "ratio_dep_over_jitter_std": 2.2163406709320306, "pass": true}, "C22": {"dep_sigma_obs": 4.470186720489121e-05, "jitter_std_sigma_obs": 1.5230287458440196e-05, "ratio_dep_over_jitter_std": 2.935063919631976, "pass": true}, "C30": {"dep_sigma_obs": 9.42354663519344e-05, "jitter_std_sigma_obs": 3.9046555913095144e-05, "ratio_dep_over_jitter_std": 2.4134130180820996, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.0006181754583824602, "jitter_std_sigma_obs": 0.00021034859132795225, "ratio_dep_over_jitter_std": 2.938814348505284, "pass": true}}
- (58,12): {"C20": {"dep_sigma_obs": 1.049816820465821e-05, "jitter_std_sigma_obs": 8.294298051362797e-05, "ratio_dep_over_jitter_std": 0.1265709061773263, "pass": true}, "C22": {"dep_sigma_obs": 2.8113747334953237e-06, "jitter_std_sigma_obs": 2.2291148773257847e-05, "ratio_dep_over_jitter_std": 0.12612067516538503, "pass": true}, "C30": {"dep_sigma_obs": 7.919714206981315e-06, "jitter_std_sigma_obs": 6.380254616153615e-05, "ratio_dep_over_jitter_std": 0.1241284977394174, "pass": true}, "libration_deg": {"dep_sigma_obs": 1.842633794620537e-05, "jitter_std_sigma_obs": 0.00013590368842221725, "ratio_dep_over_jitter_std": 0.13558379584930433, "pass": true}}
- (64,26): {"C20": {"dep_sigma_obs": 7.381891525292632e-05, "jitter_std_sigma_obs": 8.849115682180374e-05, "ratio_dep_over_jitter_std": 0.8341953919935392, "pass": true}, "C22": {"dep_sigma_obs": 5.869222744400489e-05, "jitter_std_sigma_obs": 2.428412396253233e-05, "ratio_dep_over_jitter_std": 2.416897044940159, "pass": true}, "C30": {"dep_sigma_obs": 8.485202011383336e-05, "jitter_std_sigma_obs": 6.668923147767654e-05, "ratio_dep_over_jitter_std": 1.2723496467677335, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.0004196260196450326, "jitter_std_sigma_obs": 0.00016221364114290266, "ratio_dep_over_jitter_std": 2.5868725755028312, "pass": true}}
- (84,4): {"C20": {"dep_sigma_obs": 0.0001370763910176124, "jitter_std_sigma_obs": 4.071908280823843e-05, "ratio_dep_over_jitter_std": 3.36639191170285, "pass": false}, "C22": {"dep_sigma_obs": 3.6170638273332294e-05, "jitter_std_sigma_obs": 1.0744205607132186e-05, "ratio_dep_over_jitter_std": 3.3665251388452218, "pass": false}, "C30": {"dep_sigma_obs": 0.00010011077118761819, "jitter_std_sigma_obs": 2.9744258554292958e-05, "ratio_dep_over_jitter_std": 3.365717488129127, "pass": false}, "libration_deg": {"dep_sigma_obs": 0.00016554431728673288, "jitter_std_sigma_obs": 4.9235212108153794e-05, "ratio_dep_over_jitter_std": 3.3623155095399144, "pass": false}}
- (86,22): {"C20": {"dep_sigma_obs": 0.0002236592770246481, "jitter_std_sigma_obs": 7.975147306100152e-05, "ratio_dep_over_jitter_std": 2.8044532400495243, "pass": true}, "C22": {"dep_sigma_obs": 0.00010179575269867734, "jitter_std_sigma_obs": 2.1317356026223598e-05, "ratio_dep_over_jitter_std": 4.775252267375609, "pass": false}, "C30": {"dep_sigma_obs": 4.710331784189017e-05, "jitter_std_sigma_obs": 5.790504422550518e-05, "ratio_dep_over_jitter_std": 0.8134579374199455, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.000599946817905761, "jitter_std_sigma_obs": 0.00010641114362804576, "ratio_dep_over_jitter_std": 5.638007425264047, "pass": false}}

TODO (later step): posterior-median check, sigma_model on vs off (<0.2 posterior sigma).
