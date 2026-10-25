# Enceladus mass-closure jitter gate

Overall: **PASS (all rho_s; all but posterior-median check, TODO)**

builds=90 deltaP/[1, 2, 4] nPoros=[8, 12, 16] mean build 80.2 s, total 1268.9 s; cache=/Users/svance/scratch/enceladus_zbw_production_cache.pkl

## Criteria

- rho_s=2200: C20_spread<=0.1sigma_obs_every_node: True
- rho_s=2200: C22_spread<=0.1sigma_obs_every_node: True
- rho_s=2200: every_node_build_feasible: True
- rho_s=2200: closed_stack_residual<=1e-12: True
- rho_s=2200: libration_spread<=0.05sigma_all_nodes: True
- rho_s=2200: local_quadratic<=max(3x_jitter_std,0.01sigma_obs): True
- rho_s=2200: posterior_median_shift<0.2_sigma: TODO: later step (sigma_model on vs off posterior runs)
- rho_s=2400: C20_spread<=0.1sigma_obs_every_node: True
- rho_s=2400: C22_spread<=0.1sigma_obs_every_node: True
- rho_s=2400: every_node_build_feasible: True
- rho_s=2400: closed_stack_residual<=1e-12: True
- rho_s=2400: libration_spread<=0.05sigma_all_nodes: True
- rho_s=2400: local_quadratic<=max(3x_jitter_std,0.01sigma_obs): True
- rho_s=2400: posterior_median_shift<0.2_sigma: TODO: later step (sigma_model on vs off posterior runs)
- rho_s=2600: C20_spread<=0.1sigma_obs_every_node: True
- rho_s=2600: C22_spread<=0.1sigma_obs_every_node: True
- rho_s=2600: every_node_build_feasible: True
- rho_s=2600: closed_stack_residual<=1e-12: True
- rho_s=2600: libration_spread<=0.05sigma_all_nodes: True
- rho_s=2600: local_quadratic<=max(3x_jitter_std,0.01sigma_obs): True
- rho_s=2600: posterior_median_shift<0.2_sigma: TODO: later step (sigma_model on vs off posterior runs)

## sigma_model = max(std, |bias|)

{
 "C20": {
  "max_std_over_sigma_obs": 0.00039067843207510827,
  "max_abs_bias_over_sigma_obs": 0.0006602298543879599,
  "sigma_model_over_sigma_obs": 0.0006602298543879599,
  "sigma_model_abs": 2.442190231381064e-08,
  "recommended_sigma_model_add_abs": 5.300056266617613e-06,
  "existing_sigma_model_add": 5.3e-06
 },
 "C22": {
  "max_std_over_sigma_obs": 0.00016558063222391517,
  "max_abs_bias_over_sigma_obs": 0.00021329956213515412,
  "sigma_model_over_sigma_obs": 0.00021329956213515412,
  "sigma_model_abs": 3.1355035633867656e-09,
  "recommended_sigma_model_add_abs": 1.7000028915806572e-06,
  "existing_sigma_model_add": 1.7e-06
 },
 "C30": {
  "max_std_over_sigma_obs": 9.888812470843593e-05,
  "max_abs_bias_over_sigma_obs": 0.00023902002784763698,
  "sigma_model_over_sigma_obs": 0.00023902002784763698,
  "sigma_model_abs": 7.988049330668028e-09,
  "recommended_sigma_model_add_abs": 4.400007251009038e-06,
  "existing_sigma_model_add": 4.4e-06
 },
 "libration_deg": {
  "max_std_over_sigma_obs": 0.0004560242417406011,
  "max_abs_bias_over_sigma_obs": 0.0006388222979492042,
  "sigma_model_over_sigma_obs": 0.0006388222979492042,
  "sigma_model_abs": 1.9164668938476126e-06,
  "recommended_sigma_model_add_abs": 0.00025000734558279524,
  "existing_sigma_model_add": 0.00025
 }
}

## C20-C22 jitter correlation

{"rho_obs_in_config": 0.47, "pooled": -0.9883801343637441, "per_node": {"(70,18)": -0.9999999932307817, "(77,10)": -0.9999999889201727, "(12,6)": -0.9999999997533298, "(25,30)": -0.9999999991788635, "(38,18)": -0.9999999984019166, "(50,34)": -0.9999999996696802, "(58,12)": -0.9999999962025097, "(64,26)": -0.999999995184182, "(84,4)": -0.9999999818555192, "(86,22)": -0.9999999430801347}, "flag_differs_strongly": true, "note": "runner applies rho_obs to the INFLATED sigmas; flag = |pooled - rho_obs| > 0.3 (reporting only, not gating)"}

## Failed builds

{}

## Bias (production - finest), sigma_obs

- (70,18): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.00013759619150432168, "C22": -2.48746908074722e-05, "C30": -7.152262416636016e-05, "libration_deg": -0.00012907441565568334}
- (77,10): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.0002958768695996653, "C22": -4.041725327612294e-05, "C30": -0.00016804484592390734, "libration_deg": -0.0002728969657519083}
- (12,6): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -0.0002797286610511508, "C22": 0.00011148777548769547, "C30": 7.978326076708154e-05, "libration_deg": 0.00036064452576164047}
- (25,30): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.00024018011913597514, "C22": -0.00010177047071187942, "C30": -6.0821390469989945e-05, "libration_deg": -0.0002803048823838787}
- (38,18): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -0.0006602298543879599, "C22": 0.00021329956213515412, "C30": 0.00023902002784763698, "libration_deg": 0.0006388222979492042}
- (50,34): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.00035169093499432967, "C22": -0.00015171217760785548, "C30": -8.50928515579414e-05, "libration_deg": -0.0004397924740680755}
- (58,12): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -3.057441907693546e-05, "C22": 7.475908488847628e-06, "C30": 1.3682904824200905e-05, "libration_deg": 2.770080200405663e-05}
- (64,26): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.00023654651650674514, "C22": -5.7393628215689664e-05, "C30": -0.0001052895350797566, "libration_deg": -0.00024288980925620446}
- (84,4): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -9.503996848698456e-05, "C22": 1.0720695268202576e-05, "C30": 5.627779387214614e-05, "libration_deg": 8.774228258308832e-05}
- (86,22): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.00016702101510206977, "C22": -1.9505211115484666e-05, "C30": -9.758927253861733e-05, "libration_deg": -0.00016890813501803845}

## mass_rho_rock closure per node and rho_s (range over builds in brackets; spreads per rho_s are in the criteria; the per-node spread table below is the middle rho_s)

### rho_s = 2200 kg/m3: PASS (all but posterior-median check, TODO)

| node | feasible | R_s km | R_o-R_s km | dR_s vs PP km | derived C/MR2 | C20 sp | C22 sp | lib sp |
|---|---|---|---|---|---|---|---|---|
| (70,18) | True | [204.153, 204.153] | [8.449, 8.449] | [6.175, 7.006] | [0.34248, 0.34248] | 0.0002 | 0.0000 | 0.0004 |
| (77,10) | True | [204.356, 204.356] | [5.246, 5.246] | [6.199, 6.982] | [0.34245, 0.34245] | 0.0002 | 0.0000 | 0.0003 |
| (12,6) | True | [202.534, 202.535] | [32.565, 32.566] | [7.369, 8.258] | [0.34336, 0.34336] | 0.0005 | 0.0002 | 0.0010 |
| (25,30) | True | [202.648, 202.649] | [25.951, 25.952] | [7.058, 9.701] | [0.34311, 0.34311] | 0.0009 | 0.0003 | 0.0016 |
| (38,18) | True | [203.249, 203.250] | [22.100, 22.101] | [6.671, 8.314] | [0.34284, 0.34284] | 0.0005 | 0.0001 | 0.0008 |
| (50,34) | True | [202.914, 202.915] | [19.436, 19.436] | [6.633, 8.152] | [0.34286, 0.34286] | 0.0003 | 0.0001 | 0.0007 |
| (58,12) | True | [203.758, 203.759] | [14.842, 14.842] | [6.416, 8.052] | [0.34260, 0.34260] | 0.0003 | 0.0000 | 0.0005 |
| (64,26) | True | [203.887, 203.888] | [11.714, 11.714] | [6.085, 7.643] | [0.34254, 0.34254] | 0.0003 | 0.0000 | 0.0005 |
| (84,4) | True | [204.466, 204.466] | [3.388, 3.388] | [6.207, 7.222] | [0.34244, 0.34244] | 0.0001 | 0.0000 | 0.0002 |
| (86,22) | True | [204.488, 204.488] | [2.865, 2.865] | [5.640, 7.631] | [0.34244, 0.34244] | 0.0002 | 0.0000 | 0.0004 |

### rho_s = 2400 kg/m3: PASS (all but posterior-median check, TODO)

| node | feasible | R_s km | R_o-R_s km | dR_s vs PP km | derived C/MR2 | C20 sp | C22 sp | lib sp |
|---|---|---|---|---|---|---|---|---|
| (70,18) | True | [193.894, 193.894] | [18.707, 18.708] | [-4.084, -3.253] | [0.33234, 0.33234] | 0.0004 | 0.0001 | 0.0003 |
| (77,10) | True | [194.097, 194.098] | [15.504, 15.504] | [-4.060, -3.276] | [0.33226, 0.33226] | 0.0003 | 0.0000 | 0.0003 |
| (12,6) | True | [192.374, 192.375] | [42.725, 42.726] | [-2.791, -1.902] | [0.33362, 0.33362] | 0.0007 | 0.0003 | 0.0009 |
| (25,30) | True | [192.356, 192.358] | [36.242, 36.244] | [-3.234, -0.591] | [0.33335, 0.33335] | 0.0013 | 0.0006 | 0.0015 |
| (38,18) | True | [193.038, 193.039] | [32.311, 32.312] | [-3.540, -1.897] | [0.33292, 0.33292] | 0.0007 | 0.0002 | 0.0007 |
| (50,34) | True | [192.473, 192.474] | [29.876, 29.877] | [-3.807, -2.288] | [0.33305, 0.33305] | 0.0005 | 0.0002 | 0.0006 |
| (58,12) | True | [193.530, 193.531] | [25.070, 25.071] | [-3.813, -2.177] | [0.33256, 0.33256] | 0.0005 | 0.0001 | 0.0004 |
| (64,26) | True | [193.595, 193.596] | [22.005, 22.006] | [-4.207, -2.649] | [0.33248, 0.33248] | 0.0005 | 0.0001 | 0.0005 |
| (84,4) | True | [194.204, 194.204] | [13.650, 13.650] | [-4.055, -3.040] | [0.33222, 0.33222] | 0.0002 | 0.0000 | 0.0002 |
| (86,22) | True | [194.195, 194.196] | [13.157, 13.157] | [-4.653, -2.662] | [0.33222, 0.33222] | 0.0003 | 0.0000 | 0.0003 |

### rho_s = 2600 kg/m3: PASS (all but posterior-median check, TODO)

| node | feasible | R_s km | R_o-R_s km | dR_s vs PP km | derived C/MR2 | C20 sp | C22 sp | lib sp |
|---|---|---|---|---|---|---|---|---|
| (70,18) | True | [185.428, 185.429] | [27.173, 27.173] | [-12.549, -11.718] | [0.32437, 0.32437] | 0.0005 | 0.0001 | 0.0003 |
| (77,10) | True | [185.631, 185.631] | [23.970, 23.971] | [-12.526, -11.743] | [0.32424, 0.32424] | 0.0004 | 0.0001 | 0.0003 |
| (12,6) | True | [183.988, 183.989] | [51.111, 51.112] | [-11.177, -10.289] | [0.32595, 0.32595] | 0.0009 | 0.0004 | 0.0009 |
| (25,30) | True | [183.880, 183.883] | [44.717, 44.720] | [-11.709, -9.067] | [0.32570, 0.32570] | 0.0017 | 0.0008 | 0.0015 |
| (38,18) | True | [184.611, 184.613] | [40.738, 40.739] | [-11.966, -10.324] | [0.32512, 0.32512] | 0.0010 | 0.0003 | 0.0007 |
| (50,34) | True | [183.896, 183.897] | [38.453, 38.455] | [-12.385, -10.866] | [0.32539, 0.32539] | 0.0007 | 0.0003 | 0.0006 |
| (58,12) | True | [185.088, 185.089] | [33.511, 33.512] | [-12.254, -10.619] | [0.32466, 0.32466] | 0.0007 | 0.0002 | 0.0004 |
| (64,26) | True | [185.110, 185.111] | [30.491, 30.492] | [-12.692, -11.134] | [0.32457, 0.32457] | 0.0007 | 0.0002 | 0.0005 |
| (84,4) | True | [185.734, 185.735] | [22.120, 22.120] | [-12.524, -11.510] | [0.32418, 0.32418] | 0.0003 | 0.0000 | 0.0002 |
| (86,22) | True | [185.705, 185.705] | [21.647, 21.648] | [-13.143, -11.152] | [0.32418, 0.32418] | 0.0005 | 0.0001 | 0.0003 |


## Frozen two-layer stack (M_stack-M)/M
{"n": 39, "min": -3.03345614782174e-16, "max": 3.03345614782174e-16, "median": 0.0, "min_abs": 0.0, "max_abs": 3.03345614782174e-16}


## Per node: spread / std in sigma_obs (closure ON | OFF)

| node (i,j) | zb | w | n | C20 sp | C22 sp | C30 sp | lib sp | C20 sp OFF | C22 sp OFF | lib sp OFF |
|---|---|---|---|---|---|---|---|---|---|---|
| (70,18) | 39.50 | 2.424 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.229 | 0.891 | 0.068 |
| (77,10) | 42.50 | 0.588 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.256 | 0.902 | 0.065 |
| (12,6) | 17.00 | 0.289 | 9 | 0.001 | 0.000 | 0.000 | 0.001 | 0.817 | 0.596 | 0.157 |
| (25,30) | 23.50 | 20.309 | 9 | 0.001 | 0.001 | 0.000 | 0.002 | 1.547 | 1.164 | 0.154 |
| (38,18) | 26.75 | 2.424 | 9 | 0.001 | 0.000 | 0.000 | 0.001 | 1.415 | 1.052 | 0.119 |
| (50,34) | 29.75 | 41.246 | 9 | 0.001 | 0.000 | 0.000 | 0.001 | 1.701 | 1.262 | 0.149 |
| (58,12) | 33.50 | 0.838 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.638 | 1.210 | 0.103 |
| (64,26) | 36.50 | 10.000 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.704 | 1.258 | 0.100 |
| (84,4) | 44.25 | 0.203 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 0.892 | 0.663 | 0.036 |
| (86,22) | 44.75 | 4.924 | 9 | 0.000 | 0.000 | 0.000 | 0.000 | 1.660 | 1.225 | 0.070 |

## Local quadratic check (cache neighbourhood, closure ON)

- (70,18): {"C20": {"dep_sigma_obs": 0.00015886926177478822, "jitter_std_sigma_obs": 0.00013642705481027332, "ratio_dep_over_jitter_std": 1.1644996807687806, "pass": true}, "C22": {"dep_sigma_obs": 3.762060436829796e-05, "jitter_std_sigma_obs": 2.466904511297621e-05, "ratio_dep_over_jitter_std": 1.525012589502667, "pass": true}, "C30": {"dep_sigma_obs": 8.664907421374195e-05, "jitter_std_sigma_obs": 7.090729742868357e-05, "ratio_dep_over_jitter_std": 1.2220050313000714, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.00024251707759122737, "jitter_std_sigma_obs": 0.00012800535847266368, "ratio_dep_over_jitter_std": 1.8945853555264904, "pass": true}}
- (77,10): {"C20": {"dep_sigma_obs": 0.0001776758155647313, "jitter_std_sigma_obs": 0.00012266685187930728, "ratio_dep_over_jitter_std": 1.4484419616437838, "pass": true}, "C22": {"dep_sigma_obs": 1.3097718844993252e-05, "jitter_std_sigma_obs": 1.6756589001353133e-05, "ratio_dep_over_jitter_std": 0.7816458853252044, "pass": true}, "C30": {"dep_sigma_obs": 8.952940135536943e-05, "jitter_std_sigma_obs": 6.966915709403482e-05, "ratio_dep_over_jitter_std": 1.2850650860398465, "pass": true}, "libration_deg": {"dep_sigma_obs": 3.061642233718557e-05, "jitter_std_sigma_obs": 0.00011314031253534429, "ratio_dep_over_jitter_std": 0.27060577835703964, "pass": true}}
- (12,6): {"C20": {"dep_sigma_obs": 8.846250471171657e-05, "jitter_std_sigma_obs": 0.00023317862590592982, "ratio_dep_over_jitter_std": 0.37937655892785216, "pass": true}, "C22": {"dep_sigma_obs": 5.113632437996903e-05, "jitter_std_sigma_obs": 9.292730876549174e-05, "ratio_dep_over_jitter_std": 0.5502830659716506, "pass": true}, "C30": {"dep_sigma_obs": 7.228015651138498e-05, "jitter_std_sigma_obs": 6.651542908550851e-05, "ratio_dep_over_jitter_std": 1.0866675221844493, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.005631349182309379, "jitter_std_sigma_obs": 0.00030061058740941786, "ratio_dep_over_jitter_std": 18.7330367530926, "pass": true}}
- (25,30): {"C20": {"dep_sigma_obs": 0.0005170998246965296, "jitter_std_sigma_obs": 0.00039067843207510827, "ratio_dep_over_jitter_std": 1.3235945018769726, "pass": true}, "C22": {"dep_sigma_obs": 0.00020185689422112126, "jitter_std_sigma_obs": 0.00016558063222391517, "ratio_dep_over_jitter_std": 1.219085176267171, "pass": true}, "C30": {"dep_sigma_obs": 0.00014563413271450302, "jitter_std_sigma_obs": 9.888812470843593e-05, "ratio_dep_over_jitter_std": 1.4727160935035843, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.0004790151794757206, "jitter_std_sigma_obs": 0.0004560242417406011, "ratio_dep_over_jitter_std": 1.0504160429001872, "pass": true}}
- (38,18): {"C20": {"dep_sigma_obs": 0.00033955756827790107, "jitter_std_sigma_obs": 0.00021870052306733435, "ratio_dep_over_jitter_std": 1.5526143399911156, "pass": true}, "C22": {"dep_sigma_obs": 0.0001071585975011497, "jitter_std_sigma_obs": 7.0652835484431e-05, "ratio_dep_over_jitter_std": 1.5166921011225811, "pass": true}, "C30": {"dep_sigma_obs": 0.00012416881579824932, "jitter_std_sigma_obs": 7.917791162955915e-05, "ratio_dep_over_jitter_std": 1.56822544624799, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.00030119591072630403, "jitter_std_sigma_obs": 0.00021160385886310963, "ratio_dep_over_jitter_std": 1.4233951703175372, "pass": true}}
- (50,34): {"C20": {"dep_sigma_obs": 0.0004195342899341774, "jitter_std_sigma_obs": 0.00016798737912150047, "ratio_dep_over_jitter_std": 2.4974155328105945, "pass": true}, "C22": {"dep_sigma_obs": 0.0001868826203334644, "jitter_std_sigma_obs": 7.247121340706289e-05, "ratio_dep_over_jitter_std": 2.578715210462465, "pass": true}, "C30": {"dep_sigma_obs": 9.843844536745927e-05, "jitter_std_sigma_obs": 4.0639994469749186e-05, "ratio_dep_over_jitter_std": 2.4222061703461346, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.0006172539011722877, "jitter_std_sigma_obs": 0.000210079270048306, "ratio_dep_over_jitter_std": 2.9381951918928286, "pass": true}}
- (58,12): {"C20": {"dep_sigma_obs": 1.584693172677003e-05, "jitter_std_sigma_obs": 0.00014853578672062467, "ratio_dep_over_jitter_std": 0.10668763451986102, "pass": true}, "C22": {"dep_sigma_obs": 2.923204847707116e-06, "jitter_std_sigma_obs": 3.632261506963381e-05, "ratio_dep_over_jitter_std": 0.08047892042197573, "pass": true}, "C30": {"dep_sigma_obs": 8.254390800169248e-06, "jitter_std_sigma_obs": 6.646929186705207e-05, "ratio_dep_over_jitter_std": 0.12418352247048442, "pass": true}, "libration_deg": {"dep_sigma_obs": 1.8210320270236242e-05, "jitter_std_sigma_obs": 0.0001345867767825032, "ratio_dep_over_jitter_std": 0.1353054193404508, "pass": true}}
- (64,26): {"C20": {"dep_sigma_obs": 0.00019493786984452338, "jitter_std_sigma_obs": 0.000156057848217406, "ratio_dep_over_jitter_std": 1.2491385218444968, "pass": true}, "C22": {"dep_sigma_obs": 6.95430562537399e-05, "jitter_std_sigma_obs": 3.787306102280659e-05, "ratio_dep_over_jitter_std": 1.8362143004987665, "pass": true}, "C30": {"dep_sigma_obs": 8.980010945183964e-05, "jitter_std_sigma_obs": 6.945187774778863e-05, "ratio_dep_over_jitter_std": 1.2929831757457255, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.00041586110314610164, "jitter_std_sigma_obs": 0.000160275732458849, "ratio_dep_over_jitter_std": 2.594660444012474, "pass": true}}
- (84,4): {"C20": {"dep_sigma_obs": 0.00017596927522516847, "jitter_std_sigma_obs": 5.235060541930231e-05, "ratio_dep_over_jitter_std": 3.3613608441725193, "pass": true}, "C22": {"dep_sigma_obs": 1.9722011913910507e-05, "jitter_std_sigma_obs": 5.906108106223704e-06, "ratio_dep_over_jitter_std": 3.339256843796672, "pass": true}, "C30": {"dep_sigma_obs": 0.00010433158616668777, "jitter_std_sigma_obs": 3.0998174197092694e-05, "ratio_dep_over_jitter_std": 3.3657332687830688, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.00016252142467458133, "jitter_std_sigma_obs": 4.833746376457723e-05, "ratio_dep_over_jitter_std": 3.362224908326296, "pass": true}}
- (86,22): {"C20": {"dep_sigma_obs": 0.00017815793749643475, "jitter_std_sigma_obs": 0.00010325566205484832, "ratio_dep_over_jitter_std": 1.7254059869549734, "pass": true}, "C22": {"dep_sigma_obs": 7.823723228996967e-05, "jitter_std_sigma_obs": 1.2061658952995695e-05, "ratio_dep_over_jitter_std": 6.486440430363708, "pass": true}, "C30": {"dep_sigma_obs": 4.626174507238288e-05, "jitter_std_sigma_obs": 6.032737676016178e-05, "ratio_dep_over_jitter_std": 0.7668449642075705, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.00059575173881841, "jitter_std_sigma_obs": 0.00010444802259908893, "ratio_dep_over_jitter_std": 5.7038106035298615, "pass": true}}

TODO (later step): posterior-median check, sigma_model on vs off (<0.2 posterior sigma).
