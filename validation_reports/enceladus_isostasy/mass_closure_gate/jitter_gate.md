# Enceladus mass-closure jitter gate

Overall: **FAIL**

builds=90 deltaP/[1, 2, 4] nPoros=[8, 12, 16] mean build 76.4 s, total 1212.8 s; cache=/Users/svance/scratch/enceladus_zbw_production_cache.pkl

## Criteria

- C20_sigma_model<=0.5sigma_obs: False
- C22_sigma_model<=0.5sigma_obs: False
- libration_spread<=0.1sigma_all_nodes: False
- local_quadratic<=3x_jitter_std: True
- posterior_median_shift<0.2_sigma: TODO: later step (sigma_model on vs off posterior runs)

## sigma_model = max(std, |bias|)

{
 "C20": {
  "max_std_over_sigma_obs": 0.6589542713496925,
  "max_abs_bias_over_sigma_obs": 1.272673081763123,
  "sigma_model_over_sigma_obs": 1.272673081763123,
  "sigma_model_abs": 4.707617729441793e-05,
  "recommended_sigma_model_add_abs": 4.737358408074558e-05,
  "existing_sigma_model_add": 5.3e-06
 },
 "C22": {
  "max_std_over_sigma_obs": 0.49113685452171285,
  "max_abs_bias_over_sigma_obs": 0.9457947635486433,
  "sigma_model_over_sigma_obs": 0.9457947635486433,
  "sigma_model_abs": 1.3903183024165057e-05,
  "recommended_sigma_model_add_abs": 1.4006730460868854e-05,
  "existing_sigma_model_add": 1.7e-06
 },
 "C30": {
  "max_std_over_sigma_obs": 0.004021759692977765,
  "max_abs_bias_over_sigma_obs": 0.00894214230107829,
  "sigma_model_over_sigma_obs": 0.00894214230107829,
  "sigma_model_abs": 2.988463957020364e-07,
  "recommended_sigma_model_add_abs": 4.410137091772103e-06,
  "existing_sigma_model_add": 4.4e-06
 },
 "libration_deg": {
  "max_std_over_sigma_obs": 0.07167452817479618,
  "max_abs_bias_over_sigma_obs": 0.10714809614017125,
  "sigma_model_over_sigma_obs": 0.10714809614017125,
  "sigma_model_abs": 0.00032144428842051376,
  "recommended_sigma_model_add_abs": 0.0004072179153207413,
  "existing_sigma_model_add": 0.00025
 }
}

## C20-C22 jitter correlation

{"rho_obs_in_config": 0.47, "pooled": -0.9999992040114171, "per_node": {"(70,18)": -0.9999999982353154, "(77,10)": -0.9999999968048304, "(12,6)": -0.9999999995155819, "(25,30)": -0.9999999979476428, "(38,18)": -0.9999999988385293, "(50,34)": -0.9999999986940443, "(58,12)": -0.999999998334641, "(64,26)": -0.9999999976067588, "(84,4)": -0.999999998555853, "(86,22)": -0.9999999971854545}, "flag_differs_strongly": true, "note": "runner applies rho_obs to the INFLATED sigmas; flag = |pooled - rho_obs| > 0.3 (reporting only, not gating)"}

## Failed builds

{}

## Bias (production - finest), sigma_obs

- (70,18): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.03928073563417198, "C22": -0.029266475789528207, "C30": -0.0002975025437052846, "libration_deg": -0.001932356000211577}
- (77,10): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -0.44892182388572005, "C22": 0.3359707290230956, "C30": 0.0024013790859699156, "libration_deg": 0.01847868573139298}
- (12,6): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.31382606877534586, "C22": -0.2321016803140091, "C30": -0.002580853021207368, "libration_deg": -0.054578864929608985}
- (25,30): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -0.21415075531772448, "C22": 0.15969622818276658, "C30": 0.001277640348139635, "libration_deg": 0.023208327858995997}
- (38,18): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 1.272673081763123, "C22": -0.9457947635486433, "C30": -0.00894214230107829, "libration_deg": -0.10714809614017125}
- (50,34): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -0.6062398510050953, "C22": 0.4546567541832657, "C30": 0.002811144765482189, "libration_deg": 0.048295917834950686}
- (58,12): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.38865807692106524, "C22": -0.28952124815482777, "C30": -0.002558244325280501, "libration_deg": -0.023214436370239122}
- (64,26): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -0.2548403535143451, "C22": 0.19067636359209572, "C30": 0.001349150932380369, "libration_deg": 0.013627623562208665}
- (84,4): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": 0.8247015584144082, "C22": -0.6172415258769515, "C30": -0.004528653253612847, "libration_deg": -0.0321225330393235}
- (86,22): {"finest": {"deltaP": 0.025, "nPoros": 16}, "C20": -0.125245988091446, "C22": 0.09399147719024549, "C30": 0.0005615218451284078, "libration_deg": 0.004757038516027083}

## Frozen two-layer stack (M_stack-M)/M
{"n": 39, "min": -3.03345614782174e-16, "max": 3.03345614782174e-16, "median": 0.0, "min_abs": 0.0, "max_abs": 3.03345614782174e-16}


## Per node: spread / std in sigma_obs (closure ON | OFF)

| node (i,j) | zb | w | n | C20 sp | C22 sp | C30 sp | lib sp | C20 sp OFF | C22 sp OFF | lib sp OFF |
|---|---|---|---|---|---|---|---|---|---|---|
| (70,18) | 39.50 | 2.424 | 9 | 0.685 | 0.512 | 0.004 | 0.032 | 1.229 | 0.891 | 0.068 |
| (77,10) | 42.50 | 0.588 | 9 | 0.641 | 0.480 | 0.004 | 0.027 | 1.256 | 0.902 | 0.065 |
| (12,6) | 17.00 | 0.289 | 9 | 0.791 | 0.585 | 0.006 | 0.138 | 0.817 | 0.596 | 0.157 |
| (25,30) | 23.50 | 20.309 | 9 | 2.227 | 1.660 | 0.014 | 0.242 | 1.547 | 1.164 | 0.154 |
| (38,18) | 26.75 | 2.424 | 9 | 1.407 | 1.045 | 0.010 | 0.118 | 1.415 | 1.052 | 0.119 |
| (50,34) | 29.75 | 41.246 | 9 | 1.238 | 0.928 | 0.006 | 0.099 | 1.701 | 1.262 | 0.149 |
| (58,12) | 33.50 | 0.838 | 9 | 1.371 | 1.022 | 0.009 | 0.082 | 1.638 | 1.210 | 0.103 |
| (64,26) | 36.50 | 10.000 | 9 | 1.286 | 0.962 | 0.007 | 0.069 | 1.704 | 1.258 | 0.100 |
| (84,4) | 44.25 | 0.203 | 9 | 0.825 | 0.617 | 0.005 | 0.032 | 0.892 | 0.663 | 0.036 |
| (86,22) | 44.75 | 4.924 | 9 | 1.611 | 1.207 | 0.008 | 0.063 | 1.660 | 1.225 | 0.070 |

## Local quadratic check (cache neighbourhood, closure ON)

- (70,18): {"C20": {"dep_sigma_obs": 0.5692035298566026, "jitter_std_sigma_obs": 0.3016169240254694, "ratio_dep_over_jitter_std": 1.887173711142739, "pass": true}, "C22": {"dep_sigma_obs": 0.42538244131272324, "jitter_std_sigma_obs": 0.2254593050997305, "ratio_dep_over_jitter_std": 1.8867371259064158, "pass": true}, "C30": {"dep_sigma_obs": 0.00327715683542271, "jitter_std_sigma_obs": 0.001715493153143381, "ratio_dep_over_jitter_std": 1.9103293005966344, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.026534150883393098, "jitter_std_sigma_obs": 0.014053826663503416, "ratio_dep_over_jitter_std": 1.8880374376823656, "pass": true}}
- (77,10): {"C20": {"dep_sigma_obs": 0.32196601156504223, "jitter_std_sigma_obs": 0.2592376378295495, "ratio_dep_over_jitter_std": 1.241972478459077, "pass": true}, "C22": {"dep_sigma_obs": 0.24094587302517043, "jitter_std_sigma_obs": 0.19395609794234012, "ratio_dep_over_jitter_std": 1.2422701610381932, "pass": true}, "C30": {"dep_sigma_obs": 0.0017362752734891782, "jitter_std_sigma_obs": 0.0014230125458270975, "ratio_dep_over_jitter_std": 1.2201405241160421, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.013179853223521665, "jitter_std_sigma_obs": 0.010715636866404579, "ratio_dep_over_jitter_std": 1.2299645264055972, "pass": true}}
- (12,6): {"C20": {"dep_sigma_obs": 0.08735208797452118, "jitter_std_sigma_obs": 0.25709232852377173, "ratio_dep_over_jitter_std": 0.3397693290814948, "pass": true}, "C22": {"dep_sigma_obs": 0.06437520315100186, "jitter_std_sigma_obs": 0.19016168583595078, "ratio_dep_over_jitter_std": 0.3385287781185177, "pass": true}, "C30": {"dep_sigma_obs": 0.0008448625594746763, "jitter_std_sigma_obs": 0.002112129803248499, "ratio_dep_over_jitter_std": 0.4000050367052538, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.016504330071664963, "jitter_std_sigma_obs": 0.04472528913851885, "ratio_dep_over_jitter_std": 0.3690156148694611, "pass": true}}
- (25,30): {"C20": {"dep_sigma_obs": 0.8106775921596484, "jitter_std_sigma_obs": 0.6589542713496925, "ratio_dep_over_jitter_std": 1.23024863394419, "pass": true}, "C22": {"dep_sigma_obs": 0.6044074080438135, "jitter_std_sigma_obs": 0.49113685452171285, "ratio_dep_over_jitter_std": 1.2306293092836775, "pass": true}, "C30": {"dep_sigma_obs": 0.004897678814936822, "jitter_std_sigma_obs": 0.004021759692977765, "ratio_dep_over_jitter_std": 1.2177949924478244, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.08808754780866475, "jitter_std_sigma_obs": 0.07167452817479618, "ratio_dep_over_jitter_std": 1.2289937590358646, "pass": true}}
- (38,18): {"C20": {"dep_sigma_obs": 0.5304187429971371, "jitter_std_sigma_obs": 0.43374617166831964, "ratio_dep_over_jitter_std": 1.2228782122894257, "pass": true}, "C22": {"dep_sigma_obs": 0.3942259567928446, "jitter_std_sigma_obs": 0.32234858757784784, "ratio_dep_over_jitter_std": 1.2229802517674697, "pass": true}, "C30": {"dep_sigma_obs": 0.003693709769821791, "jitter_std_sigma_obs": 0.00304993703085531, "ratio_dep_over_jitter_std": 1.2110773869931157, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.04407292307766671, "jitter_std_sigma_obs": 0.03652319172521349, "ratio_dep_over_jitter_std": 1.2067106130607235, "pass": true}}
- (50,34): {"C20": {"dep_sigma_obs": 0.9671398346762898, "jitter_std_sigma_obs": 0.4241802994018888, "ratio_dep_over_jitter_std": 2.2800206328299444, "pass": true}, "C22": {"dep_sigma_obs": 0.7252579198508873, "jitter_std_sigma_obs": 0.31805274896450697, "ratio_dep_over_jitter_std": 2.28030703149754, "pass": true}, "C30": {"dep_sigma_obs": 0.004501211424672505, "jitter_std_sigma_obs": 0.0019875951402417084, "ratio_dep_over_jitter_std": 2.2646520579262033, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.07762275318390355, "jitter_std_sigma_obs": 0.03387290538434908, "ratio_dep_over_jitter_std": 2.291588285773945, "pass": true}}
- (58,12): {"C20": {"dep_sigma_obs": 0.04874171518521957, "jitter_std_sigma_obs": 0.42739073595770466, "ratio_dep_over_jitter_std": 0.11404485657836798, "pass": true}, "C22": {"dep_sigma_obs": 0.03641805225341721, "jitter_std_sigma_obs": 0.31841682288370954, "ratio_dep_over_jitter_std": 0.11437226187863074, "pass": true}, "C30": {"dep_sigma_obs": 0.00029599383060214255, "jitter_std_sigma_obs": 0.0027689191504195526, "ratio_dep_over_jitter_std": 0.10689869025511017, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.0023609952316862768, "jitter_std_sigma_obs": 0.025464091766750803, "ratio_dep_over_jitter_std": 0.09271861149860826, "pass": true}}
- (64,26): {"C20": {"dep_sigma_obs": 0.7726648550958078, "jitter_std_sigma_obs": 0.4560714104696576, "ratio_dep_over_jitter_std": 1.694175160640141, "pass": true}, "C22": {"dep_sigma_obs": 0.5778691133655, "jitter_std_sigma_obs": 0.3410068468605445, "ratio_dep_over_jitter_std": 1.694596805564496, "pass": true}, "C30": {"dep_sigma_obs": 0.004266256151127648, "jitter_std_sigma_obs": 0.00253519814784029, "ratio_dep_over_jitter_std": 1.6828097459608942, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.04136736886469391, "jitter_std_sigma_obs": 0.024618051741147503, "ratio_dep_over_jitter_std": 1.6803672889983812, "pass": true}}
- (84,4): {"C20": {"dep_sigma_obs": 0.5806310916897595, "jitter_std_sigma_obs": 0.3074486828851689, "ratio_dep_over_jitter_std": 1.8885463624074894, "pass": true}, "C22": {"dep_sigma_obs": 0.4346752335140641, "jitter_std_sigma_obs": 0.23010445620866354, "ratio_dep_over_jitter_std": 1.8890344006197408, "pass": true}, "C30": {"dep_sigma_obs": 0.003124123699442452, "jitter_std_sigma_obs": 0.0016879220416473697, "ratio_dep_over_jitter_std": 1.8508696624361782, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.02249612916897098, "jitter_std_sigma_obs": 0.011975501311597431, "ratio_dep_over_jitter_std": 1.8785125218252918, "pass": true}}
- (86,22): {"C20": {"dep_sigma_obs": 0.6900840328424032, "jitter_std_sigma_obs": 0.4541647728649148, "ratio_dep_over_jitter_std": 1.5194574173802322, "pass": true}, "C22": {"dep_sigma_obs": 0.5170250064185494, "jitter_std_sigma_obs": 0.34038495639540034, "ratio_dep_over_jitter_std": 1.5189419999453773, "pass": true}, "C30": {"dep_sigma_obs": 0.003638427565204752, "jitter_std_sigma_obs": 0.002320235999322068, "ratio_dep_over_jitter_std": 1.568128227588847, "pass": true}, "libration_deg": {"dep_sigma_obs": 0.026488285381429572, "jitter_std_sigma_obs": 0.017645734774641507, "ratio_dep_over_jitter_std": 1.5011154661292767, "pass": true}}

TODO (later step): posterior-median check, sigma_model on vs off (<0.2 posterior sigma).
