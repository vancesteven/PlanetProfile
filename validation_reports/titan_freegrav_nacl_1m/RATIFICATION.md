# Titan free-gravity NaCl JOINT artifact — ratification

**Campaign verdict (2026-08-10): SPLIT ratification.** The gravity/structure
sector is trusted; Re_k2 is informative with a model-data-tension caveat; and
the Im_k2/dissipation sector remains quarantined, with the reference MCMC
authoritative there. [S1] [S2]

- **VERIFIED for gravity/structure:** layer thicknesses, salinity, Tb, core
  properties, C20/C22, and dC20_nh/dC22_nh. [S1] [S2]
- **INFORMATIVE WITH CAVEAT for Re_k2:** both samplers under-predict the datum;
  the deployed SBI marginal is conservative relative to the reference MCMC
  (tension-leaning, not a tightened bound). [S1] [S2]
- **NOT VERIFIED for Im_k2/dissipation:** do not quote SBI Im_k2, log10_zeta,
  or eta posteriors; use the reference MCMC. The `log10_eta_V` nuisance
  marginal is separately do-not-cite because its crosscheck FAIL is
  `FAIL-ADJUDICATED-ACCEPTABLE`. [S1] [S3]

Scientific-reviewer gate-set interpretation: agent `a719438103ced11d0`,
2026-08-10, PASS WITH CONCERNS and no STOP. Re_k2 pushforward adjudication:
agent `ad7ae4436e53cdc1`, 2026-08-10, PASS for split-status deployment. [S2]

The underlying campaign was independently checked and manager-countersigned on
2026-08-10. This copy-only consolidation awaits the manager's post-commit
countersign. [S1] [S3]

## Artifact record

Artifact: `PlanetProfile/Inference/sbi_artifacts/`
`titan_freegrav_nacl_posterior_1m.pt`; config:
`PlanetProfile/Inference/configs/test54_titan_nacl_freegrav.json`; SBI config
hash `5b3d2502546548c3`; training seed 74; 631,214 of 1,000,000 simulations
kept. [S1] [S4] [S5]

The artifact is the 13-parameter JOINT no-ocean+ocean NaCl model. Frozen nodes
are retained as real no-ocean structures. [S1]

## Gate outcomes

| Gate | Raw verdict | Recorded reading | Source |
|---|---|---|---|
| SBC | **PASS** | Requested `n_sbc=1500`; 971 kept pairs; all 13 parameters PASS after BH-FDR (minimum adjusted p = 0.881). | [S6] [S2] |
| crosscheck | **FAIL** | 12/13 parameters PASS. `log10_eta_V` alone fails shape+median: median difference 0.352 dex vs 0.30-dex tolerance; all observable-relevant parameters and primary `log10_eta_Ih` PASS. | [S7] [S3] |
| limits | **FAIL** | Containment PASS; maximum median shift 0.00722 sigma. The decreasing-monotonicity clause is a gated FAIL (`monotone_pass=false`, `monotone_gated=true`, verdict FAIL), not N/A. | [S8] [S3] |
| observable pushforward | **SPLIT** | Reference-MCMC C20/C22 medians center within 0.05 sigma of the data. Re_k2 median is 0.575 vs 0.608 observed; Im_k2 is 1.34 sigma low, so tidal quarantine stands. | [S9] [S2] |

## Adjudications and required caveats

The `log10_eta_V` crosscheck result remains
`FAIL-ADJUDICATED-ACCEPTABLE`: it is a poorly identified high-pressure ice-V
nuisance, while every observable-relevant parameter passes. The marginal is
do-not-cite, and the raw FAIL is not relabeled PASS. [S1] [S3] [S7]

The limits FAIL is retained as a FAIL and overridden for deployment because the
monotone-decreasing premise is falsified for the Titan joint mixture. The real
datum, Im_k2 = 0.135, is below the 0.15 falsified-window boundary; containment
is the binding read. The gate was not tuned or relabeled PASS. [S3] [S8]

The MCMC Re_k2 pushforward median is 0.575 and the SBI-vs-MCMC median-to-median
gap is 0.53 measurement sigma. The reviewer adjudicated the slight threshold
breach as mechanical because both samplers under-predict in the same direction;
the shared model-data-tension reading survives. [S2] [S9]

## Sources

[S1]: ../../PlanetProfile/Inference/sbi_artifacts/INDEX.md
[S2]: ../../plans/MACHINE-B-HANDOFF.md
[S3]: ../../plans/STATUS.md
[S4]: gen_manifest.json
[S5]: sbc/sbc_report.json
[S6]: sbc/sbc_report.json
[S7]: crosscheck/crosscheck_report.json
[S8]: limits/limits_report.json
[S9]: rek2_pushforward/mcmc_rek2_pushforward_report.json

## Manager countersign

Countersigned 2026-08-12 (Machine A manager, Fable 5) after file-level
verification: every gate-table number re-checked against the committed
per-gate JSON reports (SBC kept-pair counts and minimum BH-adjusted p
recomputed from sbc_report.json; crosscheck, limits-containment, and
pushforward figures against their reports) and against the 2026-08-10
INDEX/STATUS adjudication records. Verdict text matches the ratified
split-status scope, including the FAIL-ADJUDICATED vocabulary. This file
is the consolidated citation target for the campaign (methods
extractability, doc-doctor item E11).
