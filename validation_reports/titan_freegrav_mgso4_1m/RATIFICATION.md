# Titan free-gravity MgSO4 JOINT artifact — ratification

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
  or eta posteriors; use the reference MCMC. [S1] [S2]

Scientific-reviewer gate-set interpretation: agent `a719438103ced11d0`,
2026-08-10, PASS WITH CONCERNS and no STOP. Re_k2 pushforward adjudication:
agent `ad7ae4436e53cdc1`, 2026-08-10, PASS for split-status deployment. [S2]

The underlying campaign was independently checked and manager-countersigned on
2026-08-10. This copy-only consolidation awaits the manager's post-commit
countersign. [S1] [S3]

## Artifact record

Artifact: `PlanetProfile/Inference/sbi_artifacts/`
`titan_freegrav_mgso4_posterior_1m.pt`; config:
`PlanetProfile/Inference/configs/test54_titan_mgso4_freegrav.json`; SBI config
hash `d1bfdf1ee55bc737`; training seed 73; 691,075 of 1,000,000 simulations
kept. [S1] [S4] [S5]

The artifact is the 13-parameter JOINT no-ocean+ocean MgSO4 model. Frozen nodes
are retained as real no-ocean structures. [S1]

## Gate outcomes

| Gate | Raw verdict | Recorded reading | Source |
|---|---|---|---|
| SBC | **PASS** | Requested `n_sbc=1500`; 1,024 kept pairs; all 13 parameters PASS after BH-FDR (minimum adjusted p = 0.156). | [S6] [S2] |
| crosscheck | **PASS** | 13/13 parameters PASS against the repaired pooled B3 reference; `log10_eta_V` median difference is 0.209 dex, below the 0.30-dex tolerance. | [S7] [S3] |
| limits | **FAIL** | Containment PASS; maximum median shift 0.00406 sigma. The decreasing-monotonicity clause is a gated FAIL (`monotone_pass=false`, `monotone_gated=true`, verdict FAIL), not N/A. | [S8] [S3] |
| observable pushforward | **SPLIT** | Reference-MCMC C20/C22 medians center within 0.05 sigma of the data. Re_k2 median is 0.570 vs 0.608 observed; Im_k2 is 1.60 sigma low, so tidal quarantine stands. | [S9] [S2] |

## Adjudications and required caveats

The limits FAIL is retained as a FAIL and overridden for deployment because the
monotone-decreasing premise is falsified for the Titan joint mixture. The real
datum, Im_k2 = 0.135, is below the 0.15 falsified-window boundary; containment
is the binding read. The gate was not tuned or relabeled PASS. [S3] [S8]

The MCMC Re_k2 pushforward median is 0.570 and the SBI-vs-MCMC median-to-median
gap is 0.24 measurement sigma, within the 0.5-sigma rule. Both samplers
under-predict in the same direction, establishing a shared model-data tension
rather than a flow-only offset. [S2] [S9]

The high-salinity MgSO4 tail (approximately w >= 160 ppt) uses linearly
extrapolated EOS values above 800 MPa; any high-w claim must carry that caveat.
[S1]

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
