# Park et al. 2024 (JGR Planets) — ERRATUM, authoritative version of record
# (user-supplied text, filed 2026-08-17; closes the provenance gap flagged
# in r5/r6/r7: papers/park2024global.pdf is the Jan-2024 ORIGINAL)

Campaign-relevant corrections (full text below):
- Libration: 0.091 -> 0.092 deg everywhere; Eq. 6 amplitude 0.092119 deg
  (physical libration; longitude libration 7.120600 deg corrected too).
  == the config's conditioned 0.092 +/- 0.003 and recorded 0.092119. OK.
- Shell thickness 27-33 -> 25-29 km; ocean 21-26 -> 26-30 km; rho_core
  2270-2330 -> 2290-2350 kg/m3. Preferred: shell 27, ocean 28, R_core
  197, rho_core 2320, C/MR2 0.337 (1-sigma 0.335-0.338).
- Best-fit triaxial ellipsoid: 256.21 / 251.39 / 248.63 km, mean 252.08.
- Table 8 note added: C_nm scaling (R_old/R_new)^n — matches
  isostasy.rescale_coeff convention.
- Table 8 corrections affect the SHAPE-DERIVED harmonic columns only;
  the measured gravity-field rows (our C20/C22/C30 source, R_ref 256.6)
  are NOT corrected.
- Table 9 comparison column: Hemingway & Mittal (2019).

NOTES FOR THE CAMPAIGN (manager, 2026-08-17):
1. All conditioned observables VERIFIED against the erratum: libration
   0.092/0.092119 exact; gravity rows untouched. No config change needed.
2. Park PREFERRED C/MR2 = 0.337 (range 0.335-0.338) vs the template
   Cmeasured = 0.335 (Iess): inside the ocean 0.08 window trivially;
   the frozen I-F6 NON-conditioning check tests the template 0.335
   +/- 0.001 — unaffected in mechanics, but record that Park-preferred
   0.337 lies OUTSIDE that band (strengthens I-F6's discriminating
   value). Display references should quote Park 25-29 km shell /
   26-30 km ocean / 2290-2350 rho_core as the comparison bands.
3. The r6-noted "0.33 sigma" original-vs-erratum libration difference is
   now fully documented; the campaign conditions on the erratum values.

---- FULL ERRATUM TEXT (verbatim, user-supplied) ----
"During the further development of our topography model, we identified several inconsistencies that required additional investigation. The comparison of image limbs against synthetics from our model revealed a few problematic areas, with discrepancies exceeding the overall model root-mean-square error. We traced these outliers to a combination of loose constraints during early topography development and the limited ability of the orbital data to provide robust position fixes, allowing heights in certain areas to diverge". In the first key point, 0.091deg corrected to 0.092deg. Third key point: 27-33 km -> 25-29 km; 21-26 km -> 26-30 km; 2,270-2,330 kg/m3 -> 2,290-2,350 kg/m3. [Abstract, Plain Language Summary, Introduction: same corrections.] Table 5 corrected: best-fit tri-axial ellipsoid 256.21, 251.39, 248.63, mean 252.08; best-fit oblate 253.78, -, 248.63, 252.06. Section 2.5: 0.3895->0.63, 0.7018->0.75, 0.3187->0.14; 256.14->256.21, 251.16->251.39, 248.68->248.63; -3.6 km -> -2.7 km; paragraphs 3-4 deleted; Figure 9 removed, subsequent renumbered. Section 2.6: -3.5..+3.5 km -> -2.7..+2.3 km; sixth sentence deleted. Section 3: -0.091 -> -0.092 deg; Table 6: -0.091->-0.092, -0.093->-0.094, -0.092->-0.094. Equation 6: 8.325383 deg -> 7.120600 deg; 0.091295 deg -> 0.092119 deg. Table 8 notes addition: C_nm_new = (R_old/R_new)^n C_nm_old, S likewise. Table 8 "Gravity harmonics from shape assuming constant-density interior" Cnm corrected: -7905.88, 38.46, 1863.77, 894.36, 300.35; Snm: -, 258.86, -359.04, -, -. "Topography harmonics from shape" Cnm: -13425.77, 67.50, 3110.24, 2230.69, 599.79; Snm: -, 453.24, -604.38, -, -. Figure 18 removed. Sections 5.1/5.2: sentences/paragraphs deleted as listed. Section 5.3: 27-33->25-29, 21-26->26-30, 2270-2330->2290-2350; 30 km -> 27 km; 18-28 -> 20-30. Table 9: second column heading -> Hemingway & Mittal (2019); "This work" Preferred: 27, 28, 197, 2,320, 0.337; 1-sigma: 25-29, 26-30, 196-198, 2,290-2,350, 0.335-0.338. Section 6: same libration/range corrections. Data Availability: "(see version 2)" added for the dsk shape models. Supporting information replaced along with Figures 9, 10, 11, 16, 17. This may be considered the authoritative version of record.
