# Abstract Draft

**Target Length:** 200-250 words  
**Format:** Single paragraph for AASTeX

---

We implement and validate the Kalousová & Sotin (2018) high-pressure ice convection model in Planet Profile, enabling prediction of partial melt formation in Ice III, V, and VI layers through two-phase convection parameterization. The model evaluates a modified Rayleigh number (Ra*) against a critical threshold (Ra*c) that depends on heat flux from the silicate mantle, determining whether temperate (partially molten) layers form at phase boundaries. Application to Titan under ocean-present and ocean-absent scenarios reveals that ocean presence is critical for high-pressure ice melting: only the ocean configuration achieves supercritical convection (Ra*/Ra*c ~ 11) in Ice VI, forming a 21.8 km temperate layer with ~1% partial melt. The ocean-ice boundary condition reduces Ra*c by a factor ~680 compared to direct ice-silicate contact, demonstrating that phase boundary topology—not just convective vigor—controls melt formation. In contrast, the ocean-free scenario proposed by Petricca et al. (2025) maintains all high-pressure ice layers below the melting threshold despite similar absolute Rayleigh numbers. We find that incorporation of a 4 km clathrate hydrate cap reduces the ice Ih stagnant lid thickness by 93% (48 km → 3.4 km), enabling efficient tidal dissipation in the upper ice shell consistent with Titan's observed heat flux without requiring a subsurface ocean. Both scenarios require fluid-bearing high-pressure ice viscosities (10¹⁰–10¹² Pa s), two orders of magnitude below laboratory values for pure crystalline ice, consistent with small amounts of interstitial fluid enhancing creep rates. Our results support the paradigm shift toward ocean-free Titan interiors while demonstrating that high-pressure ice partial melting remains viable for ocean-bearing configurations. The implementation provides a framework for evaluating temperate layer formation across diverse ocean world scenarios, with applications to Europa, Enceladus, Ganymede, and exoplanetary ice-rich bodies.

---

**Word Count:** 281 words (needs trimming to 250)

**Key Points Covered:**
✓ Model implementation and validation  
✓ Application to Titan (two scenarios)  
✓ Key finding: Ocean presence critical (680× Ra*c factor)  
✓ Clathrate effect quantified (93% reduction)  
✓ Viscosity requirements (fluid-bearing)  
✓ Connection to Petricca et al. (2025) paradigm  
✓ Broader implications (other ocean worlds)

**Trimming Strategy:**
- Reduce methodology description (first 2 sentences → 1 sentence)
- Condense clathrate result (one number instead of two)
- Shorten final applications sentence

---

## Revised Abstract (247 words)

We implement the Kalousová & Sotin (2018) two-phase convection model in PlanetProfile, enabling prediction of partial melt formation in high-pressure ice (Ice III, V, VI) through parameterized boundary layer theory that evaluates convective vigor (modified Rayleigh number Ra*) against a critical threshold (Ra*c) for temperate layer formation. Application to Titan reveals that ocean presence is critical for high-pressure ice melting: the ocean-present scenario achieves supercritical convection in Ice VI (Ra*/Ra*c ~ 11), forming a 21.8 km temperate layer with ~1% partial melt, while the ocean-free configuration maintains all high-pressure ice layers subcritical despite similar absolute Rayleigh numbers. The ocean-ice boundary condition reduces Ra*c by a factor ~680 compared to direct ice-silicate contact, demonstrating that phase boundary topology—not just convective intensity—controls melt formation. For the ocean-free scenario proposed by Petricca et al. (2025), we find that a 4 km clathrate hydrate cap reduces ice Ih stagnant lid thickness by 93%, enabling efficient tidal dissipation without requiring subsurface ocean heat transport. Both scenarios require fluid-bearing high-pressure ice viscosities (10¹⁰–10¹² Pa s), two orders of magnitude below pure ice laboratory values, consistent with small amounts of interstitial fluid enhancing grain boundary creep. Our results support the paradigm shift toward ocean-free Titan interiors while demonstrating that high-pressure ice partial melting remains viable for ocean-bearing moons. The framework enables evaluation of temperate layer formation across ocean worlds including Europa, Enceladus, Ganymede, and exoplanetary ice-rich bodies.

---

**Status:** Ready for LaTeX integration ✅
