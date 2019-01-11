This is a read me file to describe models "CV_hhph_DEW17_fluid_nomelt_685", "CM_hhph_DEW17_fluid_nomelt_685", "CI_hhph_DEW17_fluid_nomelt_685", and "H_hhph_DEW17_fluid_nomelt_685"

These models represent generic undifferentiated bodies of CV, M, CI and H chondrite composition respectively, including volatile elements (composition simplified from Lodders and Fegley, 1998, The Planetary Scientist's Companion book):

Bulk CV chondrite composition:
Element	Weight%
H2	0.29
C	0.55
Mg	14.74
Al	1.73
Si	16.18
S2	2.27
Ca	1.9
Fe	24.22
O2	38.13
SUM	100.01

Bulk CM chondrite composition
Element	Weight%
H	1.44
C	2.26
Mg	11.8
Al	1.16
Si	13.04
S	2.77
Ca	1.32
Fe	21.86
O	44.34
SUM	99.99

Bulk CI chondrite composition:
H	2.07
C	3.53
Mg	9.94
Al	0.89
Si	10.9
S	5.54
Ca	0.95
Fe	18.65
O	47.54
SUM	100.01

Bulk H chondrite composition:
H	0.02
C	0.21
Mg	14.3
Al	1.07
Si	17.34
S	2.03
Ca	1.24
Fe	27.58
O	36.2
SUM	99.99

The thermodynamic database used was DEW17HP622ver_elements.dat.

Solution models included, and phases excluded, were based on the klb-1_hhph and morb_hhph benchmarks for Holland et al. (2013, J. Metamorph. Petrol.)'s thermodynamic set, implemented in Perple_X 6.8.5. Namely, no melt phases were allowed, among others.

Excluded phases:
stv, perL, limL, corL, hemL, qL, h2oL, foL, faL, woL, enL, diL, silL, anL, sil8L, fo8L, fa8L, q8L, cess, qjL, dijL, ctjL, fojL, fajL, hmjL, foTL, faTL.

Also excluded for H chondrite composition (because Holland et al., 2018 J. Pet. definitions were included in Perple_X 6.8.5 when the H chondrite model was calculated):
foHL, faHL, qHL

Solutions included:
Gt(H), Fper(H), Mpv(H), Cpv(H), CFer(H), Aki(H), O(JH), Wad(H), Ring(H), Cpx(JH), Opx(JH), Hpx(H), NAl(H), Cor(H).

Note that fluids (not melt) are allowed and that WERAMI tables show bulk (aggregate) properties of fluids and solids, unless specified.