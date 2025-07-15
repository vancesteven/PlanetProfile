import numpy as np
import logging
from copy import deepcopy
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
from scipy.optimize import root_scalar as GetZero
from seafreeze.seafreeze import seafreeze as SeaFreeze
from seafreeze.seafreeze import whichphase as WhichPhase
from PlanetProfile.Thermodynamics.Clathrates.ClathrateProps import ClathProps, ClathStableSloan1998, \
    ClathStableNagashima2017, ClathSeismic
from PlanetProfile.Utilities.DataManip import ResetNearestExtrap, ReturnZeros, EOSwrapper, ReturnConstantSpecies
from PlanetProfile.Thermodynamics.InnerEOS import GetphiFunc, GetphiCalc
from PlanetProfile.Thermodynamics.MgSO4.MgSO4Props import MgSO4Props, MgSO4PhaseMargules, MgSO4PhaseLookup, \
    MgSO4Seismic, MgSO4Conduct, Ppt2molal
from PlanetProfile.Thermodynamics.Seawater.SwProps import SwProps, SwPhase, SwSeismic, SwConduct
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
from PlanetProfile.Utilities.Indexing import PhaseConv, PhaseInv, MixedPhaseSeparator
from PlanetProfile.Thermodynamics.Reaktoro.reaktoroProps import RktPhaseLookup, RktPhaseOnDemand,  \
    SpeciesParser, RktProps, RktSeismic, RktConduct, RktHydroSpecies, RktRxnAffinity, EOSLookupTableLoader
# Assign logger
log = logging.getLogger('PlanetProfile')

def GetOceanEOS(compstr, wOcean_ppt, P_MPa, T_K, elecType, rhoType=None, scalingType=None, phaseType=None,
                EXTRAP=False, FORCE_NEW=False, MELT=False, PORE=False, sigmaFixed_Sm=None, LOOKUP_HIRES=False,
                etaFixed_Pas=None):
    oceanEOS = OceanEOSStruct(compstr, wOcean_ppt, P_MPa, T_K, elecType, rhoType=rhoType, scalingType=scalingType,
                              phaseType=phaseType, EXTRAP=EXTRAP, FORCE_NEW=FORCE_NEW, MELT=MELT, PORE=PORE,
                              sigmaFixed_Sm=sigmaFixed_Sm, LOOKUP_HIRES=LOOKUP_HIRES, etaFixed_Pas=etaFixed_Pas)
    if oceanEOS.ALREADY_LOADED and not FORCE_NEW:
        log.debug(f'{wOcean_ppt} ppt {compstr} EOS already loaded. Reusing existing EOS.')
        oceanEOS = EOSlist.loaded[oceanEOS.EOSlabel]

    # Ensure each EOSlabel is included in EOSlist, in case we have reused EOSs with
    # e.g. a smaller range that can reuse the larger-range already-loaded EOS.
    if oceanEOS.EOSlabel not in EOSlist.loaded.keys() or FORCE_NEW:
        EOSlist.loaded[oceanEOS.EOSlabel] = oceanEOS

    oceanEOSwrapper = EOSwrapper(oceanEOS.EOSlabel)

    return oceanEOSwrapper

class OceanEOSStruct:
    def __init__(self, compstr, wOcean_ppt, P_MPa, T_K, elecType, rhoType=None, scalingType=None,
                 phaseType=None, EXTRAP=False, FORCE_NEW=False, MELT=False, PORE=False,
                 sigmaFixed_Sm=None, LOOKUP_HIRES=False, etaFixed_Pas=None):
        if elecType is None:
            self.elecType = 'Vance2018'
        else:
            self.elecType = elecType
        if rhoType is None:
            self.rhoType = 'Millero'
        else:
            self.rhoType = rhoType
        if scalingType is None:
            self.scalingType = 'Vance2018'
        else:
            self.scalingType = scalingType
        if phaseType is None or phaseType == 'lookup':
             self.PHASE_LOOKUP = True
        else:
            self.PHASE_LOOKUP = False
        # Add ID for melting curve EOS
        if MELT:
            meltStr = f'melt{np.max(P_MPa)}'
            meltPrint = 'melting curve '
        else:
            meltStr = ''
            meltPrint = ''

        self.EOSlabel = f'{meltStr}Comp{compstr}wppt{wOcean_ppt}elec{elecType}rho{rhoType}' + \
                        f'scaling{scalingType}phase{phaseType}extrap{EXTRAP}pore{PORE}' + \
                        f'hires{LOOKUP_HIRES}etaFixed{etaFixed_Pas}'
        self.ALREADY_LOADED, self.rangeLabel, P_MPa, T_K, self.deltaP, self.deltaT \
            = CheckIfEOSLoaded(self.EOSlabel, P_MPa, T_K, FORCE_NEW=FORCE_NEW)

        if not self.ALREADY_LOADED or FORCE_NEW:
            self.comp = compstr
            self.w_ppt = wOcean_ppt
            self.EXTRAP = EXTRAP
            self.EOStype = 'ocean'

            self.Pmin = np.min(P_MPa)
            self.Pmax = np.max(P_MPa)
            self.Tmin = np.min(T_K)
            self.Tmax = np.max(T_K)
            if self.w_ppt is None:
                wStr = '0.0'
            else:
                wStr = f'{self.w_ppt:.1f}'
            log.debug(f'Loading {meltPrint}EOS for {wStr} ppt {self.comp} with ' +
                      f'P_MPa = [{self.Pmin:.1f}, {self.Pmax:.1f}, {self.deltaP:.3f}], ' +
                      f'T_K = [{self.Tmin:.1f}, {self.Tmax:.1f}, {self.deltaT:.3f}], ' +
                      f'for [min, max, step] with EXTRAP = {self.EXTRAP}.')

            # Get tabular data from the appropriate source for the specified ocean composition
            if self.comp == 'none':
                self.ufn_phase = ReturnZeros(1)
                self.type = 'No H2O'
                self.m_gmol = np.nan
                rho_kgm3 = np.zeros((np.size(P_MPa), np.size(T_K)))
                Cp_JkgK = rho_kgm3
                alpha_pK = rho_kgm3
                kTherm_WmK = rho_kgm3
                self.ufn_Seismic = ReturnZeros(2)
                self.ufn_sigma_Sm = ReturnZeros(1)
                self.EOSdeltaP = None
                self.EOSdeltaT = None
                self.propsPmax = 0
            elif self.comp in ['PureH2O', 'NH3', 'NaCl']:
                self.type = 'SeaFreeze'
                self.m_gmol = Constants.m_gmol[self.comp]

                # Set extrapolation boundaries to limits defined in SeaFreeze
                Pmax = {'PureH2O':   2300.6, 'NH3': 2228.4, 'NaCl': 1000.1}
                #Pmax = {'PureH2O':   2300.6, 'NH3': 2228.4, 'NaCl': 5000.1}
                Tmin = {'PureH2O':    239,   'NH3':  241,   'NaCl':  229.0}
                Tmax = {'PureH2O':    501,   'NH3':  399.2, 'NaCl':  501.0}
                wMax = {'PureH2O': np.nan,   'NH3':  290.1, 'NaCl':  290.3} # upper NaCl concentration is 7mol/kgH2O
                self.Pmax = np.minimum(self.Pmax, Pmax[self.comp])
                self.Tmin = np.maximum(self.Tmin, Tmin[self.comp])
                self.Tmax = np.minimum(self.Tmax, Tmax[self.comp])
                self.propsPmax = self.Pmax
                if np.size(P_MPa) == np.size(T_K):
                    log.warning(f'Both P and T inputs have length {np.size(P_MPa)}, but they are organized to be ' +
                                 'used as a grid. This will cause an error in SeaFreeze. P list will be adjusted slightly.')
                    P_MPa = np.linspace(P_MPa[0], P_MPa[-1], np.size(P_MPa)+1)
                if np.max(P_MPa) > self.Pmax:
                    log.warning(f'Input Pmax greater than SeaFreeze limit for {self.comp}. Resetting to SF max of {self.Pmax} MPa.')
                    P_MPa = np.linspace(np.min(P_MPa), self.Pmax, np.size(P_MPa))
                if np.min(T_K) < self.Tmin:
                    log.warning(f'Input Tmin less than SeaFreeze limit for {self.comp}. Resetting to SF min of {self.Tmin} K.')
                    T_K = np.linspace(self.Tmin, np.max(T_K), np.size(T_K))
                if np.max(T_K) > self.Tmax:
                    log.warning(f'Input Tmax greater than SeaFreeze limit for {self.comp}. Resetting to SF max of {self.Tmax} K.')
                    T_K = np.linspace(np.min(T_K), self.Tmax, np.size(T_K))

                if self.comp == 'PureH2O':
                    SFcomp = 'water1'
                    PTmGrid = sfPTgrid(P_MPa, T_K)
                    self.ufn_sigma_Sm = H2Osigma_Sm(sigmaFixed_Sm)
                else:
                    if self.w_ppt > wMax[self.comp]:
                        log.warning(f'Input wOcean_ppt greater than SeaFreeze limit for {self.comp}. Resetting to SF max.')
                        self.w_ppt = wMax[self.comp]
                    if self.comp == 'NaCl':
                        SFcomp = 'NaClaq'
                        self.ufn_sigma_Sm = H2Osigma_Sm(sigmaFixed_Sm)
                    else:
                        SFcomp = self.comp
                        self.ufn_sigma_Sm = H2Osigma_Sm(sigmaFixed_Sm)  # Placeholder until lab data can be implemented
                    if self.w_ppt > wMax[self.comp]:
                        log.warning(f'Input wOcean_ppt greater than SeaFreeze limit for {self.comp}. Resetting to SF max.')
                        self.w_ppt = wMax[self.comp]
                    PTmGrid = sfPTmGrid(P_MPa, T_K, Ppt2molal(self.w_ppt, self.m_gmol))
                seaOut = SeaFreeze(deepcopy(PTmGrid), SFcomp)
                rho_kgm3 = seaOut.rho
                Cp_JkgK = seaOut.Cp
                alpha_pK = seaOut.alpha
                kTherm_WmK = np.zeros_like(alpha_pK) + Constants.kThermWater_WmK  # Placeholder until we implement a self-consistent calculation

                if self.PHASE_LOOKUP:
                    if self.comp == 'PureH2O':
                        self.phase = WhichPhase(deepcopy(PTmGrid))  # FOR COMPATIBILITY WITH SF v0.9.2: Use default comp of water1 here. This is not robust, but allows support for in-development updates to SeaFreeze.
                    else:
                        self.phase = WhichPhase(deepcopy(PTmGrid), solute=SFcomp)
                    # Create phase finder -- note that the results from this function must be cast to int after retrieval
                    self.ufn_phase = RGIwrap(RegularGridInterpolator((P_MPa, T_K), self.phase, method='nearest'),
                                             self.deltaP, self.deltaT)
                    # Save EOS grid resolution in lookup table
                    self.EOSdeltaP = self.deltaP
                    self.EOSdeltaT = self.deltaT
                else:
                    self.ufn_phase = SFphase(self.w_ppt, self.comp)
                    # Lookup table is not used -- flag with nan for grid resolution.
                    self.EOSdeltaP = np.nan
                    self.EOSdeltaT = np.nan

                self.ufn_Seismic = SFSeismic(self.comp, P_MPa, T_K, seaOut, self.w_ppt, self.EXTRAP)
                if self.comp == 'PureH2O':
                    Ocean_Speciation_Info = Constants.KnownCompositions['PureH2O']
                elif self.comp == 'NaCl':
                    Ocean_Speciation_Info = Constants.KnownCompositions['NaCl']
                else:
                    # Placeholder until we get other species
                    Ocean_Speciation_Info = {'ppt_reference_g_kg': None, 'pH': None, 'species': None}
                self.ufn_species = ReturnConstantSpecies(wOcean_ppt, Ocean_Speciation_Info['ppt_reference_g_kg'],
                    Ocean_Speciation_Info['pH'], Ocean_Speciation_Info['species'])

            elif self.comp == 'Seawater':
                self.type = 'GSW'
                self.m_gmol = Constants.m_gmol['H2O']
                if((self.Tmin <= 250) or (self.Pmax > Constants.PminHPices_MPa)):
                    log.warning('GSW handles only ice Ih for determining phases in the ocean. At ' +
                                'low temperatures or high pressures, this model will be wrong as no ' +
                                'high-pressure ice phases will be found.')
                    self.Pmax = Constants.PminHPices_MPa
                if self.Tmax > 350:
                    log.warning('GSW yields physically valid properties only up to about 350 K. ' +
                                'Maximum temperature for this Seawater EOS will be set to that value.')
                    self.Tmax = 350

                self.ufn_phase = SwPhase(self.w_ppt)
                # Lookup table is not used -- flag with nan for grid resolution.
                self.EOSdeltaP = np.nan
                self.EOSdeltaT = np.nan
                rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK = SwProps(P_MPa, T_K, self.w_ppt)
                self.ufn_Seismic = SwSeismic(self.w_ppt, self.EXTRAP)
                if sigmaFixed_Sm is not None:
                    self.ufn_sigma_Sm = H2Osigma_Sm(sigmaFixed_Sm)
                else:
                    self.ufn_sigma_Sm = SwConduct(self.w_ppt)
                self.propsPmax = self.Pmax
                Ocean_Speciation_Info = Constants.KnownCompositions['Seawater']
                self.ufn_species = ReturnConstantSpecies(wOcean_ppt, Ocean_Speciation_Info['ppt_reference_g_kg'],
                    Ocean_Speciation_Info['pH'], Ocean_Speciation_Info['species'])
            elif self.comp == 'MgSO4':
                if self.elecType == 'Pan2020' and round(self.w_ppt) != 100:
                    log.warning('elecType "Pan2020" behavior is defined only for Ocean.wOcean_ppt = 100. ' +
                                'Defaulting to elecType "Vance2018".')
                    self.elecType = 'Vance2018'
                self.type = 'ChoukronGrasset2010'
                self.m_gmol = Constants.m_gmol['MgSO4']
                P_MPa, T_K, rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK \
                    = MgSO4Props(P_MPa, T_K, self.w_ppt, self.EXTRAP)
                if self.PHASE_LOOKUP:
                    self.ufn_phase = MgSO4PhaseLookup(self.w_ppt, HIRES=LOOKUP_HIRES)
                    self.phasePmax = self.ufn_phase.Pmax
                    # Save EOS grid resolution from MgSO4 lookup table loaded from disk
                    self.EOSdeltaP = self.ufn_phase.deltaP
                    self.EOSdeltaT = self.ufn_phase.deltaT
                else:
                    Margules = MgSO4PhaseMargules(self.w_ppt)
                    self.ufn_phase = Margules.arrays
                    self.phasePmax = Margules.Pmax
                    # Lookup table is not used -- flag with nan for grid resolution.
                    self.EOSdeltaP = np.nan
                    self.EOSdeltaT = np.nan
                # self.ufn_species =
                self.ufn_Seismic = MgSO4Seismic(self.w_ppt, self.EXTRAP)
                if sigmaFixed_Sm is not None:
                    self.ufn_sigma_Sm = H2Osigma_Sm(sigmaFixed_Sm)
                else:
                    self.ufn_sigma_Sm = MgSO4Conduct(self.w_ppt, self.elecType, rhoType=self.rhoType,
                                                    scalingType=self.scalingType)
                self.propsPmax = self.ufn_Seismic.Pmax
                self.Pmax = np.min([self.Pmax, self.phasePmax])
                Ocean_Speciation_Info = Constants.KnownCompositions['MgSO4']
                self.ufn_species = ReturnConstantSpecies(wOcean_ppt, Ocean_Speciation_Info['ppt_reference_g_kg'],
                    Ocean_Speciation_Info['pH'], Ocean_Speciation_Info['species'])
            elif self.comp.startswith("CustomSolution"):
                # Parse out the species list and ratio into a format compatible with Reaktoro and create a CustomSolution EOS label
                self.aqueous_species_string, self.speciation_ratio_mol_kg, self.ocean_solid_phases, self.EOS_lookup_label = SpeciesParser(self.comp, self.w_ppt)

                EOSLookupTable = EOSLookupTableLoader(self.aqueous_species_string, self.speciation_ratio_mol_kg, self.ocean_solid_phases, self.EOS_lookup_label)
                self.type = 'Reaktoro'
                P_MPa, T_K, rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK, self.EOSdeltaP, self.EOSdeltaT = (
                    RktProps(EOSLookupTable, P_MPa, T_K, self.EXTRAP))
                self.ufn_Seismic = RktSeismic(EOSLookupTable, self.EXTRAP)


                self.ufn_phase = RktPhaseLookup(EOSLookupTable, P_MPa, T_K, self.deltaP, self.deltaT)

                self.ufn_species = RktHydroSpecies(self.aqueous_species_string, self.speciation_ratio_mol_kg, self.ocean_solid_phases)
                self.ufn_rxn_affinity = RktRxnAffinity(self.aqueous_species_string, self.speciation_ratio_mol_kg, self.ocean_solid_phases)
                if sigmaFixed_Sm is not None:
                    self.ufn_sigma_Sm = H2Osigma_Sm(sigmaFixed_Sm)
                else:
                    # ions = {'Na_p1': {'mols': 0.1}, 'Cl_m1': {'mols': 0.1}}
                    # self.ufn_sigma_Sm = elecCondMcCleskey2012(T_K,ions) # see McCleskeyFig1 benchmark for example usage. this is a placeholder that doesn't have the inputs set up correctly. Has no pressure dependence currently
                    self.ufn_sigma_Sm = RktConduct(self.aqueous_species_string, self.speciation_ratio_mol_kg, self.ocean_solid_phases, self.ufn_species)
                self.propsPmax = self.Pmax
            else:
                raise ValueError(f'Unable to load ocean EOS. self.comp="{self.comp}" but options are "Seawater", "NH3", "MgSO4", ' +
                                 '"NaCl", "CustomSolution", and "none" (for waterless bodies).')

            self.ufn_rho_kgm3 = RectBivariateSpline(P_MPa, T_K, rho_kgm3)
            self.ufn_Cp_JkgK = RectBivariateSpline(P_MPa, T_K, Cp_JkgK)
            self.ufn_alpha_pK = RectBivariateSpline(P_MPa, T_K, alpha_pK)
            self.ufn_kTherm_WmK = RectBivariateSpline(P_MPa, T_K, kTherm_WmK)
            self.ufn_eta_Pas = ViscOceanUniform_Pas(etaSet_Pas=etaFixed_Pas, comp=compstr)

            # Include placeholder to overlap infrastructure with other EOS classes
            self.fn_porosCorrect = None

            # Store complete EOSStruct in global list of loaded EOSs,
            # but only if we weren't forcing a recalculation. This allows
            # us to use a finer step in getting the ice shell thickness while
            # not slowing down ocean calculations.
            if not FORCE_NEW:
                EOSlist.loaded[self.EOSlabel] = self
                EOSlist.ranges[self.EOSlabel] = self.rangeLabel

    # Limit extrapolation to use nearest value from evaluated fit
    def fn_phase(self, P_MPa, T_K, grid=False):
        # Phase stability cannot be extrapolated for some compositions. Therefore, prevent it.
        P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_phase(P_MPa, T_K, grid=grid)
    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_rho_kgm3(P_MPa, T_K, grid=grid)
    def fn_Cp_JkgK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_Cp_JkgK(P_MPa, T_K, grid=grid)
    def fn_alpha_pK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_alpha_pK(P_MPa, T_K, grid=grid)
    def fn_kTherm_WmK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_kTherm_WmK(P_MPa, T_K, grid=grid)
    def fn_Seismic(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_Seismic(P_MPa, T_K, grid=grid)
    def fn_sigma_Sm(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_sigma_Sm(P_MPa, T_K, grid=grid)
    def fn_eta_Pas(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_eta_Pas(P_MPa, T_K, grid=grid)
    def fn_species(self, P_MPa, T_K, grid = False):
        """
        Returns speciation at provided P_MPa and T_K
        """
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_species(P_MPa, T_K, grid=grid)
    def fn_rxn_affinity(self, P_MPa, T_K, reaction, concentrations, grid = False):
        """
        Calculates the affinity of a reaction whose species are at prescribed concentrations at disequilibrium

        ONLY APPLICABLE TO CUSTOMSOLUTION FOR NOW
        """
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_rxn_affinity(P_MPa, T_K, reaction, concentrations, grid=grid)


def GetIceEOS(P_MPa, T_K, phaseStr, porosType=None, phiTop_frac=0, Pclosure_MPa=0, phiMin_frac=0,
              EXTRAP=False, ClathDissoc=None, minPres_MPa=None, minTres_K=None,
            ICEIh_DIFFERENT=False, etaFixed_Pas=None, TviscTrans_K=None, mixParameters=None):
    # Check if this is a mixed EOS
    try:
        # Check if MixedPhaseSeparator can process this phaseStr (will raise error if not mixed)
        phaseOne, phaseTwo = MixedPhaseSeparator(phaseStr)
    except Exception:
        # MixedPhaseSeparator raised an error, so this is not a mixed EOS - create normal ice EOS
        iceEOS = IceEOSStruct(P_MPa, T_K, phaseStr, porosType=porosType, phiTop_frac=phiTop_frac,
                            Pclosure_MPa=Pclosure_MPa, phiMin_frac=phiMin_frac, EXTRAP=EXTRAP,
                            ClathDissoc=ClathDissoc, minPres_MPa=minPres_MPa, minTres_K=minTres_K,
                            ICEIh_DIFFERENT=ICEIh_DIFFERENT, etaFixed_Pas=etaFixed_Pas,
                            TviscTrans_K=TviscTrans_K)
    else:
        # If no error is thrown, it's a mixed phase - create mixed EOS (allow MixedEOSStruct errors to propagate)
        iceEOS = MixedEOSStruct(P_MPa, T_K, phaseStr, mixParameters,
                                  porosType=porosType, phiTop_frac=phiTop_frac, 
                                  Pclosure_MPa=Pclosure_MPa, phiMin_frac=phiMin_frac,
                                  EXTRAP=EXTRAP, ClathDissoc=ClathDissoc, 
                                  minPres_MPa=minPres_MPa, minTres_K=minTres_K,
                                  ICEIh_DIFFERENT=ICEIh_DIFFERENT, etaFixed_Pas=etaFixed_Pas,
                                  TviscTrans_K=TviscTrans_K)
    if iceEOS.ALREADY_LOADED:
        log.debug(f'Ice {phaseStr} EOS already loaded. Reusing existing EOS.')
        iceEOS = EOSlist.loaded[iceEOS.EOSlabel]

    # Ensure each EOSlabel is included in EOSlist, in case we have reused EOSs with
    # e.g. a smaller range that can reuse the larger-range already-loaded EOS.
    if iceEOS.EOSlabel not in EOSlist.loaded.keys():
        EOSlist.loaded[iceEOS.EOSlabel] = iceEOS

    iceEOSwrapper = EOSwrapper(iceEOS.EOSlabel)

    return iceEOSwrapper

class IceEOSStruct:
    def __init__(self, P_MPa, T_K, phaseStr, porosType=None, phiTop_frac=0, Pclosure_MPa=0,
                 phiMin_frac=0, EXTRAP=False, ClathDissoc=None, minPres_MPa=None, minTres_K=None,
                 ICEIh_DIFFERENT=False, etaFixed_Pas=None, TviscTrans_K=None):
        self.EOSlabel = f'phase{phaseStr}poros{porosType}phi{phiTop_frac}Pclose{Pclosure_MPa}' + \
                        f'phiMin{phiMin_frac}extrap{EXTRAP}etaFixed{etaFixed_Pas}' + \
                        f'TviscTrans{TviscTrans_K}'
        self.ALREADY_LOADED, self.rangeLabel, P_MPa, T_K, self.deltaP, self.deltaT \
            = CheckIfEOSLoaded(self.EOSlabel, P_MPa, T_K, minPres_MPa=minPres_MPa, minTres_K=minTres_K)
        self.Pmin = np.min(P_MPa)
        self.Pmax = np.max(P_MPa)
        self.Tmin = np.min(T_K)
        self.Tmax = np.max(T_K)
        if not self.ALREADY_LOADED:
            self.EOSdeltaP = self.deltaP
            self.EOSdeltaT = self.deltaT
            self.EXTRAP = EXTRAP
            self.EOStype = 'ice'
            log.debug(f'Loading EOS for {phaseStr} with ' +
                      f'P_MPa = [{self.Pmin:.1f}, {self.Pmax:.1f}, {self.deltaP:.3f}], ' +
                      f'T_K = [{self.Tmin:.1f}, {self.Tmax:.1f}, {self.deltaT:.3f}], ' +
                      f'for [min, max, step] with EXTRAP = {self.EXTRAP}.')

            # Make sure arrays are long enough to interpolate
            nPs = np.size(P_MPa)
            nTs = np.size(T_K)
            if(nPs <= 3):
                P_MPa = np.linspace(P_MPa[0], P_MPa[-1], nPs*3)
            if(nTs <= 3):
                T_K = np.linspace(T_K[0], T_K[-1], nTs*3)
            # If input arrays are equal length, repeat final T value due to a quirk of numpy arrays
            # combined with SeaFreeze's particular implementation that requires gridded P,T values
            # to have different array lengths
            if(nPs == nTs):
                T_K = np.append(T_K, T_K[-1]*1.00001)

            # Assign phase ID and string for convenience in functions where iceEOS is passed
            self.phaseStr = phaseStr
            self.phaseID = PhaseInv(phaseStr)

            if phaseStr == 'Clath':
                # Special functions for clathrate properties
                rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK \
                    = ClathProps(P_MPa, T_K)
                if ClathDissoc is not None and ClathDissoc.NAGASHIMA:
                    self.phase = ClathStableNagashima2017(P_MPa, T_K)
                else:
                    self.phase = ClathStableSloan1998(P_MPa, T_K)

                # Create phase finder -- note that the results from this function must be cast to int after retrieval
                # Returns either Constants.phaseClath (stable) or 0 (not stable), making it compatible with GetTfreeze
                self.ufn_phase = RGIwrap(RegularGridInterpolator((P_MPa, T_K), self.phase, method='nearest'),
                                         self.deltaP, self.deltaT)
                self.ufn_Seismic = ClathSeismic()
            else:
                # Get tabular data from SeaFreeze for all other ice phases
                # Set extrapolation boundaries to limits defined in SeaFreeze
                Pmin = {'Ih':   0, 'II': 0, 'III': 0, 'V': 0, 'VI': 0}
                Pmax = {'Ih':   400, 'II': 900, 'III': 500, 'V': 1000, 'VI': 3000}
                Tmin = {'Ih':    1.0,   'II':  0,   'III':  0, 'V': 0, 'VI': 0}
                Tmax = {'Ih':    300.99999999999983,   'II':  270.00000000000006, 'III':  270.00000000000006,
                        'V': 300.0000000000001, 'VI': 400.0000000000001}
                self.Pmin = np.maximum(self.Pmin, Pmin[self.phaseStr])
                self.Pmax = np.minimum(self.Pmax, Pmax[self.phaseStr])
                self.Tmin = np.maximum(self.Tmin, Tmin[self.phaseStr])
                self.Tmax = np.minimum(self.Tmax, Tmax[self.phaseStr])
                if np.min(P_MPa) < self.Pmin:
                    log.warning(f'Input Pmin less than SeaFreeze limit for ice {self.phaseStr}. Resetting to SF min of'
                                f' {self.Pmin} MPa.')
                    P_MPa = np.linspace(self.Pmin, np.max(P_MPa), np.size(P_MPa))
                if np.min(P_MPa) > self.Pmax:
                    # Sometimes when querying for HP ices, the pressure input does not make sense for the given ice and sfz returns np.nan for this phase.
                    # In this case, we should reset np.min(P_MPa) to just slightly below Pmax so we do not get np.nan returned
                    log.warning(f'Input Pmin is greater than the SeaFreeze limit for ice {self.phaseStr}. Resetting to SF max of'
                                f' {self.Pmax} MPa.')
                    P_MPa = np.linspace(self.Pmax*0.99999, np.max(P_MPa), np.size(P_MPa))
                if np.max(P_MPa) > self.Pmax:
                    log.warning(f'Input Pmax greater than SeaFreeze limit for ice {self.phaseStr}. Resetting to SF '
                                f'max of {self.Pmax} MPa.')
                    P_MPa = np.linspace(np.min(P_MPa), self.Pmax, np.size(P_MPa))
                if np.min(T_K) < self.Tmin:
                    log.warning(f'Input Tmin less than SeaFreeze limit for ice {self.phaseStr}. Resetting to SF min of'
                                f' {self.Tmin} K.')
                    T_K = np.linspace(self.Tmin, np.max(T_K), np.size(T_K))

                if np.max(T_K) > self.Tmax:
                    log.warning(f'Input Tmax greater than SeaFreeze limit for ice {self.phaseStr}. Resetting to SF '
                                f'max of'
                                f' {self.Tmax} K.')
                    T_K = np.linspace(np.min(T_K), self.Tmax, np.size(T_K))
                if (T_K[-1] - T_K[0]) < self.EOSdeltaT:
                    # Sometimes when querying for HP ices, we reset the input arrays below the EOS deltas specified by user. In this case, we should reset to the EOS delta value
                    T_K = np.linspace(T_K[0] - self.EOSdeltaT, T_K[-1], 4)
                    self.Tmin = np.min(T_K)
                if (P_MPa[-1] - P_MPa[0]) < self.EOSdeltaP:
                    # Sometimes when querying for HP ices, we reset the input arrays below the EOS deltas specified by user. In this case, we should reset to the EOS delta value
                    P_MPa = np.linspace(P_MPa[0] - self.EOSdeltaT, P_MPa[-1], 4)
                    self.Pmin = np.min(P_MPa)
                    log.warning(f'Input T_K or P_MPa range is less than SeaFreeze limit for ice {self.phaseStr}. Resetting to SF min of'
                                f' {self.Tmin} K and SF max of {self.Pmax} MPa.')
                    T_K = np.linspace(self.Tmin, self.Tmax, np.size(T_K))
                    P_MPa = np.linspace(self.Pmin, self.Pmax, np.size(P_MPa))
                if (T_K[0] >= T_K[-1]) or (P_MPa[0] >= P_MPa[-1]):
                    # Sometimes when querying for HP ices, the temperature or pressure reset makes the array no longer strictly increasing and rectbivariatespline requires strictly increasing inputs.
                    # In this case, we are outside the bounds of plausible range of this HP ice forming, so let's set all to np.nan
                    T_K = np.linspace(np.min(T_K)*0.9999, np.max(T_K), np.size(T_K))
                    P_MPa = np.linspace(np.min(P_MPa)*0.9999, np.max(P_MPa), np.size(P_MPa))
                    rho_kgm3 = np.zeros((np.size(P_MPa), np.size(T_K))) + np.nan
                    Cp_JkgK = rho_kgm3
                    alpha_pK = rho_kgm3
                    kTherm_WmK = rho_kgm3
                else:
                    PTgrid = sfPTgrid(P_MPa, T_K)
                    iceOut = SeaFreeze(PTgrid, phaseStr)
                    rho_kgm3 = iceOut.rho
                    Cp_JkgK = iceOut.Cp
                    alpha_pK = iceOut.alpha
                    if ICEIh_DIFFERENT and phaseStr == 'Ih':
                        kTherm_WmK = np.array([kThermIceIhWolfenbarger2021(T_K) for _ in P_MPa])
                    else:
                        kTherm_WmK = np.array([kThermIsobaricAnderssonInaba2005(T_K, PhaseInv(phaseStr)) for _ in P_MPa])
                self.ufn_Seismic = IceSeismic(phaseStr, self.EXTRAP)
                self.ufn_phase = returnVal(self.phaseID)
            # Interpolate functions for this ice phase that can be queried for properties
            self.ufn_rho_kgm3 = RectBivariateSpline(P_MPa, T_K, rho_kgm3)
            self.ufn_Cp_JkgK = RectBivariateSpline(P_MPa, T_K, Cp_JkgK)
            self.ufn_alpha_pK = RectBivariateSpline(P_MPa, T_K, alpha_pK)
            self.ufn_kTherm_WmK = RectBivariateSpline(P_MPa, T_K, kTherm_WmK)
            self.ufn_eta_Pas = ViscIceUniform_Pas(etaSet_Pas=etaFixed_Pas, TviscTrans_K=TviscTrans_K)

            if porosType is None or porosType == 'none':
                self.ufn_phi_frac = ReturnZeros(1)
                self.POROUS = False
            else:
                self.ufn_phi_frac = GetphiCalc(phiTop_frac,
                                              GetphiFunc(porosType, phiTop_frac, Pclosure_MPa, None, P_MPa, T_K),
                                              phiMin_frac)
                self.POROUS = True

            # Store complete EOSStruct in global list of loaded EOSs
            EOSlist.loaded[self.EOSlabel] = self
            EOSlist.ranges[self.EOSlabel] = self.rangeLabel

    def fn_porosCorrect(self, propBulk, propPore, phi, J):
        # Combine pore fluid properties with matrix properties in accordance with
        # Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
        return (propBulk**J * (1 - phi) + propPore**J * phi) ** (1/J)

    # Limit extrapolation to use nearest value from evaluated fit
    def fn_phase(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_phase(P_MPa, T_K, grid=grid)
    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_rho_kgm3(P_MPa, T_K, grid=grid)
    def fn_Cp_JkgK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_Cp_JkgK(P_MPa, T_K, grid=grid)
    def fn_alpha_pK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_alpha_pK(P_MPa, T_K, grid=grid)
    def fn_kTherm_WmK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_kTherm_WmK(P_MPa, T_K, grid=grid)
    def fn_phi_frac(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_phi_frac(P_MPa, T_K, grid=grid)
    def fn_Seismic(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_Seismic(P_MPa, T_K, grid=grid)
    def fn_sigma_Sm(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_sigma_Sm(P_MPa, T_K, grid=grid)
    def fn_eta_Pas(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_eta_Pas(P_MPa, T_K, grid=grid)


class MixedEOSStruct:
    """
    Mixed equation of state that combines properties from two component EOS phases 
    using specified mixing rules for each property type. Generates each component EOS and applies mixing rules to return the mixed phase properties.
    """
    def __init__(self, P_MPa, T_K, phaseStr, mixParameters, porosType=None, phiTop_frac=0, 
                 Pclosure_MPa=0, phiMin_frac=0, EXTRAP=False, ClathDissoc=None, 
                 minPres_MPa=None, minTres_K=None, ICEIh_DIFFERENT=False, 
                 etaFixed_Pas=None, TviscTrans_K=None):
        # Get the name of each EOS
        self.mixFrac = mixParameters['mixFrac']
        self.JmixedRheologyConstant = mixParameters['JmixedRheologyConstant']
        self.EOSlabel = f'phase{phaseStr}poros{porosType}phi{phiTop_frac}Pclose{Pclosure_MPa}' + \
                        f'phiMin{phiMin_frac}extrap{EXTRAP}etaFixed{etaFixed_Pas}' + \
                        f'TviscTrans{TviscTrans_K}' + \
                        f'mixFrac{self.mixFrac}' + \
                        f'JrheologyConstant{self.JmixedRheologyConstant}'
        self.ALREADY_LOADED, self.rangeLabel, P_MPa, T_K, self.deltaP, self.deltaT \
            = CheckIfEOSLoaded(self.EOSlabel, P_MPa, T_K, minPres_MPa=minPres_MPa, minTres_K=minTres_K)
        if not self.ALREADY_LOADED:
            self.phaseStr = phaseStr
            self.phaseID = PhaseInv(phaseStr)
            self.phaseOne, self.phaseTwo = MixedPhaseSeparator(phaseStr)
            # Create the component EOS structures
            self.firstEOS = IceEOSStruct(P_MPa, T_K, self.phaseOne, porosType=porosType,
                                        phiTop_frac=phiTop_frac, Pclosure_MPa=Pclosure_MPa,
                                        phiMin_frac=phiMin_frac, EXTRAP=EXTRAP,
                                        ClathDissoc=ClathDissoc, minPres_MPa=minPres_MPa,
                                        minTres_K=minTres_K, ICEIh_DIFFERENT=ICEIh_DIFFERENT,
                                        etaFixed_Pas=etaFixed_Pas, TviscTrans_K=TviscTrans_K)
            if self.firstEOS.ALREADY_LOADED:
                log.debug(f'Ice {self.phaseOne} EOS already loaded. Reusing existing EOS.')
                self.firstEOS = EOSlist.loaded[self.firstEOS.EOSlabel]
            
            self.secondEOS = IceEOSStruct(P_MPa, T_K, self.phaseTwo, porosType=porosType,
                                            phiTop_frac=phiTop_frac, Pclosure_MPa=Pclosure_MPa,
                                            phiMin_frac=phiMin_frac, EXTRAP=EXTRAP,
                                            ClathDissoc=ClathDissoc, minPres_MPa=minPres_MPa,
                                            minTres_K=minTres_K, ICEIh_DIFFERENT=ICEIh_DIFFERENT,
                                            etaFixed_Pas=etaFixed_Pas, TviscTrans_K=TviscTrans_K)
            if self.secondEOS.ALREADY_LOADED:
                log.debug(f'Ice {self.phaseTwo} EOS already loaded. Reusing existing EOS.')
                self.secondEOS = EOSlist.loaded[self.secondEOS.EOSlabel]
            # Set the limits to the minimal bounds of the component EOSs
            self.Pmin = self.firstEOS.Pmin if self.firstEOS.Pmin > self.secondEOS.Pmin else self.secondEOS.Pmin
            self.Pmax = self.firstEOS.Pmax if self.firstEOS.Pmax < self.secondEOS.Pmax else self.secondEOS.Pmax
            self.Tmin = self.firstEOS.Tmin if self.firstEOS.Tmin > self.secondEOS.Tmin else self.secondEOS.Tmin
            self.Tmax = self.firstEOS.Tmax if self.firstEOS.Tmax < self.secondEOS.Tmax else self.secondEOS.Tmax
            self.deltaP = np.maximum(self.firstEOS.deltaP, self.secondEOS.deltaP)
            self.deltaT = np.maximum(self.firstEOS.deltaT, self.secondEOS.deltaT)
            self.EOSdeltaP = self.deltaP
            self.EOSdeltaT = self.deltaT
            self.EXTRAP = EXTRAP
            self.EOStype = 'ice'
            self.POROUS = self.firstEOS.POROUS and self.secondEOS.POROUS

                    
            log.debug(f'Created mixed EOS "{phaseStr}": {self.mixFrac}% {self.phaseTwo} + ' +
                    f'{1-self.mixFrac}% {self.phaseOne}')
            
            # Store complete EOSStruct in global list of loaded EOSs
            EOSlist.loaded[self.EOSlabel] = self
            EOSlist.ranges[self.EOSlabel] = self.rangeLabel
    
    def _mixingRule(self, prop1, prop2, rule):
        if rule == 'arithmetic':
            if isinstance(prop1, tuple):
                return tuple(array * self.mixFrac + array2 * (1 - self.mixFrac) for array, array2 in zip(prop1, prop2))
            else:
                return (prop1 * self.mixFrac + prop2 * (1 - self.mixFrac))
        elif rule == 'VRHAverage':
            if isinstance(prop1, tuple):
                M_V = tuple(array * self.mixFrac + array2 * (1 - self.mixFrac) for array, array2 in zip(prop1, prop2))
                M_R_Inverse = tuple(self.mixFrac / (array) + (1 - self.mixFrac) / (array2) for array, array2 in zip(prop1, prop2))
                M_R = tuple(1 / inv for inv in M_R_Inverse)
                return tuple((v + r) / 2 for v, r in zip(M_V, M_R))
            else:
                M_V = prop1 * self.mixFrac + prop2 * (1 - self.mixFrac)
                M_R_Inverse = self.mixFrac / (prop1) + (1 - self.mixFrac) / (prop2)
                M_R = 1 / M_R_Inverse
                return (M_V + M_R) / 2
        elif rule == 'Carahan2004Averaging':
            if self.JmixedRheologyConstant == 0:
                # Geometric mean for J = 0
                if isinstance(prop1, tuple):
                    return tuple(array ** self.mixFrac * array2 ** (1 - self.mixFrac) for array, array2 in zip(prop1, prop2))
                else:
                    return prop1 ** self.mixFrac * prop2 ** (1 - self.mixFrac)
            else:
                if isinstance(prop1, tuple):
                    return tuple((array ** self.JmixedRheologyConstant * self.mixFrac + array2 ** self.JmixedRheologyConstant * (1 - self.mixFrac)) ** (1 / self.JmixedRheologyConstant) for array, array2 in zip(prop1, prop2))
                else:
                    return (prop1 ** self.JmixedRheologyConstant * self.mixFrac + prop2 ** self.JmixedRheologyConstant * (1 - self.mixFrac)) ** (1 / self.JmixedRheologyConstant)
        else:
            raise ValueError(f'Invalid mixing rule: {rule}')

    
    def fn_phase(self, P_MPa, T_K, grid=False):
        """ For mixed phases, we call the phase function from both component EOS and add them together. This takes advantage of the two facts:
        1) that the phaseID of the mixed phase is the sum of the phaseIDs of the component phases.
        2)  When one phase becomes unstable, the fn_phase returns zero for the unstable phase and the phaseID of the mixed phase will be the phaseID of the only stable phase."""
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        #TODO Fix this implementation so that it uses oceanEOS.fn_phase
        phase1 = self.firstEOS.fn_phase(P_MPa, T_K, grid)
        phase2 = self.secondEOS.fn_phase(P_MPa, T_K, grid)
        if phase1 == 0:
            return phase1
        else:
            return phase1 + phase2
    
    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        rho1 = self.firstEOS.fn_rho_kgm3(P_MPa, T_K, grid)
        rho2 = self.secondEOS.fn_rho_kgm3(P_MPa, T_K, grid)
        return self._mixingRule(rho1, rho2, 'arithmetic')
    
    def fn_Cp_JkgK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        Cp1 = self.firstEOS.fn_Cp_JkgK(P_MPa, T_K, grid)
        Cp2 = self.secondEOS.fn_Cp_JkgK(P_MPa, T_K, grid)
        return self._mixingRule(Cp1, Cp2, 'arithmetic')
    
    def fn_alpha_pK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        alpha1 = self.firstEOS.fn_alpha_pK(P_MPa, T_K, grid)
        alpha2 = self.secondEOS.fn_alpha_pK(P_MPa, T_K, grid)
        return self._mixingRule(alpha1, alpha2, 'arithmetic')
    
    def fn_kTherm_WmK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        kTherm1 = self.firstEOS.fn_kTherm_WmK(P_MPa, T_K, grid)
        kTherm2 = self.secondEOS.fn_kTherm_WmK(P_MPa, T_K, grid)
        return self._mixingRule(kTherm1, kTherm2, 'arithmetic')
    
    def fn_phi_frac(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        phi1 = self.firstEOS.fn_phi_frac(P_MPa, T_K, grid)
        phi2 = self.secondEOS.fn_phi_frac(P_MPa, T_K, grid)
        return self._mixingRule(phi1, phi2, 'arithmetic')
    
    def fn_Seismic(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        seismic1 = self.firstEOS.fn_Seismic(P_MPa, T_K, grid)
        seismic2 = self.secondEOS.fn_Seismic(P_MPa, T_K, grid)
        return self._mixingRule(seismic1, seismic2, 'VRHAverage')
    
    def fn_sigma_Sm(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        sigma1 = self.firstEOS.fn_sigma_Sm(P_MPa, T_K, grid)
        sigma2 = self.secondEOS.fn_sigma_Sm(P_MPa, T_K, grid)
        return self._mixingRule(sigma1, sigma2, 'arithmetic')
    
    def fn_eta_Pas(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        eta1 = self.firstEOS.fn_eta_Pas(P_MPa, T_K, grid)
        eta2 = self.secondEOS.fn_eta_Pas(P_MPa, T_K, grid)
        return self._mixingRule(eta1, eta2, 'Carahan2004Averaging')
    
    def fn_averageValuesAccordingtoRule(self, prop1, prop2, rule):
        if rule == 'arithmetic':
            return self._mixingRule(prop1, prop2, 'arithmetic')
        elif rule == 'VRHAverage':
            return self._mixingRule(prop1, prop2, 'VRHAverage')
        elif rule == 'Carahan2004Averaging':
            return self._mixingRule(prop1, prop2, 'Carahan2004Averaging')


class returnVal:
    def __init__(self, val):
        self.val = val
    def __call__(self, P, T, grid=False):
        if grid:
            P, _ = np.meshgrid(P, T, indexing='ij')
        return (np.ones_like(P) * self.val).astype(np.int_)

def CheckIfEOSLoaded(EOSlabel, P_MPa, T_K, FORCE_NEW=False, minPres_MPa=None, minTres_K=None):
    """ Determine if we need to load a new EOS, or if we can reuse one that's already been
        loaded within this session.

        Args:
            EOSlabel (str): A unique identifier containing all the settings passed to the
                EOSStruct class instantiator.
            P_MPa (float, shape N): Pressures to desired for constructing the EOS in MPa.
            T_K (float, shape N): Temperatures to desired for constructing the EOS in K.
            FORCE_NEW = False (bool): Whether to force a reload each time, instead of checking.
                Overwrites any previously loaded EOS in the EOSlist that has the same EOSlabel.
        Returns:
            ALREADY_LOADED (bool): Whether we can make use of an EOSStruct in EOSlist.loaded
                with a label matching the one we wish to load now.
            rangeLabel (str): A string identifying the min/max/step values for P and T.
            outP_MPa (float, shape N): Pressures to use for constructing the EOS in MPa.
            outT_K (float, shape N): Temperatures to use for constructing the EOS in K.
    """

    # Create label for identifying P, T arrays
    deltaP = np.maximum(np.round(np.mean(np.diff(P_MPa)), 2), 0.001)
    deltaT = np.maximum(np.round(np.mean(np.diff(T_K)), 2), 0.001)
    Pmin = np.min(P_MPa)
    Pmax = np.max(P_MPa)
    Tmin = np.min(T_K)
    Tmax = np.max(T_K)
    rangeLabel = f'{Pmin:.2f},{Pmax:.2f},{deltaP:.2e},' + \
                 f'{Tmin:.3f},{Tmax:.3f},{deltaT:.2e}'
    if (not FORCE_NEW) and EOSlabel in EOSlist.loaded.keys():
        if EOSlist.ranges[EOSlabel] == rangeLabel:
            # This exact EOS has been loaded already. Reuse the one in memory
            ALREADY_LOADED = True
            outP_MPa = np.array([Pmin, Pmax])
            outT_K = np.array([Tmin, Tmax])
        else:
            # Check if we can reuse an already-loaded EOS because the
            # P, T ranges are contained within the already-loaded EOS
            EOSPmin = EOSlist.loaded[EOSlabel].Pmin
            EOSPmax = EOSlist.loaded[EOSlabel].Pmax
            EOSTmin = EOSlist.loaded[EOSlabel].Tmin
            EOSTmax = EOSlist.loaded[EOSlabel].Tmax
            nopeP = np.min(P_MPa) < EOSPmin * 0.9 or \
                    np.max(P_MPa) > EOSPmax * 1.1 or \
                    deltaP < EOSlist.loaded[EOSlabel].deltaP * 0.9
            nopeT = np.min(T_K) < EOSTmin - 0.1 or \
                    np.max(T_K) > EOSTmax + 0.1 or \
                    deltaT < EOSlist.loaded[EOSlabel].deltaT * 0.9
            if nopeP or nopeT:
                # The new inputs have at least one min/max value outside the range
                # of the previously loaded EOS, so we have to load a new one.
                ALREADY_LOADED = False
                # Set P and T ranges to include the outer bounds from
                # the already-loaded EOS and the one we want now
                minPmin = np.minimum(np.min(P_MPa), EOSlist.loaded[EOSlabel].Pmin)
                maxPmax = np.maximum(np.max(P_MPa), EOSlist.loaded[EOSlabel].Pmax)
                minTmin = np.minimum(np.min(T_K), EOSlist.loaded[EOSlabel].Tmin)
                maxTmax = np.maximum(np.max(T_K), EOSlist.loaded[EOSlabel].Tmax)
                deltaP = np.round(np.minimum(np.mean(np.diff(P_MPa)), EOSlist.loaded[EOSlabel].deltaP), 2)
                deltaT = np.round(np.minimum(np.mean(np.diff(T_K)), EOSlist.loaded[EOSlabel].deltaT), 2)
                if deltaP == 0: deltaP = 0.01
                if deltaT == 0: deltaT = 0.01
                if minPres_MPa is not None and deltaP < minPres_MPa:
                    log.warning(f'deltaP of {deltaP:.2f} MPa less than minimum res setting of {minPres_MPa}. Resetting to {minPres_MPa}.')
                    deltaP = minPres_MPa
                if minTres_K is not None and deltaT < minTres_K:
                    log.warning(f'deltaT of {deltaT:.2f} K less than minimum res setting of {minTres_K}. Resetting to {minTres_K}.')
                    deltaT = minTres_K
                nPs = int((maxPmax - minPmin) / deltaP)
                nTs = int((maxTmax - minTmin) / deltaT)
                # Ensure we don't have too few points that we error later
                if nPs < 5:
                    nPs = 5
                if nTs < 5:
                    nTs = 5
                outP_MPa = np.linspace(minPmin, maxPmax, nPs)
                outT_K = np.linspace(minTmin, maxTmax, nTs)
                rangeLabel = f'{np.min(outP_MPa):.2f},{np.max(outP_MPa):.2f},{deltaP:.2e},' + \
                             f'{np.min(outT_K):.3f},{np.max(outT_K):.3f},{deltaT:.2e}'
            else:
                # A previous EOS has been loaded that has a wider P or T range than the inputs,
                # so we will use the previously loaded one.
                ALREADY_LOADED = True
                outP_MPa = np.array([EOSPmin, EOSPmax])
                outT_K = np.array([EOSTmin, EOSTmax])
    else:
        # This EOS has not been loaded, so we need to load it with the input parameters
        ALREADY_LOADED = False
        # Override min resolution if set
        if minPres_MPa is not None and deltaP < minPres_MPa:
            log.warning(f'deltaP of {deltaP:.2f} MPa less than minimum res setting of {minPres_MPa}. Resetting to {minPres_MPa}.')
            deltaP = minPres_MPa
        if minTres_K is not None and deltaT < minTres_K:
            log.warning(f'deltaT of {deltaT:.2f} K less than minimum res setting of {minTres_K}. Resetting to {minTres_K}.')
            deltaT = minTres_K
        rangeLabel = f'{Pmin:.2f},{Pmax:.2f},{deltaP:.2e},' + \
                     f'{Tmin:.3f},{Tmax:.3f},{deltaT:.2e}'
        # Ensure that P_MPa is strictly ascending, namely that {Pmin and Pmax are not the same
        # If so, then increment Pmax just slightly
        if Pmax - Pmin == 0:
            Pmax = Pmax + 0.001
            deltaP = 0.001
        if Tmax - Tmin == 0:
            Tmax = Tmin + 0.001
            deltaT = 0.001
        # Use of np.arange would be simpler here, but can cause errors when loading an EOS for a thin layer, e.g. for
        # some cases with ice VI inside pores.
        nPs = np.maximum(round(abs(Pmax-Pmin)/deltaP), 10)
        outP_MPa = np.linspace(Pmin, Pmax, nPs)
        # outP_MPa = np.arange(Pmin, Pmax, np.maximum(deltaP)
        nTs = np.maximum(round(abs(Tmax-Tmin)/deltaT), 11)
        # To prevent Seafreeze errors, ensure that nTs and nPs don't have same size
        if nTs == nPs:
            nTs = nTs + 1
        outT_K = np.linspace(Tmin, Tmax, nTs)
    return ALREADY_LOADED, rangeLabel, outP_MPa, outT_K, deltaP, deltaT


# Create a function that can pack up (P,T) pairs that are compatible with SeaFreeze
def sfPTpairs(P_MPa, T_K):
    return np.array([(P, T) for P, T in zip(P_MPa, T_K)], dtype='f,f').astype(object)
# Same for PTm triplets
def sfPTmTrips(P_MPa, T_K, m_molal):
    return np.array([(P, T, m_molal) for P, T in zip(P_MPa, T_K)], dtype='f,f,f').astype(object)
# Same as above but for a grid
def sfPTgrid(P_MPa, T_K):
    return np.array([P_MPa, T_K], dtype=object)
# Same for PTm grid
def sfPTmGrid(P_MPa, T_K, m_molal):
    return np.array([P_MPa, T_K, np.array([m_molal])], dtype=object)

# Create callable class to act as a wrapper for SeaFreeze phase lookup
class SFphase:
    def __init__(self, w_ppt, comp):
        self.w_ppt = w_ppt
        if comp == 'PureH2O':
            self.comp = 'water1'
            self.m_molal = np.nan
            self.path = None
        else:
            self.comp = comp
            self.m_molal = Ppt2molal(self.w_ppt, Constants.m_gmol[self.comp])
            self.path = SFmatPath[self.comp]

    def PTpairs(self, Pin, Tin):
        if np.size(Pin) == 1 and np.size(Tin) == 1:
            return np.array([(Pin, Tin)], dtype='f,f').astype(object)
        elif np.size(Pin) == 1:
            return np.array([(Pin, T) for T in Tin], dtype='f,f').astype(object)
        elif np.size(Tin) == 1:
            return np.array([(P, Tin) for P in Pin], dtype='f,f').astype(object)
        elif np.size(Pin) == np.size(Tin):
            return sfPTpairs(Pin, Tin)
        else:
            log.warning('2D array as input to SeaFreeze phase finder when 1D array was expected. A 2D array will be output.')
            return sfPTgrid(Pin, Tin)

    def PTmTrips(self, Pin, Tin):
        if np.size(Pin) == 1 and np.size(Tin) == 1:
            return np.array([(Pin, Tin, self.m_molal)], dtype='f,f').astype(object)
        elif np.size(Pin) == 1:
            return np.array([(Pin, T, self.m_molal) for T in Tin], dtype='f,f').astype(object)
        elif np.size(Tin) == 1:
            return np.array([(P, Tin, self.m_molal) for P in Pin], dtype='f,f').astype(object)
        elif np.size(Pin) == np.size(Tin):
            return sfPTmTrips(Pin, Tin, self.m_molal)
        else:
            log.warning('2D array as input to SeaFreeze phase finder when 1D array was expected. A 2D array will be output.')
            return sfPTmGrid(Pin, Tin, self.m_molal)

    def __call__(self, P_MPa, T_K, grid=False):
        if self.comp == 'water1':
            if grid:
                PT = sfPTgrid(P_MPa, T_K)
            else:
                PT = self.PTpairs(P_MPa, T_K)
            return WhichPhase(PT).astype(np.int_)
        else:
            if grid:
                PTm = sfPTmGrid(P_MPa, T_K, self.m_molal)
            else:
                PTm = self.PTmTrips(P_MPa, T_K)
            return WhichPhase(PTm, solute=self.comp).astype(np.int_)
    
class RGIwrap:
    """ Creates a wrapper for passing P, T arguments as a tuple to RegularGridInterpolator """
    def __init__(self, RGIfunc, deltaP, deltaT):
        self.RGIfunc = RGIfunc
        self.deltaP = deltaP
        self.deltaT = deltaT
        
    def __call__(self, P_MPa, T_K, grid=False):
        if grid:
            P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
        return (self.RGIfunc((P_MPa, T_K)))
        

class SFSeismic:
    """ Creates a function call for returning seismic properties of depth profile for SeaFreeze solutions. """
    def __init__(self, compstr, P_MPa, T_K, seaOut, wOcean_ppt, EXTRAP):
        self.comp = compstr
        self.w_ppt = wOcean_ppt
        self.EXTRAP = EXTRAP

        self.ufn_VP_kms = RectBivariateSpline(P_MPa, T_K, seaOut.vel * 1e-3)
        self.ufn_KS_GPa = RectBivariateSpline(P_MPa, T_K, seaOut.Ks * 1e-3)

    def __call__(self, P_MPa, T_K, grid=False):
        return self.ufn_VP_kms(P_MPa, T_K, grid=grid), self.ufn_KS_GPa(P_MPa, T_K, grid=grid)


class H2Osigma_Sm:
    def __init__(self, sigmaSet_Sm=None):
        if sigmaSet_Sm is None:
            self.sigma_Sm = Constants.sigmaH2O_Sm
        else:
            self.sigma_Sm = sigmaSet_Sm

    def __call__(self, P_MPa, T_K, grid=False):
        if grid:
            return np.zeros((np.size(P_MPa), np.size(T_K))) + self.sigma_Sm
        else:
            return np.zeros_like(P_MPa) + self.sigma_Sm


class IceSeismic:
    def __init__(self, phaseStr, EXTRAP):
        self.phase = phaseStr
        self.EXTRAP = EXTRAP

    def __call__(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            # Set extrapolation boundaries to limits defined in SeaFreeze
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, 0, 3000, 0, 400)
        if grid:
            PT = sfPTgrid(P_MPa, T_K)
        else:
            PT = sfPTpairs(P_MPa, T_K)
        seaOut = SeaFreeze(PT, self.phase)
        return seaOut.Vp * 1e-3, seaOut.Vs * 1e-3,  seaOut.Ks * 1e-3, seaOut.shear * 1e-3


def GetPfreeze(oceanEOS, phaseTop, Tb_K, PLower_MPa=0.1, PUpper_MPa=300, PRes_MPa=0.1, UNDERPLATE=None, HPNOOCEAN = None,
               ALLOW_BROKEN_MODELS=False, DO_EXPLOREOGRAM=False):
    """ Returns the pressure at which ice changes phase based on temperature, salinity, and composition

        Args:
            oceanEOS (OceanEOSStruct): Interpolator functions for evaluating the ocean EOS
            Tb_K (float): Temperature of the phase transition in K
            HPNOOCEAN (bool): Whether to find the pressure transition to high pressure ice (just like we do for underplating)
        Returns:
            Pfreeze_MPa (float): Pressure at the phase change interface consistent with Tb_K
    """
    phaseChangeUnderplateOrHPNoOcean = lambda P: 0.5 + (phaseTop - oceanEOS.fn_phase(P, Tb_K))
    if (UNDERPLATE or HPNOOCEAN) is None:
        raise ValueError('UNDERPLATE or HPNOOCEAN is not set. Please set one of these to True or False.')
        TRY_BOTH = True
        UNDERPLATEOrHPNOOCEAN = True
    else:
        TRY_BOTH = False
        UNDERPLATEOrHPNOOCEAN = UNDERPLATE or HPNOOCEAN
    if UNDERPLATEOrHPNOOCEAN:
        phaseChange = phaseChangeUnderplateOrHPNoOcean
    else:
        phaseChange = lambda P: 0.5 - (phaseTop - oceanEOS.fn_phase(P, Tb_K))

    Pfreeze_MPa = None
    if phaseChange(PLower_MPa) * phaseChange(PUpper_MPa) > 0:
        # GetZero will error out in this case. Use a brute force strategy instead
        Pvals = np.arange(PLower_MPa, PUpper_MPa, PRes_MPa)
        changes = phaseChange(Pvals)
        if np.any(changes > 0):
            iChange = np.where(changes < 0)[0]
            if np.size(iChange) != 0:
                Pfreeze_MPa = Pvals[iChange[0]] + PRes_MPa/5

    if Pfreeze_MPa is None:
        try:
            Pfreeze_MPa = GetZero(phaseChange, bracket=[PLower_MPa, PUpper_MPa]).root + PRes_MPa/5
        except ValueError:
            if UNDERPLATEOrHPNOOCEAN:
                msg = f'Tb_K of {Tb_K:.3f} is not consistent with explicit underplating or automatic high pressure no ocean ice III; ' + \
                      f'the phases at the top and bottom of this range are ' + \
                      f'{PhaseConv(oceanEOS.fn_phase(PLower_MPa, Tb_K))} and ' + \
                      f'{PhaseConv(oceanEOS.fn_phase(PUpper_MPa, Tb_K))}, respectively.'
                if ALLOW_BROKEN_MODELS:
                    if DO_EXPLOREOGRAM:
                        log.info(msg)
                    else:
                        log.error(msg)
                    Pfreeze_MPa = np.nan
                else:
                    raise ValueError(msg)
            elif TRY_BOTH:
                try:
                    Pfreeze_MPa = GetZero(phaseChangeUnderplateOrHPNoOcean, bracket=[PLower_MPa, PUpper_MPa]).root + PRes_MPa / 5
                except ValueError:
                    msg = f'No transition pressure was found below {PUpper_MPa:.3f} MPa ' + \
                          f'for ice {PhaseConv(phaseTop)}. Increase PUpper_MPa until one is found.'
                    if ALLOW_BROKEN_MODELS:
                        if DO_EXPLOREOGRAM:
                            log.info(msg)
                        else:
                            log.error(msg)
                        Pfreeze_MPa = np.nan
                    else:
                        raise ValueError(msg)
            else:
                msg = f'No transition pressure was found below {PUpper_MPa:.3f} MPa ' + \
                      f'for ice {PhaseConv(phaseTop)} and with explicit underplating and automatic high pressure is explicitly set to False.'
                if DO_EXPLOREOGRAM:
                    log.info(msg)
                else:
                    log.warning(msg)
                Pfreeze_MPa = np.nan

    return Pfreeze_MPa


def GetTfreeze(oceanEOS, P_MPa, T_K, TfreezeRange_K=50, TRes_K=0.05):
    """ Returns the temperature at which a solid layer melts based on temperature, salinity, and composition

        Args:
            oceanEOS (OceanEOSStruct): Interpolator functions for evaluating the ocean EOS
            P_MPa (float): Pressure of the fluid in MPa
            T_K (float): Temperature of the fluid in K
        Returns:
            Tfreeze_K (float): Temperature of nearest higher-temperature solid-liquid phase transition
                at this pressure
    """
    topPhase = oceanEOS.fn_phase(P_MPa, T_K)
    if topPhase == 0:
        log.warning('Attempting to get phase change from liquid to solid, not solid to liquid as expected.')
    phaseChange = lambda T: 0.5 - (1 - int(oceanEOS.fn_phase(P_MPa, T) > 0))

    try:
        Tfreeze_K = GetZero(phaseChange, bracket=[T_K, T_K+TfreezeRange_K]).root + TRes_K/5
    except ValueError:
        raise ValueError(f'No melting temperature was found above {T_K:.3f} K ' +
                         f'for ice {PhaseConv(topPhase)} at pressure {P_MPa:.3f} MPa. ' +
                          'Check to see if T_K is close to default Ocean.THydroMax_K value. ' +
                          'If so, increase Ocean.THydroMax_K. Otherwise, increase TfreezeRange_K ' +
                          'until a melting temperature is found.')

    return Tfreeze_K


def kThermIsobaricAnderssonInaba2005(T_K, phase):
    """ Calculate thermal conductivity of ice at a fixed pressure according to
        Andersson and Inaba (2005) as a function of temperature.
        See https://doi.org/10.1039/B500373C.
        Range of validity is as follows:
        Phase:  P (MPa):    T range (K):
        Ih      0.1         40-180*
        II      240         120-240
        III     240         180-250
        V       530         240-270
        VI      1000        135-250
        *Andersson and Inaba give an alternate equation that accounts for the range 180-273 K
        for ice Ih at 0.1 MPa, but as this was not included in the Matlab version, it's
        skipped here too. This implementation does not apply at the relevant T and P values
        for icy moon shells except at specific points, so a more versatile and accurate
        model should be found and used to replace this.

        Args:
            T_K (float, shape N): Temperatures to evaluate in K
            phase (int): Phase ID
        Returns:
            kTherm_WmK (float, shape N): Thermal conductivity of desired phase at specified temperatures
                in W/(m K)
    """
    D = np.array([np.nan, 630, 695, 93.2, np.nan, 38.0, 50.9])
    X = np.array([np.nan, 0.995, 1.097, 0.822, np.nan, 0.612, 0.612])

    kTherm_WmK = D[abs(phase)] * T_K**(-X[abs(phase)])

    return kTherm_WmK


def kThermIsothermalAnderssonInaba2005(P_MPa, phase):
    """ Calculate thermal conductivity of ice at a fixed temperature according to
        Andersson and Inaba (2005) as a function of pressure.
        See https://doi.org/10.1039/B500373C.
        Range of validity is as follows:
        Phase:  P range (GPa):  T (K):
        Ih      0-0.5           130
        II      0-0.24          120
        III     0.2-0.35        240
        V       0.35-0.6        246
        VI      0.7-2.0         246
        This implementation does not apply at the relevant T and P values for icy moon
        shells except at specific points, so a more versatile and accurate model should
        be found and used to replace this.

        Args:
            P_MPa (float, shape N): Pressure to evaluate in MPa
            phase (int, shape N): Phase index
        Returns:
            kTherm_WmK (float, shape N): Thermal conductivity of desired phase at specified pressures
                in W/(m K)
    """
    E = np.array([np.nan, 1.60, 1.25, -0.02, np.nan, 0.16, 0.37])
    F = np.array([np.nan, -0.44, 0.2, 0.2, np.nan, 0.2, 0.16])

    # Note the 1e-3 factor because F has units of 1/GPa
    kTherm_WmK = np.exp(E[abs(phase)] + F[abs(phase)] * P_MPa * 1e-3)

    return kTherm_WmK


def kThermMelinder2007(T_K, Tmelt_K, ko_WmK=2.21, dkdT_WmK2=-0.012):
    """ Calculate thermal conductivity of ice Ih according to Melinder (2007).

        Args:
            T_K (float, shape N): Temperature in K
            Tmelt_K (float, shape N): Melting temperature at the evaluated pressure in K
            ko_WmK = 2.21 (float): Thermal conductivity at the melting temperature in W/(m K)
            dkdT_WmK2 = -0.012 (float): Constant temperature derivative of k in W/(mK^2)
        Returns:
            kTherm_WmK (float, shape N): Thermal conductivity of ice Ih at specified temperature
                in W/(m K)
    """

    kTherm_WmK = ko_WmK + dkdT_WmK2 * (T_K - Tmelt_K)
    return kTherm_WmK


def kThermHobbs1974(T_K):
    """ Calculate thermal conductivity of ice Ih according to Hobbs (1974), as
        reported by Ojakangas and Stevenson (1989).

        Args:
            T_K (float, shape N): Temperature value(s) in K
        Returns:
            kTherm_WmK (float, shape N): Thermal conductivities in W/(m K)
    """
    a0 = 4.68e4  # Units of ergs/(K cm s)
    a1 = 4.88e7  # Units of ergs/(cm s)
    a0_SI = a0 * Constants.erg2J * 1e2
    a1_SI = a1 * Constants.erg2J * 1e2
    kTherm_WmK = a1_SI/T_K + a0_SI

    return kTherm_WmK


def kThermIceIhWolfenbarger2021(T_K):
    """ Calculate thermal conductivity of ice Ih at a fixed pressure according to
        Wolfenbarger et al. (2021) as a function of temperature.
        See https://doi.org/10.1016/j.dib.2021.107079.
        Wolgenbarger et al. consider a variety of sources, but for only ice Ih.
        They fit a curve that grossly improves on the residuals of the data,
        primarily in the region of interest for convecting shells, around 220-250 K.

        Args:
            T_K (float, shape N): Temperatures to evaluate in K
        Returns:
            kTherm_WmK (float, shape N): Thermal conductivity of desired phase at specified temperatures
                in W/(m K)
    """

    return 612 / T_K


class ViscOceanUniform_Pas:
    def __init__(self, etaSet_Pas=None, comp=None):
        if etaSet_Pas is None:
            if comp == 'Seawater':
                self.eta_Pas = Constants.etaSeawater_Pas
            else:
                self.eta_Pas = Constants.etaH2O_Pas
        else:
            self.eta_Pas = etaSet_Pas

    def __call__(self, P_MPa, T_K, grid=False):
        if grid:
            return np.zeros((np.size(P_MPa), np.size(T_K))) + self.eta_Pas
        else:
            return np.zeros_like(P_MPa) + self.eta_Pas


class ViscIceUniform_Pas:
    def __init__(self, etaSet_Pas=None, TviscTrans_K=None):
        if etaSet_Pas is None:
            self.eta_Pas = Constants.etaIce_Pas
        else:
            self.eta_Pas = etaSet_Pas

        if TviscTrans_K is None:
            self.TviscTrans_K = Constants.TviscIce_K
        else:
            self.TviscTrans_K = TviscTrans_K

    def __call__(self, P_MPa, T_K, grid=False):
        Ttrans_K = np.insert([0.0, np.inf], 1, self.TviscTrans_K)
        if grid:
            eta_Pas = np.zeros((np.size(P_MPa), np.size(T_K)))
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[:, np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas
        else:
            eta_Pas = np.zeros_like(P_MPa)
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas

        return eta_Pas

