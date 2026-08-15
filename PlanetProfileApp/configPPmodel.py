""" Default ice model parameters for alternative (non-self-consistent) ice layer modeling """
from PlanetProfile.Utilities.defineStructs import ModelSubstruct

configModelVersion = 1  # Integer number for config file version. Increment when new settings are added to the default config file.

def modelAssign():
    ModelParams = ModelSubstruct()

    # Specify parameters for each ice phase
    # Each parameter can be specified as a constant value, a profile (list/array), or None (to use default calculation)
    # Keys must match the phase names used in the codebase: "Ih", "III", "V"
    
    # Temperature profiles (in K)
    # If None, the default conductive profile will be used
    ModelParams.T_K = {
        "Ih": None,   # Default conductive profile: T_K = Tb_K**(Pratios) * Tsurf_K**(1 - Pratios)
        "III": None,  # If specified, provide array of length nIceIIILitho or a constant value
        "V": None     # If specified, provide array of length nIceVLitho or a constant value
    }
    
    # Pressure profiles (in MPa)
    # If None, the default linear profile will be used
    ModelParams.P_MPa = {
        "Ih": None,   # Default: linear from Psurf_MPa to PbI_MPa
        "III": None,  # Default: linear from PbI_MPa to PbIII_MPa
        "V": None     # Default: linear from PbIII_MPa to PbV_MPa
    }
    
    # Density profiles (in kg/m³)
    # If None, values will be calculated from the EOS based on P and T
    ModelParams.rho_kgm3 = {
        "Ih": None,
        "III": None,
        "V": None
    }
    
    # Heat capacity profiles (in J/kg/K)
    # If None, values will be calculated from the EOS based on P and T
    ModelParams.Cp_JkgK = {
        "Ih": None,
        "III": None,
        "V": None
    }
    
    # Thermal expansivity profiles (in 1/K)
    # If None, values will be calculated from the EOS based on P and T
    ModelParams.alpha_pK = {
        "Ih": None,
        "III": None,
        "V": None
    }
    
    # Thermal conductivity profiles (in W/m/K)
    # If None, values will be calculated based on default models
    ModelParams.kTherm_WmK = {
        "Ih": None,
        "III": None,
        "V": None
    }
    
    # Porosity profiles (as fractions, not percentages)
    # If None, values will be calculated based on default porosity models if POROUS_ICE is True
    ModelParams.phi_frac = {
        "Ih": None,
        "III": None,
        "V": None
    }
    
    # Electrical conductivity profiles (in S/m)
    # If None, values will be calculated based on default conductivity models if CALC_CONDUCT is True
    ModelParams.sigma_Sm = {
        "Ih": None,
        "III": None,
        "V": None
    }
    
    # Tidal heating profiles (in W/m³)
    # If None, tidal heating values will be calculated or applied from Planet settings
    ModelParams.Htidal_Wm3 = {
        "Ih": None,
        "III": None,
        "V": None
    }
    
    # Boundary temperature overrides (in K)
    # If None, use values from Planet.Bulk
    ModelParams.Tb_K = {
        "Ih": None,  # Same as Planet.Bulk.Tb_K
        "III": None, # Same as Planet.Bulk.TbIII_K
        "V": None    # Same as Planet.Bulk.TbV_K
    }
    
    # Boundary pressure overrides (in MPa)
    # If None, use the calculated values from phase transitions
    ModelParams.Pb_MPa = {
        "Ih": None,  # If set, overrides the calculated PbI_MPa
        "III": None, # If set, overrides the calculated PbIII_MPa
        "V": None    # If set, overrides the calculated PbV_MPa
    }
    
    # Shell thickness overrides (in km)
    # If provided, these will be used to determine layer thicknesses directly
    # instead of calculating from phase transition pressures
    ModelParams.thickness_km = {
        "Ih": None,
        "III": None,
        "V": None
    }
    
    return ModelParams 
