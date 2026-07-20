from abc import ABC, abstractmethod

class OceanEOSBase(ABC):
    """
    Abstract base class for custom ocean EOS objects.

    Subclass this, implement all abstract methods, and set the required scalar
    attributes in __init__. Then assign the instance to Planet.Ocean.EOS and set
    Planet.Do.CustomEOS['Ocean'] = True before running PlanetProfile.

    Required scalar attributes (must be set in subclass __init__):
        EOSlabel   (str)  : Unique string identifier used as the EOSlist cache key.
        EOStype    (str)  : Must be 'ocean'.
        comp       (str)  : Composition description string (e.g. species template label).
        w_ppt      (float or None): Salinity in parts per thousand, or None.
        Pmin, Pmax (float): Pressure range of this EOS in MPa.
        Tmin, Tmax (float): Temperature range of this EOS in K.
        deltaP, deltaT (float): Input grid step sizes (MPa, K). May equal EOSdeltaP/EOSdeltaT.
        EOSdeltaP, EOSdeltaT (float): EOS lookup table grid resolution (MPa, K).
        propsPmax  (float): Maximum pressure for property function evaluation in MPa.
        EXTRAP     (bool) : Whether extrapolation beyond table bounds is permitted.
        fn_porosCorrect: Set to None for non-porous ocean EOS.
    """

    fn_porosCorrect = None  # Set to None for standard ocean; override for porous EOS.

    @abstractmethod
    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        """Return density in kg/m³.
        Shape is (nP, nT) when grid=True, (n,) when grid=False.
        """

    @abstractmethod
    def fn_Cp_JkgK(self, P_MPa, T_K, grid=False):
        """Return isobaric heat capacity in J/(kg·K).
        Shape is (nP, nT) when grid=True, (n,) when grid=False.
        """

    @abstractmethod
    def fn_alpha_pK(self, P_MPa, T_K, grid=False):
        """Return thermal expansivity in 1/K.
        Shape is (nP, nT) when grid=True, (n,) when grid=False.
        """

    @abstractmethod
    def fn_kTherm_WmK(self, P_MPa, T_K, grid=False):
        """Return thermal conductivity in W/(m·K).
        Shape is (nP, nT) when grid=True, (n,) when grid=False.
        """

    @abstractmethod
    def fn_phase(self, P_MPa, T_K, grid=False):
        """Return integer phase at each (P, T) point per PhaseConv convention
        (0 = liquid, 1 = ice Ih, 2 = ice II, ...).
        Shape is (nP, nT) when grid=True, (n,) when grid=False.
        """

    @abstractmethod
    def fn_Seismic(self, P_MPa, T_K, grid=False):
        """Return (VP_kms, VS_kms, KS_GPa, GS_GPa).
        For liquids VS_kms and GS_GPa are arrays of zeros.
        Each array has shape (nP, nT) when grid=True, (n,) when grid=False.
        """

    @abstractmethod
    def fn_sigma_Sm(self, P_MPa, T_K, grid=False):
        """Return electrical conductivity in S/m.
        Shape is (nP, nT) when grid=True, (n,) when grid=False.
        """

    @abstractmethod
    def fn_eta_Pas(self, P_MPa, T_K, grid=False):
        """Return dynamic viscosity in Pa·s.
        Shape is (nP, nT) when grid=True, (n,) when grid=False.
        """

    @abstractmethod
    def fn_species(self, P_MPa, T_K, grid=False, reactionEquation=None):
        """Return speciation information, or None if not modeled.
        Format must match that expected by any downstream PP species consumers.
        """

    def updateConvectionViscosity(self, etaConv_Pas, Tconv_K):
        """Update viscosity for convective layers. No-op by default; override if needed."""
        pass
    