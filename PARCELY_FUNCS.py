# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 11:34:16 2022

File containing all the functions needed to run PARCELY.

@author: Dan Barthaux
"""

import numpy as np
from pandas import read_csv
import numba as nb
from scipy.optimize import fmin, toms748, fsolve
from scipy.signal import savgol_filter
import warnings
from sklearn.cluster import KMeans


# =============================================================================
# Global Variables
# =============================================================================

g       = 9.81          # Gravitational acceleration, m/s^2
rhoL    = 997.0         # Density of liquid water, kg/m^3
Ma      = 28.97e-3      # Molar mass of dry air, kg/mol
Mw      = 18.01528e-3   # Molar Mass of water, kg/mol
R       = 8.3145        # Universal molar gas constant, J/K/mol
cp      = 1005.0        # Specific heat capacity (constant pressure), J/kg/K
Rv      = 467           # Gas constant for water vapour, J/kg/K
Rd      = 287.0         # Gas constant of dry air, J/kg/K
eta     = Rd/Rv
kB      = 1.380649e-23  # Boltzmann constant, J/K
NA      = 6.02214076e23 # Avogadro's number
Rd_alt  = 8.20573e-5    # Ideal gas constant [m3 atm K-1 mol-1]


# =============================================================================
# Dynamics
# =============================================================================

@nb.njit
def TerminalVel(Rad, Temperature, Pressure, AmbSat):
    """Function to determine the terminal velocity of a droplet as a function
    of its radius, using the analytic expression from force balance and an
    empirically derived drag coefficient function of radius.
    
    Rad         : m (Droplet radius)
    Temperature : K (Ambient temperature)
    Pressure    : Pa (Ambient pressure)
    AmbSat      : % (Saturation ratio)
    Output      : Array, m/s"""

    # Get dry air density
    rhoa = AirDensity(Temperature, Pressure, AmbSat)
    # Drag coefficient
    Cd = DragCoefficient(Rad)
    # Terminal velocity
    VelTerm = np.sqrt(8*Rad*g*(rhoL-rhoa)/(3*Cd*rhoa))
    
    return VelTerm


@nb.njit
def BoundaryConds3D(Drops, ZLIM): 
    """Function applying the periodic boundary conditions to the domain.
    
    Drops   : Droplet array
    ZLIM    : m (Array of the z-direction bounds)
    Output  : Droplet array."""
    
    # Copy the data so as not to overwrite
    DropTemp = Drops.copy()
    
    # Mask for if the position is larger/smaller than the bounds
    HitTopZBoundary = DropTemp[:,1] > ZLIM[1]
    HitBottomZBoundary = DropTemp[:,1] < ZLIM[0]
    
    # Shifting the positions to the opposite boundary
    DropTemp[HitTopZBoundary, 1] = ZLIM[0]
    DropTemp[HitBottomZBoundary, 1] = ZLIM[1]

    return DropTemp


# =============================================================================
# Physics
# =============================================================================

@nb.njit
def WaterSurfaceTension(Temperature, CritTemp=647.15):
    """Calculate the surface tension of water using the equation from Kalova
    and Mares (2018).
    
    Temperature : K (Ambient temperature of air)
    CritTemp    : K (Critical temperature of water)
    Output      : J/m2 (Surface tension)"""
    
    # Dimensionless variable
    tau = 1 - Temperature/CritTemp
    
    # Surface tension, J/m2
    sigma = 241.322*(tau**1.26)*(1- 0.0589*(tau**0.5) - 0.56917*tau)*1e-3
    
    return sigma
    

@nb.njit
def LatentHeatEvap(Temperature, T0=273.16, cpv=2040, L0=2.501e6):
    """Latent heat of condensation from Ambaum (2020).
    
    Temperature : K (Ambient temperature)
    T0          : K (Temperature at triple point)
    cpv         : J/kg/K (Isobaric mass heat capacity of vapor)
    L0          : J/kg (Latent heat of vaporization at triple point)
    Output      : J/kg (Latent heat of vaporization)."""
    
    # Latent heat of vaporization, J/kg
    L = L0 - (cpl-cpv)*(Temperature- T0)
    
    return L


@nb.njit
def LatentHeatSub(Temperature, T0=273.16, cpv=1885, L0=2.260e6):
    """Latent heat of condensation from Ambaum (2020).
    
    Temperature : K (Ambient temperature)
    T0          : K (Temperature at triple point)
    cpv         : J/kg/K (Isobaric mass heat capacity of vapor)
    L0          : J/kg (Latent heat of sublimation at triple point)
    Output      : J/kg (Latent heat of sublimation)."""
    
    # Latent heat of sublimation, J/kg
    L = L0 - (cpi-cpv)*(Temperature- T0)
    
    return L


@nb.njit
def SVP(Temperature, Phase, es0=611.655, T0=273.16, cpvl=2040, cpvi=1885,
        L0l = 2.501e6, L0i = 2.260e6):
    """Calculate the saturation vapour pressure for a given temperature and
    phase of water (ice or liquid). Taken from Ambaum (2020).
    
    Temperature : K (Ambient temperature)
    Phase       : String ('ice' or 'liquid')
    es0         : Pa (Saturation vapor pressure at triple point [TRPLP])
    T0          : K (Temperature at TRPLP)
    cpvl        : J/kg/K (Isobaric mass heat capacity of vapor [liquid SVP])
    cpvi        : J/kg/K (Isobaric mass heat capacity of vapor [ice SVP])
    L0l         : J/kg (Latent heat of vaporization at TRPLP)
    L0i         : J/kg (Latent heat of sublimation at TRPLP)
    Output      : Pa (Saturation vapor pressure)."""
    
    if Phase == 'liquid':
        
        # Latent heat of vaporization, J/kg
        L = LatentHeatEvap(Temperature, T0, cpvl, L0l)
        
        Part1 = es0*(T0/Temperature)**((cpl-cpvl)/Rv)
        Part2 = np.exp(L0l/(Rv*T0) - L/(Rv*Temperature))
        
        es = Part1*Part2
        
    if Phase == 'ice':
        
        # Latent heat of vaporization, J/kg
        L = LatentHeatSub(Temperature, T0, cpvi, L0i)
        
        Part1 = es0*(T0/Temperature)**((cpi-cpvi)/Rv)
        Part2 = np.exp(L0i/(Rv*T0) - L/(Rv*Temperature))
        
        es = Part1*Part2
        
    return es


@nb.njit
def Diffusivity(Temperature, Pressure, P0 = 101325, T0=273.15):
    """Calculate the coefficient of diffusivity given temperature T.
    Equation from Seinfeld and Pandis Chapter 17, equivalently also from
    Pruppacher and Klett (1997).
    
    Temperature : K (Ambient temperature)
    Pressure    : Pa (Ambient pressure)
    Output      : m^2/s (Diffusivity)."""
    
    # Diffusivity, m^2/s
    D = 1e-4*0.211*((Temperature/T0)**1.94)*(P0/Pressure)
    
    return D


@nb.njit
def ThermConduct(Temperature):
    """Calculate the thermal conductivity of air at a given temperature.
    Equation from Pruppacher and Klett or Seinfeld and Pandis.
    
    Temperature : K (Ambient temperature)
    Output      : J/m/s/K."""
    
    # Converts Kelvin to Celsius, J/m/s/K
    K = 4.1868e-3*(5.69 + 0.017*(Temperature - 273.15))
    
    return K


@nb.njit
def DragCoefficient(Radius):
    """Calculate the drag coefficient for a droplet with radius Rad using
    a fitted power-law function to data from Shafrir and Gal-Chen 1971.
    
    Radius : m
    Output : unitless."""
    
    # Coefficient with 95% confidence bounds
    a = 9.86e-10    #(7.956e-10, 1.176e-09)
    b = -2.375      #(-2.394, -2.357)
    c = 1.014       #(0.9607, 1.068)
    
    # Drag coefficient
    Cd = a*(Radius**b) + c
    
    return Cd


@nb.njit
def KineticEffects(Radius, Temperature, Pressure, AmbSat, ac, at):
    """Calculate thermal conductivity and diffusivity adjusted for kinetic
    effects using equations from Seinfeld and Pandis (2006).
    
    Radius      : m (Droplet radius)
    Temperature : K (Ambient temperature)
    Pressure    : Pa (Ambient pressure)
    AmbSat      : % (Saturation ratio)
    ac, at      : (coefficients used in the kinetic adjustments)
    Outputs     : (1) J/m/K/s (Conductivity), (2) m^2/s (Diffusivity)."""
     
    # Conductivity, diffusivity (J/m/s/K, m^2/s)
    K = ThermConduct(Temperature)
    D = Diffusivity(Temperature, Pressure)

    # Density of dry air, kg/m3
    rhoa = AirDensity(Temperature, Pressure, AmbSat)
    
    # Vapor and thermal jump lengths (m)
    Vv = 1.3*(8e-8)
    VT = 2.16e-7
        
    # Adjusted values of conductivity and diffusivity
    Kn = K/(Radius/(Radius+VT) + \
            (K/(at*Radius*rhoa*cp))*np.sqrt(2.0*np.pi/(Rd*Temperature)))
    
    Dn = D/(Radius/(Radius+Vv) + \
            (D/(ac*Radius))*np.sqrt(2.0*np.pi/(Rv*Temperature)))
    
    return Kn, Dn


@nb.njit
def AirDensity(Temperature, Pressure, AmbSat):
    """Function to calculate the density of air as a function of pressure and
    temperature, taking into account the partial pressure of water vapour using
    the environmental saturation.
    
    Temperature : K (Ambient temperature)
    Pressure    : Pa (Ambient pressure)
    AmbSat      : %, (Saturation ratio)
    Output      : kg/m^3 (Dry air density)."""
    
    # Saturation vapour pressure, Pa
    es = SVP(Temperature, 'liquid')
    # Given water vapour partial pressure, Pa
    e = es*AmbSat
    # Dry air density
    rhoa = (Pressure-e)/(Rd*Temperature)
    
    return rhoa


@nb.njit(cache=True)
def DropSurfTension(WetData, DryData, OAParameters, 
                    Temperature, method='OrgFilm', sft=0.072):
    """Function to calculate the effective surface tension of a droplet. Four
    modes are listen below in the 'method' parameter to calculate surface
    tension to varying complexity.
    
    WetData       : Droplet array
    DryData       : Solute array
    OAParameters  : SVO parameter array
    Temperature   : K (Ambient temperature)
    method        : string : ('WaterTemp', 'Constant', 'MixRule', 'OrgFilm')
    sft           : float, used for Constant surface tension method
    Output        : J/m^2 (Surface tension)"""
    
    # Surface tension of water as a function of temperature
    if method == 'WaterTemp':
        # J/m2
        AltSigma = np.ones_like(WetData[:,0])*WaterSurfaceTension(Temperature)
        
        return AltSigma
    
    # Constant surface tension, user defined
    if method == 'Constant':
        # J/m2
        AltSigma = sft*np.ones_like(WetData[:,0])
        
        return AltSigma
    
    # Volumes from mass and density
    Vos = DryData[:,6:]/OAParameters[:,1] # Organics
    Vw  = WetData[:,3]/rhoL               # Water
    Vd  = DryData[:,2]/DryData[:,4]       # Solutes
    
    # Surface tension of water, J/m2
    SigWater = WaterSurfaceTension(Temperature)
    
    # Volume-weighted mixing rule
    if method == 'MixRule':
        # Check if any organics even exist in droplets
        if np.all(DryData[:,6:]) != 0:
            
            # Calculating for every droplet
            if DryData[:,6:].ndim > 1:
                # Total volume, m^3
                Vtot = Vw + Vos.sum(axis=1) + Vd
                # Volume-weighted surface tension contributions
                VSigW  = SigWater*Vw/Vtot
                VSigD  = SigWater*Vd/Vtot
                VSigOA = np.zeros_like(DryData[:,6:])
                # Get contribution for each OA
                for i in range(DryData[:,6:].shape[1]):
                    VSigOA[:,i] = OAParameters[i,6]*Vos[:,i]/Vtot
                        
                # Mixture surface tension, J/m2
                SigMix = VSigW + VSigOA.sum(axis=1) + VSigD
                
            # Calculating only for one droplet
            if DryData[:,6:].ndim == 1:
                # Total volume
                Vtot = Vw + Vos.sum() + Vd
                # Volume-weighted surface tension contributions
                VSigW = SigWater*Vw/Vtot
                VSigD  = SigWater*Vd/Vtot
                VSigOA = np.zeros((Vw.size, DryData[:,6:].size))
                # Get contribution for each OA
                for i in range(DryData[:,6:].size):
                    VSigOA[:,i] = OAParameters[i,6]*Vos[i]/Vtot
        
                # Mixture surface tension
                SigMix = VSigW + VSigOA.sum(axis=1) + VSigD
        
        else:
            SigMix = np.ones(Vw.size)*SigWater
        
        return SigMix
    
    # Organic-film method
    if method == 'OrgFilm':
        
        # Characteristic monolayer thickness taken from AIOMFAC
        delz = 0.2e-9
        # Volume of the monolayer, m^3
        Vmono = (4/3)*np.pi*(WetData[:,4]**3 - (WetData[:,4] - delz)**3)
        
        # Determine if there's enough organics to cover the particle
        # If not enough, define the surface coverage fraction, else = 1
        SurfCov = np.where(Vos.sum(axis=1) <= Vmono, Vos.sum(axis=1)/Vmono, 1)
        
        # Get volume weighted-mean surface tension
        VSigOA = np.zeros_like(SurfCov)
        for i in range(Vos.shape[0]):
            VSigOA[i] = np.sum(Vos[i,:]*OAParameters[:,6]/Vos[i,:].sum())
        
        # Effective surface tension, J/m2
        EffSigma = VSigOA*SurfCov + (1-SurfCov)*SigWater
        
        return EffSigma
    

@nb.njit
def KappaKoehler(WetData, SolData, OAParameters, Temperature, SurfMode, sft):
    """Function to calculate the saturation ratio of water vapor at the 
    surface of a droplet using Kappa-Koehler theory by 
    Petters and Kreidenweis (2007).
    
    WetData      : Droplet data array
    SolData      : Solute data array
    OAParameters : Organic aerosol parameter array
    Temperature  : K, (Ambient temperature)
    SurfMode     : string, method of calculating surface tension
    sft          : if method is 'Constant', value for surface tension
    Output       : Unitless, J/m^2."""
    
    # Kelvin factor, unitless
    Factor = (WetData[:,4]**3 - SolData[:,1]**3)/(WetData[:,4]**3 - \
                (SolData[:,1]**3)*(1.0- SolData[:,5]))
    
    # Check to see if there's any organics
    Check = np.all(SolData[:,6:] == 0)
    # If there are organics or a constant surface tension is requested
    if (SolData.shape[1] > 6 and Check == False) or SurfMode=='Constant':
        # Calculate surface tension of the droplet, J/m2
        Sigsw = DropSurfTension(WetData, SolData, OAParameters, 
                                Temperature, SurfMode, sft)
            
        # Droplet surface saturation ratio
        DropSat = Factor*np.exp((2.0*Sigsw*Mw)/\
                                (Temperature*R*WetData[:,4]*rhoL))
    
    else:
        # Calculate surface tension of the droplet, J/m2
        Sigsw = np.ones_like(Factor)*WaterSurfaceTension(Temperature)
        # Droplet surface saturation ratio
        DropSat = Factor*np.exp((2.0*Sigsw*Mw)/\
                                (Temperature*R*WetData[:,4]*rhoL))
    
    return DropSat, Sigsw


@nb.njit
def KapKrit(Radius, DryRadius, A, Kappa):
    """Wrapper function to find the critical radius and supersaturation of
    a droplet, using in CritParamFinder. Returns a negative supersaturation
    in order to use the fmin function.
    
    Radius      : m (Wet radius)
    DryRadius   : m (Solute radius)
    A           : (Collection of parameters in the surface tension term)
    Kappa       : unitless (Hygroscopicity)
    Output      : -% (Negative supersaturation)"""
    
    # Rename variables
    r = Radius
    rd = DryRadius
    k = Kappa
    
    # Supersaturation
    Ko = (((r**3 - rd**3)/(r**3 - rd**3*(1-k))*np.exp(A/r)) -1)*100
    
    # Return negative value
    return -Ko


@nb.njit
def KapAmbient(Radius, DryRadius, A, Kappa, AmbSat):
    """Wrapper function to find the equilibrium radius of a droplet at a 
    given saturation, used in CritParamFinder. 
    
    Radius      : m (Wet radius)
    DryRadius   : m (Solute radius)
    A           : (Collection of parameters in the surface tension term)
    Kappa       : unitless (Hygroscopicity)
    AmbSat      : % (Ambient saturation ratio)
    Output      : (Inverted supersaturation)"""
    
    # Rename variables
    r = Radius
    rd = DryRadius
    k = Kappa
    S = AmbSat
    
    # Droplet surface saturation ratio minus the ambient saturation ratio
    Ko = ((r**3 - rd**3)/(r**3 - rd**3*(1-k))*np.exp(A/r)) - S
    
    return Ko


@nb.njit
def GrowthRate(WetData, AmbSat, Temperature, Pressure, ac, at):
    """Function that calculates the new droplet radius and the mass growth 
    rate of a droplet(s).
    
    WetData     : Water data array
    AmbSat      : % (Saturation ratio)
    Temperature : K (Ambient temperature)
    Pressure    : Pa (Ambient pressure)
    Output      : kg/s (Growth rates)"""

    # Saturation vapor pressure, Pa
    es = SVP(Temperature, 'liquid')
    # Conductivity (J/m/s/K) and Diffusivity (m^2/s)
    K, D = KineticEffects(WetData[:,4], Temperature, Pressure, AmbSat, ac, at)
    # Latent head of condensation, J/kg
    L = LatentHeatEvap(Temperature)
    
    # Growth parameters, mass-based
    Fk = (L/(Rv*Temperature)-1.0)*L/(K*Temperature)
    Fd = Rv*Temperature/(D*es)
    
    # Growth rate
    GRate = (AmbSat-WetData[:,5])*((4*np.pi)**(2/3))*\
            ((3/rhoL)**(1/3))*(WetData[:,3]**(1/3))/(Fd + Fk)
    
    return GRate


@nb.njit
def NewMass(WetData, DryData, AmbSat, Temperature, Pressure,
            OAParameters, ac, at, dt, mtol, SurfMode, sft):
    """Runge-Kutta Cash-Karp 4th order method to solve the growth equation from
    the GrowthRate function, returning a new array of water masses for each
    droplet.
    
    WetData         : Water data array
    DryData         : Solute data array
    AmbSat          : % (Saturation ratio)
    Temperature     : K (Ambient temperature)
    Pressure        : Pa (Ambient pressure)
    OAParameters    : Organics data array
    ac, at          : Kinetic coefficients
    dt              : s (time-step size)
    mtol            : kg (Error tolerance)
    SurfMode        : string (Surface tension method)
    sft             : if SurfMode is 'Constant', value for it in J/m^2
    Output          : kg (New mass)."""
    
    # Calculate k1
    k1 = WetData[:,6]
    
    # Adjust solution
    StepData = WetData.copy()
    
    StepData[:,3] = WetData[:,3] + k1*dt/5
    StepData[:,4] = ((3/(4*np.pi)*(StepData[:,3]/rhoL)) +\
                      DryData[:,1]**3)**(1/3) 
    StepData[:,5], StepData[:,7] = KappaKoehler(StepData, DryData, 
                                    OAParameters, Temperature, SurfMode, sft)
    
    # Calculate k2
    k2 = GrowthRate(StepData, AmbSat, Temperature, Pressure, ac, at)
    
    # Adjust solution
    StepData[:,3] = WetData[:,3] + dt*(k1*(3/40) + k2*(9/40))
    StepData[:,4] = ((3/(4*np.pi)*(StepData[:,3]/rhoL)) + \
                      DryData[:,1]**3)**(1/3)
    StepData[:,5], StepData[:,7] = KappaKoehler(StepData, DryData, 
                                    OAParameters, Temperature, SurfMode, sft)
    
    # Calculate k3
    k3 = GrowthRate(StepData, AmbSat, Temperature, Pressure, ac, at)
    
    # Adjust solution
    StepData[:,3] = WetData[:,3] + dt*(k1*(3/10) - k2*(9/10) + k3*(6/5))
    StepData[:,4] = ((3/(4*np.pi)*(StepData[:,3]/rhoL)) + \
                      DryData[:,1]**3)**(1/3)
    StepData[:,5], StepData[:,7] = KappaKoehler(StepData, DryData, 
                                    OAParameters, Temperature, SurfMode, sft)
    
    # Calculate k4
    k4 = GrowthRate(StepData, AmbSat, Temperature, Pressure, ac, at)
    
    # Adjust solution
    StepData[:,3] = WetData[:,3] - dt*(k1*(11/54) + k2*(5/2) - k3*(70/27)\
                    + k4*(35/27))
    StepData[:,4] = ((3/(4*np.pi)*(StepData[:,3]/rhoL)) + \
                      DryData[:,1]**3)**(1/3)
    StepData[:,5], StepData[:,7] = KappaKoehler(StepData, DryData, 
                                    OAParameters, Temperature, SurfMode, sft)
    
    # Calculate k5
    k5 = GrowthRate(StepData, AmbSat, Temperature, Pressure, ac, at)
    
    # Adjust solution
    StepData[:,3] = WetData[:,3] + dt*(k1*(1631/55296) + k2*(175/512) - \
                    k3*(575/13824) + k4*(44275/110592) + k5*(253/4096))
    StepData[:,4] = ((3/(4*np.pi)*(StepData[:,3]/rhoL)) + \
                      DryData[:,1]**3)**(1/3)
    StepData[:,5], StepData[:,7] = KappaKoehler(StepData, DryData, 
                                    OAParameters, Temperature, SurfMode, sft)
    
    # Calculate k6
    k6 = GrowthRate(StepData, AmbSat, Temperature, Pressure, ac, at)
    # Mass 4th order accurate
    MassE4 = WetData[:,3] + dt*(((37/378)*k1) + ((250/621)*k3) + \
                                ((125/594)*k4) + ((512/1771)*k6))
    # Mass 5th order accurate
    MassE5 = WetData[:,3] + dt*(((2825/27648)*k1) + ((18575/48384)*k3) + \
            ((13525/55296)*k4)+((277/14336)*k5)+((1/4)*k6))
    
    # Take the largest error from the population of particles
    Error = np.max(np.abs(MassE4-MassE5))
    # Calculate new time-step size, s
    if Error != 0:
        Newdt = 0.8*dt*(mtol/Error)**(0.2)
    else:
        Newdt = 2*dt
    
    return MassE4, Newdt, Error


@nb.njit
def EnvSatTime(Temperature, Pressure, AmbSat, Updraft, GrowthRates, dt, Cubes):
    """Equation for the change in environmental saturation ratio, dependent on
    temperature, pressure, updraft velocity, the sum of the droplet growth 
    rates, the current saturation ratio and the time-step size.
    
    Temperature : K (Ambient temperature)
    Pressure    : Pa (Ambient pressure)
    AmbSat      : % (Saturation ratio)
    Updraft     : m/s (Updraft velocity)
    GrowthRates : kg/s (Droplet growth rates)
    dt          : s (Time-step size)
    Cubes       : unitless (number of cubic centimeters)
    Output      : % (New saturation ratio)."""
     
    # Saturation vapour pressure, Pa
    es = SVP(Temperature, 'liquid')
    # Latent heat, J/kg
    L = LatentHeatEvap(Temperature)
    # Density of dry air, kg/m3
    rhoa = AirDensity(Temperature, Pressure, AmbSat)
    
    # Parameters
    # For adiabatic cooling
    Q1 = (1/Temperature)*((L*g)/(Rv*Temperature*cp)) - g/(Rd*Temperature)
    # For condensation onto droplets
    Q2 = Pressure/(eta*es) + (1/Temperature)*(L**2/(Rv*Temperature*cp))

    # New saturation ratio
    Snew = AmbSat + dt*(Q1*Updraft - \
                    (Q2*GrowthRates.sum())/(rhoa*Cubes*(1e-2)**3))
    
    return Snew


@nb.njit
def VirtualTemp(Temperature, Pressure, AmbSat):
    """Calculate the virtual temperature.
    
    Temperature : K (Ambient temperature)
    Pressure    : Pa (Ambient pressure)
    AmbSat      : % (Saturation ratio)
    Output      : K (Virtual temperature."""
    
    # Saturation vapor pressure, Pa
    es = SVP(Temperature, 'liquid')
    # Vapor pressure using saturation ratio, Pa
    e = es*AmbSat
    # Mixing ratio
    w = eta*e/Pressure
    # Virtual temperature, K
    Tv = Temperature*((1 + w/eta)/(1+w))
    
    return Tv


@nb.njit
def PressEvolution(Temperature, Pressure, AmbSat, Updraft, dt):
    """Function to calculate the new pressure with time using the 
    hypsometric equation.
    
    Temperature : K (Ambient temperature)
    Pressure    : Pa (Ambient pressure)
    AmbSat      : % (Saturation ratio)
    Updraft     : m/s (Updraft velocity)
    dt          : s (TIme-step size)
    Output      : Pa (New ambient pressure)."""
    
    # Virtual temperature, K
    Tv = VirtualTemp(Temperature, Pressure, AmbSat)
    # New pressure, Pa
    Pnew = Pressure*np.exp(-g*Updraft*dt/(Rd*Tv))
    
    return Pnew


@nb.njit
def TempEvolution(Temperature, Pressure, AmbSat, Updraft,
                  GrowthRates, dt, Cubes):
    """Function to calculate the new temperature with time.
    
    Temperature : K (Ambient temperature)
    Pressure    : Pa (Ambient pressure)
    AmbSat      : % (Saturation ratio)
    Updraft     : m/s (Updraft velocity)
    dt          : s (Time-step size)
    Cubes       : unitless (Number of cubic centimeters)
    Output      : K (New ambient temperature)."""
    
    # Latent heat of condensation, J/kg
    L = LatentHeatEvap(Temperature)
    # Dry air density, kg/m3
    rhoa = AirDensity(Temperature, Pressure, AmbSat)
    # New environmental temperature, K
    Tnew = Temperature - \
        (g*Updraft/cp - (L/(cp*rhoa))*GrowthRates.sum()/(Cubes*(1e-2)**3))*dt
    
    return Tnew


# =============================================================================
# Solute Population
# =============================================================================

def SolPopInitialize(Ns, Rs, Stds, Inorgs, InorgPopFrac, Orgs, PercCond,
                        OAParameters, RNG, MassFractions, Dirichlet, 
                        KappaRand, Organics, Cubes, Distribution):
    """Function to create the solute array for the droplet population.
    
    Ns              : # (Number of particles)
    Rs              : um (Mean radii)
    Stds            : um (Standard deviation)
    Inorgs          : Array of lists of indices (selected inorganic species)
    InorgPopFrac    : Array of lists (Population percentage of each inorganic)
    Orgs            : List of lists (condensed organics in each 
                                                         sub-distribution)
    PercCond        : Integer or array (condensed percentage of each organic
                                                    in each sub-distribution)
    OAParameters    : Organics parameter data array
    RNG             : Numpy random generator, object
    MassFractions   : Array of mass fractions
    Dirichlet       : Bool Trigger, randomizes mass fractions using a Dirichlet
                                                distribution for every particle
    KappaRand       : Bool Trigger, sets hygroscopicity values around a normal
                                        distribution of the pure (input) value
    Organics        : Bool Trigger, whether organics are included in the run
    Cubes           : unitless (Number of cubic centimeters)
    Distribution    : string (Aerosol distribution: 'mono', 'normal', 
                                                          'lognormal')."""
    
    # Solute Parameters
    # 0: Tracking number
    # 1: Radius (m)
    # 2: Mass (kg)
    # 3: Molar Mass (kg/mol)
    # 4: Density (kg/m^3)
    # 5: Kappa (CCN)
    # 6-on: Condensed OA Masses (kg)
    
    # Total number of particles
    N = np.sum(Ns)

    # Solute parameters data table
    SoluteData = read_csv('CCN_List.csv', delimiter=',')
    
    # Radius distribution
    PopModes = RadiusDistributor(Ns, Rs, Stds, RNG, Distribution)
    
    # Holding lists of the lognormal components for later concatenating
    Components = []
    FiltComponents = []
    
    # If PercCond is an array, calculate mass fractions and vapor content
    # of organics alternatively
    if type(PercCond) != int:
        MassFractions, Vapor = MassFracCalc(OAParameters, PercCond, Cubes, 
                                        Ns, Rs, Stds, Inorgs, InorgPopFrac,
                                        RNG, Distribution, SoluteData)
    
    # Loop through each lognormal contribution
    for m in range(Ns.size):
        
        # Add organics columns if True
        if Organics is True:
            Mode = np.zeros(shape=(Ns[m], 6+OAParameters.shape[0]))
        if Organics is False:
            Mode = np.zeros(shape=(Ns[m], 6))
        
        # Get dry radii, m
        Mode[:,1] = PopModes[m]
        # Set any radii smaller than 5 nm to 5 nm
        Mode[:,1] = np.where(Mode[:,1] < 5e-9, 5e-9, Mode[:,1])
        
        # Allocate according to user-defined population fractions
        AllocInOrg = SoluteAllocator(Inorgs[m], InorgPopFrac[m], Ns[m])

        # Assign correct CCN parameters
        for n in range(Ns[m]):
            for CCN in Inorgs[m]:
                if AllocInOrg[n] == CCN:
                    # Assign attributes according to assigned number
                    if KappaRand is True:
                        Mode[n,3:5] = SoluteData.iloc[CCN,2:4]
                        Mode[n,5] = KappaRandomizer(SoluteData.iloc[CCN,4],
                                                    RNG, 'Normal')
                    else:
                        Mode[n,3:6] = SoluteData.iloc[CCN,2:]
                    # Calculate mass by assuming spherical CCN
                    if Organics is False or \
                        (Organics is True and \
                          (len(MassFractions[m]) == len(Inorgs[m])) is True):
                            Mode[n,2] = (4*np.pi/3)*SoluteData.iloc[CCN,3]\
                                *(Mode[n,1]**3)

        # Convert molar mass from g/mol to kg/mol
        Mode[:,3] = Mode[:,3]*1e-3

        if Organics is True:
            # Indexes for ease later
            RelOrgs = np.array(Orgs[m])+6
            
            # If a non-monotonic mass fraction distribution is wanted for the
            # particles in the domain, create a Dirichlet distribution
            if Dirichlet is True:
                MassFractions[m] = RNG.dirichlet(MassFractions[m]*100, Ns[m])
            
            # Calculate the mean density of each droplet accounting for
            # chosen condensed organics, kg/m3
            MeanDensity = np.zeros(Ns[m])
            for n in range(Ns[m]):
                if len(MassFractions[m].shape) > 1:
                    MeanDensity[n] = np.dot(MassFractions[m][:,n], 
                                        np.append(Mode[n,4], 
                                              OAParameters[Orgs[m],1]))
                else:
                    MeanDensity[n] = np.dot(MassFractions[m], 
                                        np.append(Mode[n,4], 
                                              OAParameters[Orgs[m],1]))
            
            # Calculate the mean mass of each droplet, kg
            Mode[:,2] = MeanDensity*(Mode[:,1]**3)*(4/3)*np.pi
        
            # Calculate each component's mass from mass fractions and mean mass
            for n in range(Ns[m]):
                if len(MassFractions[m].shape) > 1:
                    Mode[n,RelOrgs] = Mode[n,2]*MassFractions[m][1:,n]
                else:
                    Mode[n,RelOrgs] = Mode[n,2]*MassFractions[m][1:]
                    
            # Create array of pure inorganic CCN for future use
            ModeFilter = Mode[:,:6].copy()
            if len(MassFractions[m].shape) > 1:
                ModeFilter[:,2] = Mode[:,2]*MassFractions[m][0,:]
                
            else:
                ModeFilter[:,2] = Mode[:,2]*MassFractions[m][0]
                
            # Get radii, m
            ModeFilter[:,1] = ((ModeFilter[:,2]/Mode[:,4])*\
                                  (3/(4*np.pi)))**(1/3)
                
            # Find effective hygroscopicity
            Mode[:,5] = KappaChanger(ModeFilter[:,2], Mode[:,RelOrgs], 
                                      ModeFilter[:,4], 
                                      OAParameters[Orgs[m],1], 
                                      ModeFilter[:,5], 
                                      OAParameters[Orgs[m],5])
                
            # Adjust densities to mean densities
            Mode[:,4] = MeanDensity
        
        Components.append(Mode)
        if Organics is True:
            FiltComponents.append(ModeFilter)
    
    # Put all sub-distributions together
    Solutes = np.concatenate(Components)
    # Tracking IDs
    Solutes[:,0] = np.arange(0, N, 1, dtype=int)
    
    if Organics is True:
        SoluteFilter = np.concatenate(FiltComponents)
        SoluteFilter[:,0] = np.arange(0, N, 1, dtype=int)
        # Total condensed organics in ug/m3
        TotalCondensedOrg = (Solutes[:,6:].sum(axis=0))*1e9/(Cubes*(1e-2)**3)
        # How much is left in the gas phase
        VapConc = OAParameters[:,7] - TotalCondensedOrg
        VaporConcentrations = np.where(VapConc>0, VapConc, 0)
        OAParameters[:,7] = VaporConcentrations + TotalCondensedOrg

        return Solutes, SoluteFilter, VaporConcentrations        

    if Organics is False:
        SoluteFilter = Solutes.copy()
        return Solutes, SoluteFilter, 0


def MassFracCalc(OAParameters, PercCond, Cubes, Ns, Rs, Stds, Inorgs,
                InorgPopFrac, RNG, Distribution, SoluteData):
    """Function to calculate the solute mass fractions of the organic species
    if PercCond is given, meaning the condensed percentage of each organic
    is requested. Does not consistently return good results.
    
    OAParameters    : Organics parameter data array
    PercCond        : Array (condensed percentage of each organic in each 
                                                             sub-distribution)
    Cubes           : unitless (Number of cubic centimeters)
    Ns              : # (Number of particles)
    Rs              : um (Mean radii)
    Stds            : um (Standard deviation)
    Inorgs          : Array of lists of indices (selected inorganic species)
    InorgPopFrac    : Array of lists (Population percentage of each inorganic)
    RNG             : Numpy random generator, object
    Distribution    : string (Aerosol distribution: 'mono', 'normal', 
                                                              'lognormal')
    SoluteData      : CSV file data on Inoranic parameters."""
    
    # Get OA total concentrations into total contents, kg
    Concentrations = OAParameters[:,7].copy()*1e-15*Cubes
    # Get condensed and vapor totals, kg
    InParticles = Concentrations*PercCond
    Vapor = Concentrations*(1-PercCond)
    
    # Radius distribution
    Radii = np.concatenate(RadiusDistributor(Ns, Rs, Stds, RNG, Distribution))
    # Use radii as a size-weighting
    Weight = Radii/(Radii.sum())
    # Assign mass of OA to each particle
    OAPerParticle = np.zeros((OAParameters.shape[0], Radii.size))
    for x in range(OAParameters.shape[0]):
        OAPerParticle[x,:] = Weight*InParticles[x]
    
    # Allocate inorganic density
    SaltDens = []
    for m in range(Ns.size):
        AllocInOrg = SoluteAllocator(Inorgs[m], InorgPopFrac[m], Ns[m])
        Densities = np.zeros_like(AllocInOrg)
        
        for n in range(Ns[m]):
            for CCN in Inorgs[m]:
                if AllocInOrg[n] == CCN:
                    # Assign density according to assigned number
                    Densities[n] = SoluteData.iloc[CCN,3]
        
        SaltDens.append(Densities)
    
    # Inorganic density, kg/m3
    SaltRho = np.concatenate(SaltDens)
    
    # Shorthand for wrapper function
    a = (3/(4*np.pi*Radii**3))
    b = np.zeros_like(a)
    c = OAPerParticle.sum(axis=0)
    # Assuming the particle is wholly inorganic
    PureSalt = (4/3)*np.pi*((Radii)**3)*SaltRho
    # Initial guess for inorganic mass
    InitGuess = PureSalt - c
    # Populate shorthand b as sum mass*density of organic
    for n in range(Radii.size):
        b[n] = np.dot(OAParameters[:,1], OAPerParticle[:,n])
    
    # Numerically solve for total mass
    Mtot = np.zeros(Radii.size)
    for i in range(Radii.size):
        Mtot[i] = fsolve(WrapMassFrac, InitGuess[i], 
                               args=(a, b, c, SaltRho[i], i))
    
    # Get inorganic mass
    Msalt = Mtot - c

    # Calculate mass fractions
    MassFractions = np.zeros((OAParameters.shape[0]+1, Radii.size))
    for n in range(Radii.size):
        MassFractions[0,n] = Msalt[n]/Mtot[n]
        MassFractions[1:,n] = OAPerParticle[:,n]/Mtot[n]   
        
    # Redivide MassFractions for each mode
    AllMassFracs = []
    for m in range(Ns.size):
        if m == 0:
            Divider = MassFractions[:,:Ns[m]]
            AllMassFracs.append(Divider)
            
        if (m > 0 and m < Ns.size-1):
            Divider = MassFractions[:,Ns[:m].sum():Ns[:m+1].sum()]
            AllMassFracs.append(Divider)
            
        if m == Ns.size-1:
            Divider = MassFractions[:,Ns[:m].sum():]
            AllMassFracs.append(Divider)
    
    return AllMassFracs, Vapor


def WrapMassFrac(x, a, b, c, rho, k):
    """Wrapper function to numerically resolve amount of each solute component
    in an aerosol to fit the given radius."""
    
    return a[k]*x**2 - rho*x - b[k] + rho*c[k]


def KappaRandomizer(Kappa, Generator, Distribution):
    """Function to select a hygroscopicity value from a normal distribution
    around the pure component value.
    
    Kappa : unitless (Hygroscopicity)
    Generator : Numpy random generator object
    Distribution : String (Aerosol distirbution: 'mono', 'normal', 'lognormal')
    Output : Unitless (New hygroscopicity)."""
    
    if Distribution == 'Normal':
        NewKappa = Generator.normal(Kappa, 0.1, 1)
    
    return NewKappa


def CritParamFinder(Solutes, Temperature, AmbSat):
    """Function to calculate the critical parameters of activation of each
    droplet and determine the equilibrium radius for the initial saturation.
    This function assumes the surface tension of water even if there are
    organics present. This makes the calculation massively simpler.
    
    Solutes     : Data array containing each CCN's dry parameters.
    Temperature : K (Ambient temperature)
    AmbSat      : % (Saturation ratio)
    Output      : Equilibrium radius [m], critical parameters [m, %]."""
    
    # Number of droplets
    DropCount = Solutes.shape[0]
    
    # Surface tension of water, J/m2
    SigWater = WaterSurfaceTension(Temperature)
    
    # Collect surface tension parameters
    A = 2*SigWater*Mw/(R*Temperature*rhoL)

    # Critical parameters array, wet radius and saturation ratio
    CritParams = np.zeros(shape=(DropCount, 2))
    
    # Array to find the radii at the initial environmental saturation ratio
    EquiR = np.zeros(DropCount)
    
    for n in range(DropCount):
        # Find the critical supersaturation and critical radius
        CritParams[n,:] = fmin(KapKrit, Solutes[n,1]*1.0000001, 
                               args=(Solutes[n,1], A, Solutes[n,5]), 
                               full_output=True, disp=0)[:2]
        
        # Using the critical radius as an upper limit, find the equilibrium
        # radius for the given ambient saturation ratio.
        EquiR[n] = toms748(KapAmbient, Solutes[n,1]*1.0001, 
                                 CritParams[n,0], args=(Solutes[n,1], A, 
                                                       Solutes[n,5], AmbSat))
    
    # Reformat critical supersaturation into saturation ratio
    CritParams[:,1] = -1*CritParams[:,1]/100 + 1
    
    return EquiR, CritParams


def RadiusDistributor(Pops, Means, Stds, RNG, DistKind):
    """Function to create a distribution of CCN radii for model initialization.
    
    Pops    : Number count of each mode. Sum must be equal to total number
              of droplets specified
    Means   : um (Underlying means to be translated to geometric means)
    Stds    : um (Underlying standard deviations to be translated to geometric)
    RNG     : Numpy random generator, object
    Output  : List of arrays [m]."""
    
    # List to hold the underlying distributions
    Distributions = []
    
    if DistKind == 'mono':
        for i in range(len(Pops)):
            # Get particle count, radius of each mono-population
            Pop, mu = Pops[i], Means[i]*1e-6
            # Get population
            X = np.ones(Pop)*mu
            Distributions.append(X)
   
    if DistKind == 'normal':
       for i in range(len(Pops)):
           # Get population, mean, and standard deviation
           Pop, mu, sigma = Pops[i], Means[i]*1e-6, Stds[i]*1e-6
           # Get distribution
           X = RNG.normal(mu, sigma, Pop)
           Distributions.append(X) 
           
    if DistKind == 'lognormal':
        for i in range(len(Pops)):
            # Get population, mean, and standard deviation
            Pop, mu, sigma = Pops[i], Means[i]*1e-6, Stds[i]
            # Get distribution
            X = RNG.lognormal(np.log(mu), np.log(sigma), Pop)
            Distributions.append(X)
            
    # Return list
    return Distributions


def SoluteAllocator(SoluteNames, Proportions, NumberDrops):
    """Function to populate the CCN array with the solutes in the list
    SoluteNames and their respective count proportions (array of percentages).
    
    SoluteNames : list/array of row indices
    Proportions : List or array of solute proportions
    NumberDrops : integer
    Output      : Data array."""
    
    # Renaming input variables for ease
    S = SoluteNames
    N = NumberDrops
    P = (N*Proportions).astype(int)
    
    if P.sum() != N:
        if P.sum() > N:
            Extra = P.sum() - N
            counter = 0
            while Extra != 0:
                P[counter] += -1
                counter += 1
                Extra += -1
                if counter == P.size:
                    counter = 0
        
        if P.sum() < N:
            Remaining = N - P.sum()
            counter = 0
            while Remaining != 0:
                P[counter] += 1
                counter += 1
                Remaining += -1
                if counter == P.size:
                    counter = 0
    
    # Assign amount of drops to each CCN
    SolutePop = np.zeros(N)
    
    for i in range(S.shape[0]):
        if i == 0:
            SolutePop[:P[i]] = S[i]
            
        if (i > 0 and i < S.shape[0]-1):
            SolutePop[P[:i].sum():P[:i+1].sum()] = S[i]
            
        if i == S.shape[0]-1:
            SolutePop[P[:i].sum():] = S[i]

    return SolutePop


# =============================================================================
# Droplet Population
# =============================================================================

def ActivateDetectVisNew(DropData, Offset, Method='argmax'):
    """Alternate activation detection function. Instead of detecting activation
    from a threshold critical radius, find minimum in time series of radius
    over mass (both of water).
    
    DropData : Wet data array
    OffSet   : If there is numerical noise in the beginning, offset start by
                                                       x amount of time-steps.
    Method   : Either 'argmax' or 'kmeans' in detecting which droplets 
                                                                activated.
    Output :
    (1) Array of lists for every time-step containing which droplets
        activated at that time-step.
    (2) Number of activated droplets over time
    (3) Derivative for every droplet for analysis if suspicious results.
    (4) Points of activation on the derivative, also for analysis.
    (5) Indices of activated droplets.
    (6) Indices of inactivated droplets."""
    
    # Get radii array
    Radii = DropData[:,4,:]
    # Smooth radii
    RadFilt = savgol_filter(Radii, 51, 5)
    
    if Method == 'kmeans':
        
        # Get final time-step of every droplet's radius
        End = DropData[:,4,-1].reshape(-1,1)
        kmeans = KMeans(init='k-means++' ,n_clusters=2, n_init=3, max_iter=300)
    
        # Suppress warning that triggers when no droplets activate
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Using k-means clustering, identify indices of droplets that 
            # belong to either activated or inactivated
            Cluster = kmeans.fit(np.log10(End)).labels_
    
        # Cluster groups
        Group1 = np.argwhere(Cluster==0).flatten()
        Group2 = np.argwhere(Cluster==1).flatten()
        
        # If one group empty, set populated as activated
        if Group1.size == 0:
            Activated = Group2
            UnActivated = Group1
        if Group2.size == 0:
            Activated = Group1
            UnActivated = Group2
        
        # If two groups, assign larger mean to activated
        if Group1.size > 0 and Group2.size > 0:
            if DropData[Group1,4,-1].mean() > DropData[Group2,4,-1].mean():
                Activated = Group1
                UnActivated = Group2
            else:
                Activated = Group2
                UnActivated = Group1
    
    if Method == 'argmax':
        Activated = []
        UnActivated = []
        for d in range(Radii.shape[0]):
            if np.argmax(RadFilt[d]) == RadFilt.shape[1]-1:
                Activated.append(d)
            else:
                UnActivated.append(d)
    
        Activated = np.array(Activated)
        UnActivated = np.array(UnActivated)
    
    # Get log10 of second derivative of radius,
    # ignore first Offset steps due to equilibrium noise
    Deriv = np.log10(np.abs(np.diff(np.diff(RadFilt[:,Offset:]))))
    # Array to store index of activation for each droplet
    IndLoc = np.zeros(Deriv.shape[0])

    for d in range(Radii.shape[0]):
        if d in Activated:
            # Find maximum peaks in the 2nd-derivative
            IndLoc[d] = np.argmax(Deriv[d]) + Offset
        else:
            IndLoc[d] = -1

    # Two arrays for activation over time. First as a function of time,
    # second as a function of droplet index. Functionally same data.
    TimeAct = []
    Peaks = [[] for d in range(Radii.shape[0])]
    # Easy form of activated count over time for plotting and analysis
    CumulCount = np.zeros(Radii.shape[1])

    for t in range(Radii.shape[1]):
        Holder = []
        for d in range(Radii.shape[0]):
            if IndLoc[d] == t:
                if d in Activated:
                    Peaks[d].append(t-Offset)
                    Holder.append(d)
        TimeAct.append(Holder)
        CumulCount[t] = len(TimeAct[t])

    # Cumulative summation
    CumulCount = np.cumsum(CumulCount)

    # Convert Peaks sublists to arrays for ease
    for d in range(Radii.shape[0]):
        Peaks[d] = np.array(Peaks[d])
    
    return TimeAct, CumulCount, Deriv, Peaks, Activated, UnActivated


def ActivationDetection(DropData, CritData):
    """Function to detect which droplets in a population activate at which
    point in time using the critical parameters.
    
    DropData : Model output 3D array.
    CritData : 1D array of the critical parameters for each droplet.
    Output: 
    (1) Number of activated droplets over time.
    (2) Indices of activated droplets.
    (3) Indices of unactivated droplets."""
    
    # Get wet radii
    RadWet = DropData[:,8,:]
    # Get critical radii
    CritRad = CritData[:,0]
    
    # Lists for tracking numbers
    Activated = []
    
    for t in range(DropData.shape[2]):
        # Where the radius is larger than the critical radius
        Filter = RadWet[:,t] >= CritRad
        # Get the IDs of the activated droplets
        DropIDs = DropData[Filter,0,0].astype(int)
        Activated.append(DropIDs)
        
    # Array of droplets that never grow
    NeverActivated = np.setdiff1d(DropData[:,0,0], Activated[-1]).astype(int)
    
    Cumulative = np.zeros(len(Activated))

    for t in range(Cumulative.size):
        Cumulative[t] = Activated[t].size

    return Cumulative, Activated, NeverActivated
    

def PopInitialize(NumberDrops, Radius, Solutes, RNG,
                  Temperature, Pressure, AmbSat, Updraft,
                  BoxZLimits, ac, at, OAParameters, SurfMode, sft):
    """Function to create a droplet population array with position, velocity,
    and mass using initial conditions from solute and environmental data.
    
    NumberDrops  : unitless, integer.
    Radius       : m (Wet radii)
    Solutes      : Solute Array
    RNG          : Numpy random generator, object
    Temperature  : K (Ambient temperature)
    Pressure     : Pa (Ambient pressure)
    AmbSat       : % (Saturation ratio)
    Updraft      : m/s (Updraft velocity)
    BoxZLimits   : m (Vertical size of the simulation domain)
    ac, at       : Unitless (Kinetic effects coefficients)
    OAParameters : Organic aerosol parameters
    SurfMode     : String (Surface tension method)
    sft          : if SurfMode is 'Constant', value for it in J/m^2
    Output       : Droplet data array."""
    
    # 0:        Tracking number
    # 1:        Vertical Position
    # 2:        Vertical Velocity
    # 3:        Mass
    # 4:        Radius
    # 5:        Surface Saturation
    # 6:        Growth Rate
    # 7:        Surface Tension
    
    # Rename variables for ease
    N = NumberDrops
    T = Temperature
    P = Pressure
    S = AmbSat

    # Droplets array
    DropPop = np.zeros((N, 8), dtype=np.float64)

    # Tracking numbers
    DropPop[:,0] = np.arange(0, N, 1, dtype=np.int32)

    # Uniformly random positions with an updraft
    DropPop[:,1] = RNG.uniform(BoxZLimits[0], BoxZLimits[1], N)
    DropPop[:,2] = Updraft

    # Droplet wet radii
    DropPop[:,4] = Radius

    # Masses following m = rho V, take into account dry mass at center
    DropPop[:,3] = (4*np.pi/3)*rhoL*(Radius**3 - Solutes[:,1]**3)
    DropPop[:,3] = np.where(DropPop[:,3]>0, DropPop[:,3], 0)

    # Droplet saturations and surface tensions
    Sdrop, SurfTens = KappaKoehler(DropPop, Solutes, OAParameters, 
                                   T, SurfMode, sft)
    DropPop[:,5], DropPop[:,7] = Sdrop, SurfTens
    
    # Growth rates
    GRates = GrowthRate(DropPop, S, T, P, ac, at)
    DropPop[:,6] = GRates
    
    return DropPop


# =============================================================================
# Co-condensation of Organics
# =============================================================================

@nb.njit(cache=True)
def OASatPressure2(Temperature, MolarMass, C0):
    """Calculate the saturation pressure of an OA component using a known
    C0 value from Cappa and Jimenez (2010). This assumes that C0 is a fixed 
    constant.
    
    Temperature : K (Ambient temperature)
    MolarMass   : kg/mol (Molar mass)
    C0          : ug/m3 (Saturation concentraion)
    Output      : Pa (Saturation vapor pressure)"""

    # Saturation vapor pressure
    pSat = (R*Temperature*C0)/(MolarMass*1e9)
    
    return pSat


def DiffusivityOA(MolarMass):
    """Function to calculate the diffusivity of an organic aerosol using its
    molar mass, taken from Lim, Ho-Jin, Annmarie G. Carlton, and Barbara J. 
    Turpin (2005).
    
    MolarMass : kg/mol (Molar mass)
    Output    : m^2/s (diffusivity)."""
    
    # Diffusivity, converted into m^2/s
    Dv = (1.9*((MolarMass*1e3)**(-2/3)))*1e-4

    return Dv


@nb.njit(cache=True)
def VaporMoleFraction(OAVapor, OAParameters, Temperature, 
                      Pressure, AmbSat, Cubes):
    """Function to calculate the mole fraction of the organic aerosol in the
    domain.
    
    OAVapor             : Pa (SVO vapor pressures)
    OAParameters        : SVO parameter array
    Temperature         : K (Ambient temperature)
    Pressure            : Pa (Ambient pressure)
    AmbSat              : % (Saturation ratio)
    Cubes               : Unitless (Number of cubic centimeters)
    Output              : unitless (Mole fraction)."""
    
    # Rename variables for ease
    T = Temperature
    S = AmbSat
    P = Pressure
    
    p = np.zeros_like(OAVapor)
    nOA = np.zeros_like(OAVapor)
    
    # Convert vapor concentration (ug/m3) to partial pressure (Pa)
    for c in range(OAVapor.shape[0]):
        p[c,:] = OAVapor[c,:]*(R*T)/(OAParameters[c,0]*1e9)
        nOA[c,:] = (p[c,:]*Cubes*(1e-2)**3)/(R*T)
    
    # Saturation vapour pressure, Pa
    es = SVP(T, 'liquid')
    # Vapour pressure, Pa
    e = es*S
    # Dry air pressure, Pa
    Pdry = P - e - np.sum(p, axis=0)
    
    # Molar amounts
    nAir = Pdry*((Cubes*1e-2)**3)/(R*T)
    nWat = e*((Cubes*1e-2)**3)/(R*T)
    
    # Mole fraction
    x = np.zeros_like(OAVapor)
    for c in range(OAVapor.shape[0]):
        x[c,:] = nOA[c,:]/(nAir+nWat+np.sum(nOA, axis=0))
    
    return x


@nb.njit(cache=True)
def DropMoleFrac(MassWet, MassDry, MassOA, MolarMassDry, MolarMassOA):
    """Function to calculate the mole fraction of OA in a droplet.
    
    MassWet         : kg (Condensed water mass)
    MassDry         : kg (Inorganic mass)
    MassOA          : kg (Organic mass)
    MolarMassDry    : kg/mol (Inorganic molar mass)
    MolarMassOA     : kg/mol (Organic molar mass)
    Output          : unitless (Mole fraction)"""
    
    # Get moles
    nW = MassWet/Mw
    nD = MassDry/MolarMassDry
    
    # OA moles array
    nOs = np.zeros((MassWet.shape[0], MolarMassOA.shape[0]))
    for i in range(MolarMassOA.shape[0]):
        nOs[:,i] = MassOA[:,i]/MolarMassOA[i]
    
    # Total moles
    nTot = nW + nD + nOs.sum(axis=1)
    
    # Mole fraction array
    xs = np.zeros_like(nOs)
    # OA mole fraction
    for i in range(MolarMassOA.shape[0]):
        xs[:,i] = nOs[:,i]/nTot
    
    return xs


@nb.njit(cache=True)
def SatCondTempFunc(Temperature, Concentration, Tref):
    """Function to calculate the new saturation concentration as a function
    of temperature, based off of a Clausius-Clapeyron relation from Cappa
    and Jimenez (2010).
    
    Temperature     : K (Ambient temperature)
    Concentration   : ug/m3 (Organic saturation concentration)
    Tref            : K (Reference temperature)
    Output          : ug/m3 (New saturation concentration)."""
    
    # Delh_vap, enthalpy of vaporization term parameterization, in J/mol
    # from the same paper/Epstein et al (2009)
    DelH = (131 - 11*np.log10(Concentration))*1e3
    
    # New saturation concentration, ug/m3
    NewC = Concentration*(Tref/Temperature)*\
            np.exp(-(DelH/R)*(1/Temperature - 1/Tref))
    
    return NewC


@nb.njit(cache=True)
def Cocondense(Temperature, Tref, OAParameters, DropArray, 
                SoluteArray, SoluteDry, dt, Cubes):
    
    """Function to calculate the condensation equation for each organic species
    in the simulation.
    
    Temperature     : K (Ambient temperature)
    Tref            : K (Reference temperature, initial tmperature of 
                      simulation)
    OAParameters    : Organics parameter array
    DropArray       : Water data array
    SoluteArray     : Solute data array
    SoluteDry       : Dry solute array
    dt              : s (Time-step size)
    Cubes           : Unitless (Number of cubic centimeters)
    Output          : (1) Concentrations array (ug/m3), 
                      (2) Growth Rates array (kg/s),
                      (3) Organics data array."""
    
    # Rename variable for ease
    T = Temperature

    # Concentration, from ug/m3 to molecules/cm^3
    C = ((OAParameters[:,4]*1.0e-9)/OAParameters[:,0])*1.0e-6*NA
    # Saturation concentration, ug/m3
    CSat = SatCondTempFunc(T, OAParameters[:,3], Tref)
    # Convert C* to a pure component saturation vapour pressure [atm]
    P_sat = (CSat*Rd_alt*T)/(1.0e9*OAParameters[:,0])
    
    # Mole fractions, condensed
    x = DropMoleFrac(DropArray[:,3], SoluteDry[:,2], SoluteArray[:,6:], 
                      SoluteArray[:,3], OAParameters[:,0])
    
    # Kelvin effect factor
    KelFactor = np.zeros_like(x)
    # Equilibrium pressure at particle surface
    EquiPress = np.zeros_like(x)
    # Mean molar volume of each particle
    MolVol = np.zeros((x.shape[0], x.shape[1]+1))
    # Loop over organics
    for i in range(C.size):
        MolVol[:,i] = x[:,i]*OAParameters[i,0]
    # Include inorganic
    MolVol[:,-1] = (1-x.sum(axis=1))*SoluteDry[:,3]
    MolVol = MolVol.sum(axis=1)/SoluteArray[:,4]
    
    # Populate
    for i in range(C.size):
        # Kelvin effect
        KelFactor[:,i] = np.exp((2*DropArray[:,7]/(R*T*DropArray[:,4]))*MolVol)
        # Equilibrium pressure
        EquiPress[:,i] = x[:,i]*P_sat[i]*KelFactor[:,i]
    
    # Equilibrium concentration at particle surface in molecules/cm^3
    Csurf = EquiPress*(NA/(Rd_alt*1.0e6*T))
    
    # Condensation rates array, kg/s
    dmdt = np.zeros((DropArray[:,4].size, OAParameters.shape[0]))
    # Arrays to save means (necessary due to numba package)
    CSurfmeans = np.zeros_like(C)
    xmeans = np.zeros_like(C)
    
    for i in range(C.size):
        
        # Adjusted for mass accommodation effects, alpha=0.1 (Sahle et al 2013)           
        Dact = OAParameters[i,2]/(DropArray[:,4]/(DropArray[:,4] + \
                    (1.12e-7)/2) + (OAParameters[i,2]/(1.0*DropArray[:,4]))*\
                        np.sqrt(2.0*np.pi*OAParameters[i,0]/(R*T)))
        
        CSurfmeans[i] = np.mean(Csurf[:,i])
        xmeans[i] = np.mean(x[:,i])

        # Growth rate per droplet, in kg/s
        dmdt[:,i] = (4*np.pi*DropArray[:,4]*1e2*Dact*1e4*(C[i] - Csurf[:,i]))*\
                    OAParameters[i,0]/NA
    
    # Total growth rate for OA in kg/cm^3/s
    NewConc = np.where(OAParameters[:,4] - dt*(dmdt.sum(axis=0))*1.0e6*1.0e9 \
                       < OAParameters[:,7], 
                       OAParameters[:,4] - dt*(dmdt.sum(axis=0))*1.0e6*1.0e9, 
                       OAParameters[:,7])
    
    DatChem = np.zeros((4, C.size))
    DatChem[0,:] = C
    DatChem[1,:] = CSat
    DatChem[2,:] = CSurfmeans
    DatChem[3,:] = xmeans
    
    return NewConc, dmdt, DatChem


@nb.njit(cache=True)
def KappaChanger(MassSol, MassOa, RhoSol, RhoOA, KappaSol, KappaOA):
    """Function to calculate the new kappa of a droplet undergoing co-
    condensation, based off of Petters and Kreidenweis (2007) mixing-rule.
    
    MassSol     : kg (Inorganic mass)
    MassOA      : kg (Organic mass)
    RhoSol      : kg/m3 (Inorganic density)
    RhoOA       : kg/m3 (Organic density)
    KappaSol    : Unitless (Inorganic hygroscopicity)
    KappaOA     : Unitless (Organic hygroscopicity)
    Output      : Unitless (New hygroscopicity)"""
    
    # Volumes from mass and density, m^3
    Vd = MassSol/RhoSol
    Vos = MassOa/RhoOA
    
    # Total volume, m^3
    Vtot = Vd + Vos.sum(axis=1)
    # Volume fractions of solute and organics
    VfracD = Vd/Vtot
    Vfracos = np.zeros_like(MassOa)
    for i in range(MassOa.shape[1]):
        Vfracos[:,i] = Vos[:,i]/Vtot
    
    Vfracos = np.ascontiguousarray(Vfracos)
    KappaOA = np.ascontiguousarray(KappaOA)
    
    # New kappa mixing rule
    Mixture = np.dot(Vfracos, KappaOA)
    NewKappa = KappaSol*VfracD + Mixture
    
    return NewKappa


@nb.njit(cache=True)
def SoluteGrowth(Solutes, SolutesDry, OAGrowthRate, OAParameters, dt):
    """Function to calculate the change in dry radius and dry mass of each
    particle due to co-condensation of organic aerosols.
    
    Solutes      : Data array
    SolutesDry   : Data array, Inorganics
    OAGrowthRate : kg/s (Organic condensation growth rate)
    OAParameters : Data array
    dt           : s (Time-step size)
    Output       : (1) kg (mass), (2) kg/m3 (density), (3) m (radius), 
                    (4) kg (Organic mass) (5) unitless (hygroscopicity)."""
    
    # Initialize arrays
    CondMasses = np.zeros_like(OAGrowthRate)
    DensityFraction = np.zeros_like(OAGrowthRate)
    
    # Condensed mass of each OA so far, kg
    for i in range(OAGrowthRate.shape[1]):
        CondMasses[:,i] = Solutes[:,6+i] + OAGrowthRate[:,i]*dt
        # Prevent negative mass
        CondMasses[:,i] = np.where(CondMasses[:,i]>0, CondMasses[:,i], 
                                   Solutes[:,6+i])
        
    # New total solute mass, kg
    # Prevent total mass smaller than the original pure solute
    Md = np.where(SolutesDry[:,2] + CondMasses.sum(axis=1) < SolutesDry[:,2], 
                  SolutesDry[:,2], SolutesDry[:,2] + CondMasses.sum(axis=1))
    
    # Density fractions of each OA
    for i in range(OAGrowthRate.shape[1]):
        DensityFraction[:,i] = (CondMasses[:,i]/Md)*OAParameters[i,1]
    
    # New density, kg/m3
    Rhodry = (SolutesDry[:,2]/Md)*SolutesDry[:,4] + DensityFraction.sum(axis=1)
    
    # New radius, m
    rd = ((3/(4*np.pi))*(Md/Rhodry))**(1/3)
        
    # New kappa
    NewKappa = KappaChanger(SolutesDry[:,2], CondMasses, SolutesDry[:,4],
                            OAParameters[:,1], SolutesDry[:,5],
                            OAParameters[:,5])

    return Md, Rhodry, rd, CondMasses, NewKappa


def OrganicsArray(MolarMass, Concentration, Kappa, SurfTension, Densities, 
                  Cstar):
    """Function to initialize a parameter array for the OA components,
    using Topping et al. (2013) as a reference.
    
    MolarMass       : kg/mol (Molar mass)
    Concentration   : ug/m^3 (Total concentration)
    Kappa           : unitless (Hygroscopicity)
    SurfTension     : J/m^2 (Pure component surface tension)
    Densities       : kg/m3 (Density)
    CStar           : ug/m3 (Saturation concentration)
    Output          : Data array."""
    
    # Create the parameter file
    OAParameters = np.zeros((MolarMass.shape[0], 8))
    OAParameters[:,0] = MolarMass
    OAParameters[:,1] = Densities
    OAParameters[:,2] = DiffusivityOA(OAParameters[:,0])
    OAParameters[:,3] = Cstar
    OAParameters[:,4] = Concentration # Vapor
    OAParameters[:,5] = Kappa
    OAParameters[:,6] = SurfTension
    OAParameters[:,7] = Concentration # Total

    return OAParameters


def KappaSigmoidFit(OCRatio):
    """Function to estimate/approximate the kappa value of an organic using
    its O:C ratio. Sigmoidal fit provided by Prof. Andreas Zuend, 
    code provided by Camilo Serrano Damha."""
    
    # Parameters
    a = 0.0
    b = 27.0
    c = 0.33
    d = 1.0
    g = 0.026
    
    # Curve
    sigmoidal_kappa = d + (a-d)/(1.0+np.exp(OCRatio-c)**b)**g

    return sigmoidal_kappa


def MolecularCorridor(Vols, RNG):
    """Function to create approximate values for a given input of volatilities
    (C*), based off of the molecular corridor approach detailed in Sharaiwa et
    al. (2014).
    
    Vols : ug/m3 (Volatilities/log10 of Saturation Concentration)
    RNG : Numpy random generator, object
    Output: 
    (1) kg/mol (Molar mass) 
    (2) kg/m3 (Density)
    (3) Unitless (hygroscopicity)"""
    
    # Molar mass maxima
    Mmax = -30*Vols + 360
    # Molar mass minima
    Mmin = -10*Vols + 140
    # Create randomly selected molar masses
    MoMass = RNG.uniform(Mmin, Mmax, (1, Vols.size)).reshape(Vols.size)
    # Attain assumed O:C ratio
    OCRatio = 1 - (MoMass - Mmin)/(Mmax-Mmin)
    # Get H:C ratio using approximation from Heald et al. (2010)
    HCRatio = 2 - OCRatio
    # Get density using approximation from Kuwata et al. (2012)
    rho = 1000*(12 + HCRatio + 16*OCRatio)/(7 + 5*HCRatio + 4.15*OCRatio)
    # Approximate kappa using a sigmoidal kappa function
    k = KappaSigmoidFit(OCRatio)
    
    return MoMass, rho, k


# =============================================================================
# Data Sorting
# =============================================================================

def DataSaver(DataList, Name, CoCond):
    """Function to save output data stored in DataList as npy files:
    0 : Droplet data
    1 : Solute data
    2 : Critical radius/saturation
    3 : Organic aerosol parameters
    4 : Organic aerosol data
    5 : Saturation
    6 : Temperature
    7 : Pressure
    8 : Mean height
    9 : Time."""
    
    # Files = ['DropOutput', 'Solute', 'CritParams', 'OAParams', 
    #          'ChemData', 'SatTime', 'TempTime', 'PressTime', 
    #          'HeightTime', 'Time']
    
    if CoCond == True:
        np.savez(Name+'.npz', DropOutput=DataList[0], Solute=DataList[1],
                 SolFilter=DataList[2], CritParams=DataList[3], 
                 OAParams=DataList[4], ChemData=DataList[5], 
                 SatTime=DataList[6], TempTime=DataList[7], 
                 PressTime=DataList[8], HeightTime=DataList[9], 
                 Time=DataList[10], MontVars=DataList[11])
        
    else:
        np.savez(Name+'.npz', DropOutput=DataList[0], Solute=DataList[1],
                 CritParams=DataList[2], SatTime=DataList[3], 
                 TempTime=DataList[4], PressTime=DataList[5], 
                 HeightTime=DataList[6], Time=DataList[7],
                 MontVars=DataList[8])
    return


# =============================================================================
# Simulation
# =============================================================================

@nb.njit(cache=True)
def Simulator(DropPop, DropTrack, Solutes, dt, EnvEvolve, T, P, DomainZLims, 
             Updraft, Sdt, Pdt, Tdt, Zlims, Cubes, ac, at, S, SoluteFilter,
             OAParameters, SolTrack, ChemData, CoCond, RunTime, Instances, 
             Time, BreakTime, mtol, DTs, ErrTrack, SurfMode, sft):
    
    # Save initial time as reference temperature for organics
    Tref = T
    
    # Initialize time tracking
    time = 0
    timeNew = 0
            
    while time <= RunTime:
        
        # Update time
        time = timeNew
        timeNew += dt
        
        # If there is a NaN value, break the simulation
        if np.sum(np.isnan(DropPop)) != 0:
            # Save the timestep at which it occurs
            BreakTime[0] = int((np.abs(time-Instances)).argmin())
            
            break
        
        # Get terminal velocities
        TerminalVels = TerminalVel(DropPop[:,4], T, P, S)
        # Update velocities and positions
        DropPop[:,2] = Updraft - TerminalVels
        # Update vertical boundaries
        DomainZLims = np.median(DropPop[:,2])*dt + DomainZLims
        # Apply boundary conditions
        DropPop = BoundaryConds3D(DropPop, DomainZLims)
        
        Good = False
        
        while Good == False:
        
            if CoCond is True:
                
                # Co-condensation of organics
                TempOA, OAGrowth, Chems = \
                        Cocondense(T, Tref, OAParameters, DropPop, 
                                    Solutes, SoluteFilter, dt, Cubes)
                # Changes to dry mass and size
                TempSol2, TempSol4, TempSol1, TempSol6, TempSol5 = \
                    SoluteGrowth(Solutes, SoluteFilter, OAGrowth, 
                                  OAParameters, dt)
            
            # New droplet mass
            TempMass, Ndt, StepError = NewMass(DropPop, Solutes, S, T, P, 
                                                OAParameters, ac, at, dt, mtol,
                                                SurfMode, sft)
            
            TempDrop = DropPop.copy()
            TempDrop[:,3] = TempMass
            
            # Update radii                 
            TempRads = ((3/(4*np.pi)*(TempMass/rhoL)) + Solutes[:,1]**3)**(1/3) 
            TempRads = np.where(TempRads<Solutes[:,1], Solutes[:,1], TempRads)
            TempDrop[:,4]  = TempRads
            
            # Droplet saturation
            TempSdrop, TempSigma = KappaKoehler(TempDrop, Solutes, 
                                                OAParameters, T, SurfMode, sft)
            TempDrop[:,5], TempDrop[:,7]  = TempSdrop, TempSigma
            
            # New diffusional growth rate
            TempGRates = GrowthRate(TempDrop, S, T, P, ac, at)
            TempDrop[:,6] = TempGRates
        
            if np.sum(np.isnan(TempDrop)) == 0:
                if CoCond is True:
                    OAParameters[:,4] = TempOA
                    Solutes[:,2], Solutes[:,4] = TempSol2, TempSol4
                    Solutes[:,1] =  TempSol1
                    Solutes[:,6:], Solutes[:,5] = TempSol6, TempSol5
                    
                DropPop = TempDrop
                    
                Good = True
                
            else:
                dt = 0.0001*dt
                
        if EnvEvolve is True:
            # Environment
            # Update pressure level
            P = PressEvolution(T, P, S, Updraft, dt)
            # Update temeprature
            T = TempEvolution(T, P, S, Updraft, DropPop[:,6], dt, Cubes)
            # Update saturation ratio grid field
            S = EnvSatTime(T, P, S, Updraft, DropPop[:,6], dt, Cubes)
        
        # Update time-step size
        dt = Ndt
        
        # Data Tracking
        for i in range(Instances.size):
            if timeNew > Instances[i] and time < Instances[i]:
                
                print(timeNew)
                
                Time[i] = timeNew
                
                DTs[i]    = dt
                ErrTrack[i] = StepError
                
                Sdt[i]    = S
                Tdt[i]    = T
                Pdt[i]    = P
                Zlims[i,:] = DomainZLims
                            
                if CoCond is True:
                    ChemData[:,:,i] = Chems
                
                # Alloting data according to tracking number
                for Drop in range(DropTrack.shape[0]):
                    if DropTrack[Drop,0,i] in DropPop[:,0]:
                        DropTrack[Drop,1:,i] = \
                                        DropPop[:,1:][DropPop[:,0] == Drop]
                        
                        SolTrack[Drop,1:,i] = \
                                        Solutes[:,1:][Solutes[:,0] == Drop]
    
                                
    return   
