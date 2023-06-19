# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 11:40:20 2022


Main file to run the Parcely model.

@author: Dan Barthaux
"""

# =============================================================================
# Imports
# =============================================================================

import numpy as np
import PARCELY_FUNCS_HET as DF
import time


# =============================================================================
# Domain Set-Up
# =============================================================================

# Total cubic centimeters of the domain (e.g 5cm3, 4.3cm3 etc.)
Cubes = 1

# Domain Parameters
BoxLength = (Cubes**(1/3))*(1e-2)       # Length of boundary box side (m)
BoxHeight = (Cubes**(1/3))*(1e-2)       # Height of boundary box side (m)

# Setting a random generator seed for reproducibility
SeedNo = 1996
RNG = np.random.default_rng(SeedNo)

# Size threshold for collision-coalescence
SizeThreshold = 15e-6
# Number of nearest neighbors to keep
k = 5
# Allow collision-coalescence
CollCoal = True


# =============================================================================
# Physics Set-Up
# =============================================================================

# Allow environment to change or not
EnvEvolve = True
# Set initial pressure (Pa), from which domain height is set
P = 900.0e2
InitZ = 8.5e3*np.log(1013.25e2/P)
# Starting Temperature, K
T = 283.0
# Initial saturation, %
S = 0.96801
# Define boundaries using array
DomainXLims = np.array([0, BoxLength])
DomainYLims = np.array([0, BoxLength])
DomainZLims = np.array([0, BoxHeight], dtype='float64')+InitZ
# Initial updraft velocity, and update initial velocities
WComp = 1.0
Updraft = np.array([0, 0, WComp])
# Mass/thermal accommodation coefficients for kinetic effects
ac = 0.2
at = 0.2


# =============================================================================
# Time Set-Up
# =============================================================================

# Initial timestep, seconds
dt = 0.000001
# Total model runtime wanted in seconds
RunTime = 0.00001
# Number of evenly spaced time instances to save
M = 10
# Times of instances
Instances = np.arange(0, RunTime+RunTime/M, RunTime/M)
# Actually measured time instances
Time = np.zeros_like(Instances)
# Tolerance
mtol = 1e-23


# =============================================================================
# Organic Aerosol Initialization
# =============================================================================

# Surface tension calculation method
Modes = ['OrgFilm', 'Constant', 'MixRule', 'WaterTemp']
# Surface tension of pure water J/m^2
sft = 0.072
# Selected method
SurfMode = Modes[1]

# Total concentration (all phases), ug/m3
Scaling = 1
Concentrations = np.array([0.1, 0.1, 0.1, 0.12, 0.12, 0.12, 0.15,
                            0.24, 0.48, 1.18])*Scaling

# Volatility/Saturation Concentration, ug/m3
LogVols = np.arange(-6, 4, 1, dtype=np.float64)
Cstar = 10.0**(LogVols)

# Organic parameters
# Surface tension, J/m2
SurfTension = 0.04*np.ones(Concentrations.size)

EstimateParameters = False

if EstimateParameters is False:
    # unitless
    Kappas = np.ones(Concentrations.size)*0.1
    # kg/m3
    Densities = np.ones(Concentrations.size)*1400
    # kg/mol
    MoMass =  np.linspace(400, 200, Concentrations.size)*1e-3

if EstimateParameters is True:
    MoMass, Densities, Kappas = DF.MolecularCorridor(LogVols, RNG)

# OA parameter table
OAParameters = DF.OrganicsArray(MoMass, Concentrations, Kappas,
                                SurfTension, Densities, Cstar)    


# =============================================================================
# Solutes Set-Up
# =============================================================================

# Dirichlet mass fractions for non-monotone distribution
Dirichlet = False
# Allow organics
Organics = False
# Allow co-condensation
CoCond = False
# Allow distribution in kappa values
KappaRand = False

# Distribution type
Distribution = 'mono'
# Number concentrations
Ns = np.array([100])*Cubes
# Mean radii, um
Rs = np.array([0.05])
# Standard deviations, um
Stds = np.array([1])
# Inorganics indices from CSV
Inorgs = np.array([[0]])
# Population percentage of each inorganic
InorgPopPerc = np.array([[100]])
# Organics in condensed phase from OAParameter table
Orgs = [[0]]

# If irrelevant, set to 0
PercCond = 0
OrgBase = np.array([])

# OrgMFr = RNG.dirichlet(OrgBase, 1).reshape(10)*(1-0.6)
# Solute mass fractions (inorganic + organic)
# MassFractions = [np.append(np.array([0.6]), OrgMFr), 
#                  np.append(np.array([0.6]), OrgMFr)] 
MassFractions = [np.array([1]), np.array([1])]

# Number of droplets
NumberDrops = Ns.sum()

start = time.time()

# Get solute arrays
Solutes, SoluteFilter, VaporConcs = DF.SolPopInitialize(Ns, Rs, Stds,
                                        Inorgs, InorgPopPerc, Orgs, PercCond, 
                                        OAParameters, RNG, MassFractions, 
                                        Dirichlet, KappaRand, Organics,
                                        Cubes, Distribution)

# Get vapor concentrations
OAParameters[:,4] = VaporConcs

end = time.time()
print('Solute/OA Initialization: ', np.round(end - start, 4), 's')


# =============================================================================
# Saturation/Temperature/Pressure Field
# =============================================================================

# By how much the air parcel's side is subdivided
SatDivide = 5
SatField, SatGrid, SatInd = DF.SatFieldInit(S, SatDivide, DomainXLims, 
                                            DomainYLims, DomainZLims,
                                            RNG, 'mono')
    
Diffusion = False
DropMove = True

TempField = np.ones_like(SatField)*T
PressField = np.ones_like(SatField)*P


# =============================================================================
# Koehler Curves and Critical/Equilibrium parameters
# =============================================================================

start = time.time()

EquiR, CritParams = DF.CritParamFinder(Solutes, T, S)

end = time.time()
print('Finding Critical Parameters: ', np.round(end - start, 4), 's')


# =============================================================================
# Droplet Initialization
# =============================================================================

start = time.time()
DropPop, OrigVertVel = DF.PopInitializeF(NumberDrops, EquiR, Solutes, RNG, 
                                          TempField, PressField, 
                                          Updraft, SatField, SatGrid, SatInd,
                                          BoxLength, DomainZLims, 
                                          ac, at, OAParameters, SurfMode, sft)

end = time.time()
print('Droplets Initialization: ', np.round(end - start, 4), 's')


# =============================================================================
# Tracking Data
# =============================================================================

# Environmental trackers
Sdt = np.zeros((M+1, SatDivide, SatDivide, SatDivide)) # Saturation ratio
Tdt = Sdt.copy()
Pdt = Sdt.copy()

Sdt[:]  = np.nan
Tdt[:]  = np.nan
Pdt[:]  = np.nan

# Initial conditions
Sdt[0,:,:,:]  = SatField
Tdt[0,:,:,:] = TempField
Pdt[0,:,:,:] = PressField

# Droplet and Solute Tracker
DropTrack = np.empty(shape=(NumberDrops, 12, M+1))
SolTrack = np.empty(shape=(NumberDrops, Solutes.shape[1], M+1))

DropTrack[:,:,:] = np.nan
SolTrack[:,:,:] = np.nan

# Initial step
DropTrack[:,:,0] = DropPop
SolTrack[:,:,0] = Solutes

# Keep tracking numbers throughout
for t in range(M):
    DropTrack[:,0,t+1] = DropTrack[:,0,0]
    SolTrack[:,[0,3,5],t+1]  = SolTrack[:,[0,3,5],0]

# Height limits counter
Zlims = np.zeros(shape=(M+1,2))
Zlims[0,:] = DomainZLims

ChemData = np.zeros((4, Concentrations.size, M+1))
ChemData[:,:,0] = DF.Cocondense(T, T, OAParameters, DropPop, 
                            Solutes, SoluteFilter, dt, Cubes)[2]

BreakTime = np.zeros(1).astype(np.int32)

DTs = np.zeros_like(Tdt)
ErrTrack = np.zeros_like(Tdt)
tolTrack = np.zeros_like(Tdt)
tolTrack[0] = mtol

MontVars = np.zeros(4)
MontVars[0] = SeedNo
MontVars[1] = WComp
MontVars[2] = ac
MontVars[3] = at


# =============================================================================
# Simulation
# =============================================================================

start = time.time()

DF.SimulatorF(DropPop, DropTrack, Solutes, OrigVertVel, dt, EnvEvolve,
              TempField, PressField, DomainZLims, DomainXLims, DomainYLims, 
              Updraft, Sdt, Pdt, Tdt, Zlims, Cubes, ac, at, S, SoluteFilter,
              OAParameters, SolTrack, ChemData, CoCond, RunTime,
              Instances, Time, BreakTime, mtol, DTs, ErrTrack,
              tolTrack, SurfMode, SatGrid, SatDivide, SatInd, SatField,
              sft, Diffusion, DropMove, SizeThreshold, k, CollCoal)


end = time.time()
print('Model Run: ', end - start, 's')

if BreakTime[0] != 0:
    Time   = Time[:BreakTime[0]]
    Sdt         = Sdt[:BreakTime[0]]
    Tdt         = Tdt[:BreakTime[0]]
    Pdt         = Pdt[:BreakTime[0]]
    Zlims       = Zlims[:BreakTime[0]]
    DropTrack   = DropTrack[:,:,:BreakTime[0]]
    SolTrack    = SolTrack[:,:,:BreakTime[0]]
    
    if CoCond is True:
        ChemData    = ChemData[:,:,:BreakTime[0]]


MeanHeight = np.nanmean(DropTrack[:,3,:], axis=0)
        
if Time[-1] == 0:
    Time = Time[:-1]

if CoCond == True:

    AllData = np.array([DropTrack, SolTrack, SoluteFilter, CritParams, 
                        OAParameters, ChemData, Sdt, Tdt, Pdt, 
                        MeanHeight-InitZ, Time, MontVars], dtype=object)

else:
    AllData = np.array([DropTrack, SolTrack[:,:,0].reshape(NumberDrops*Cubes,6,1), 
                        CritParams, Sdt, Tdt, Pdt, 
                        MeanHeight-InitZ, Time, MontVars], dtype=object)


# =============================================================================
# Data Saving
# =============================================================================

# DF.DataSaver(AllData, 'File', CoCond)