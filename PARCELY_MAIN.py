# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 11:40:20 2022

Input/Main file to run the PARCELY model.

@author: Dan Barthaux, McGill University
"""

# =============================================================================
# Imports
# =============================================================================

import numpy as np
import PARCELY_FUNCS as DF
import time

# =============================================================================
# Domain Set-Up
# =============================================================================

# Total cubic centimeters of the domain (e.g 5cm3, 4.3cm3 etc.)
Cubes = 1
# Height of boundary box side (m)
BoxHeight = (Cubes**(1/3))*(1e-2)
# Setting a random generator seed for reproducibility
SeedNo = 1996
RNG = np.random.default_rng(SeedNo)


# =============================================================================
# Physics Set-Up
# =============================================================================

# Allow environment to change or not
EnvEvolve = True
# Set initial pressure (Pa), from which domain height is set
P = 900.0e2
InitZ = 8.5e3*np.log(1013.25e2/P)
# Starting Temperature, K
T = 283
# Initial saturation, %
S = 0.96801
# Define boundaries
DomainZLims = np.array([0, BoxHeight], dtype='float64')+InitZ
# Updraft velocity, m/s
Updraft = 0.32
# Mass/thermal accommodation coefficients for kinetic effects
ac = 0.2
at = 0.2


# =============================================================================
# Time Set-Up
# =============================================================================

# Initial time-step, seconds
dt = 0.000001
# Total model runtime, seconds
RunTime = 0.00001
# Number of evenly spaced time instances to save
M = 10
# Times of instances
Instances = np.arange(0, RunTime+RunTime/M, RunTime/M)
# Actually measured time instances
Time = np.zeros_like(Instances)
# Tolerance for Runge-Kutta Cash-Karp adaptive time-step scheme
mtol = 1e-23


# =============================================================================
# Organics Set-Up
# =============================================================================

# Surface tension calculation method
Modes = ['OrgFilm', 'Constant', 'MixRule', 'WaterTemp']
# Surface tension of pure water J/m^2
sft = 0.072
# Selected method
SurfMode = Modes[0]

# Total concentration (all phases), ug/m3
Scaling = 1
Concentrations = np.array([0.1])

# Volatility/Saturation Concentration, ug/m3
LogVols = np.arange(-6, -5, 1, dtype=np.float64)
Cstar = 10.0**(LogVols)

# Organic parameters
# Surface tension, J/m2
SurfTension = 0.04*np.ones(Concentrations.size)

EstimateParameters = False

if EstimateParameters is False:
    # unitless
    Kappas = np.ones(Concentrations.size)*1e-5
    # kg/m3
    Densities = np.ones(Concentrations.size)*852
    # kg/mol
    MoMass = np.ones(Concentrations.size)*256.4e-3

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
Organics = True
# Allow co-condensation
CoCond = False
# Allow distribution in kappa values
KappaRand = False

# Distribution type
Distribution = 'lognormal'
# Number concentrations
Ns = np.array([226,134])*Cubes
# Mean radii, um
Rs = np.array([0.0196, 0.0695])
# Standard deviations, um
Stds = np.array([1.71, 1.7])
# Inorganics indices from CSV
Inorgs = np.array([[0],[2]])
# Population percentage of each inorganic
InorgPopPerc = np.array([[100],[100]])
# Organics in condensed phase from OAParameter table
Orgs = [[0],[0]]

# If irrelevant, set to 0
PercCond = 0
OrgBase = np.array([100, 100, 100, 99, 99, 90, 45, 12, 1, 0.1])

OrgMFr = RNG.dirichlet(OrgBase, 1).reshape(10)*(1-0.6)
# Solute mass fractions (inorganic + organic)
MassFractions = [np.array([0.8, 0.2]), np.array([0.8, 0.2])] 

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
# Koehler Curves and Critical/Equilibrium parameters
# =============================================================================

start = time.time()
# Equilibrium radius, critical parameters for activation
EquiR, CritParams = DF.CritParamFinder(Solutes, T, S)

end = time.time()
print('Finding Critical Parameters: ', np.round(end - start, 4), 's')


# =============================================================================
# Droplet Initialization
# =============================================================================

start = time.time()
# Droplet array
DropPop = DF.PopInitialize(NumberDrops, EquiR, Solutes, RNG,
              T, P, S, Updraft, DomainZLims, ac, at, 
              OAParameters, SurfMode, sft)

end = time.time()
print('Droplets Initialization: ', np.round(end - start, 4), 's')


Test = DF.WaterSurfaceTension(T)


# =============================================================================
# Tracking Data
# =============================================================================

# Environmental trackers
Sdt     = np.zeros(M+1)     # Saturation
Tdt     = np.zeros(M+1)     # Temperature
Pdt     = np.zeros(M+1)     # Pressure

Sdt[:]  = np.nan
Tdt[:]  = np.nan
Pdt[:]  = np.nan

# Initial conditions
Sdt[0] = S
Tdt[0] = T
Pdt[0] = P

# Droplet and Solute Tracker
DropTrack = np.empty(shape=(NumberDrops, 8, M+1))
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
# Track time-step size
DTs = np.zeros_like(Tdt)
# Track error size
ErrTrack = np.zeros_like(Tdt)
# Model run parameters
RunParams = np.zeros(4)
RunParams[0] = SeedNo
RunParams[1] = Updraft
RunParams[2] = ac
RunParams[3] = at


# =============================================================================
# Simulation
# =============================================================================

start = time.time()

DF.Simulator(DropPop, DropTrack, Solutes, dt, EnvEvolve, T, P, DomainZLims, 
              Updraft, Sdt, Pdt, Tdt, Zlims, Cubes, ac, at, S, SoluteFilter,
              OAParameters, SolTrack, ChemData, CoCond, RunTime, Instances, 
              Time, BreakTime, mtol, DTs, ErrTrack, SurfMode, sft)

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

MeanHeight = np.nanmean(DropTrack[:,1,:], axis=0)
        
if Time[-1] == 0:
    Time = Time[:-1]

if CoCond == True:

    AllData = np.array([DropTrack, SolTrack, SoluteFilter, CritParams, 
                        OAParameters, ChemData, Sdt, Tdt, Pdt, 
                        MeanHeight-InitZ, Time, RunParams], dtype=object)

else:
    AllData = np.array([DropTrack, 
                        SolTrack[:,:,0].reshape(NumberDrops,Solutes.shape[1],1), 
                        CritParams, Sdt, Tdt, Pdt, 
                        MeanHeight-InitZ, Time, RunParams], dtype=object)


# =============================================================================
# Data Saving
# =============================================================================

# DF.DataSaver(AllData, 'Save1', CoCond)
