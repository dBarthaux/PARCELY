# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 10:43:35 2026

@author: Dan Barthaux
"""

import os
import warnings
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from sklearn.cluster import KMeans

# =============================================================================
# 
# =============================================================================

def ReadPARCELY(Path):

    # Read in environment file
    EnvFile = pd.read_csv(Path+'environment_output.txt',  sep='   ',  engine='python')
    # Inorganic data
    Inorganics = pd.read_csv(Path+'inorganics_output.txt',  sep='   ',  engine='python')
    
    # Droplets (water radius)
    DropRadius = open(Path+'droplet_output.txt').read().splitlines()[1:]
    for i, ival in enumerate(DropRadius):
        DropRadius[i] = ival.split()
    # Convert to micrometers
    DropRadius = np.array(DropRadius, dtype=float)
    DropRadius[:,1:] *= 1e6
    
    # Activation, just number of droplets activated
    ToActivate = DropRadius[:,1:].T
    ActRez = ActivateDetectVisNew(ToActivate, 30)
    
    # Solute properties
    SolPropNames = ['soldensity', 'solkappa', 'solmass',
                    'solmolarmass', 'solradius', 'soltension']
    
    SolData = []
    for name in SolPropNames:
        try:
            SolProperty = open(Path+f'{name}_output.txt').read().splitlines()[1:]
            for i, ival in enumerate(SolProperty):
                SolProperty[i] = ival.split()
            SolProperty = np.array(SolProperty, dtype=float)
            SolData.append(SolProperty)
        except:
            print(f'File {name} not found.')
            SolData.append([])
    
    try:
        SolData = np.array(SolData)
    except:
        pass
    
    Used = SolPropNames + ['environment', 'inorganics', 'droplet']
    
    # Any remaining files assumed to be organic
    Orgfiles = [f for f in os.listdir(Path) if f.split('_')[0] not in Used]
    
    OrgData = {}
    for org in Orgfiles:
        Organic = open(Path+org).read().splitlines()[1:]
        for i, ival in enumerate(Organic):
            Organic[i] = ival.split()
            
        try:
            OrgData[org] = np.array(Organic, dtype=float)
        except:
            for line in Organic:
                for v, val in enumerate(line):
                    if 'E' not in val:
                        line[v] = '0.000000E-00'
            OrgData[org] = np.array(Organic, dtype=float)
    
    Output = {'Env':EnvFile, 'Inorg':Inorganics, 'Drop':DropRadius,
              'Sol':SolData, 'Org':OrgData, 'Act':ActRez}

    return Output


# =============================================================================
# 
# =============================================================================

def ActivateDetectVisNew(Radii, Offset, Method='argmax'):
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
    # Radii = DropData[:,4,:]
    # Smooth radii
    RadFilt = savgol_filter(Radii, 51, 5)
    
    if Method == 'kmeans':
        
        # Get final time-step of every droplet's radius
        End = Radii[:,-1].reshape(-1,1)
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
            if Radii[Group1,-1].mean() > Radii[Group2,-1].mean():
                Activated = Group1
                UnActivated = Group2
            else:
                Activated = Group2
                UnActivated = Group1
    
    if Method == 'argmax':
        Activated = []
        UnActivated = []
        for d in range(Radii.shape[0]):
            if (np.argmax(RadFilt[d]) == RadFilt.shape[1]-1) and (RadFilt[d,-1] > 0.05):
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
