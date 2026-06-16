# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:21:19 2026

@author: Dan Barthaux
"""

import Python_Funcs_PARCELY_V2 as PV
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# =============================================================================
# Load the output data
# =============================================================================

Data = PV.ReadPARCELY('Test Data/')

N = Data['Inorg'].shape[0]
Time = Data['Drop'][:,0]

# =============================================================================
# Environmental Plots
# =============================================================================

Fig = plt.figure(figsize=(13,4))
ax = Fig.subplots(nrows=1, ncols=3)
ax = ax.flatten()

Fig.suptitle('Environmental Variables', fontsize=16)

ax[0].plot(Data['Env']['Time (s)'], 100*(Data['Env']['RH (-)']-1), linewidth=3)
ax[1].plot(Data['Env']['Time (s)'], Data['Env']['Temperature (K)'], linewidth=3)
ax[2].plot(Data['Env']['Time (s)'], Data['Env']['Pressure (Pa)']/100, linewidth=3)

ax[0].set_ylabel('Supersaturation [-]', fontsize=14)
ax[1].set_ylabel('Temperature [K]', fontsize=14)
ax[2].set_ylabel('Pressure [hPa]', fontsize=14)

Fig.supxlabel('Time [s]', fontsize=15)

plt.tight_layout()


# =============================================================================
# Inorganic Plots
# =============================================================================

Fig = plt.figure(figsize=(12,6))
ax = Fig.subplots(nrows=2, ncols=2)
ax = ax.flatten()

Fig.suptitle('Inorganic/Dry Particle Properties', fontsize=16)

for c, col in enumerate(Data['Inorg'].columns[1:5]):
    if col == 'Radius (m)':
        ax[c].hist(Data['Inorg'][col]*1e9, bins=50)
        
        if np.all(Data['Inorg'][col] == Data['Inorg'][col].iloc[0]):
            ax[c].set_xlim([Data['Inorg'][col].iloc[0]*1e9-1,
                            Data['Inorg'][col].iloc[0]*1e9+1])
        
    elif col == 'Molar Mass (kg/mol)':
        ax[c].hist(Data['Inorg'][col]*1e3, bins=50)
        
        if np.all(Data['Inorg'][col] == Data['Inorg'][col].iloc[0]):
            ax[c].set_xlim([Data['Inorg'][col].iloc[0]*1e3-1,
                            Data['Inorg'][col].iloc[0]*1e3+1])
        
    else:
        ax[c].hist(Data['Inorg'][col], bins=50)
        
        if np.all(Data['Inorg'][col] == Data['Inorg'][col].iloc[0]):
            ax[c].set_xlim([Data['Inorg'][col].iloc[0]-1,
                            Data['Inorg'][col].iloc[0]+1])

ax[0].set_xlabel('Radius [nm]', fontsize=14)
ax[1].set_xlabel('Molar Mass [$\\mathrm{g~mol^{-1}}$]', fontsize=14)
ax[2].set_xlabel('Density [$\\mathrm{kg~m^{-3}}$]', fontsize=14)
ax[3].set_xlabel('Kappa [-]', fontsize=14)

Fig.supylabel('Particle Count', fontsize=15)

plt.tight_layout()


# =============================================================================
# Activation Plots
# =============================================================================

LegLines = [Line2D([0], [0], color='blue', lw=3),
            Line2D([0], [0], color='gray', lw=3)]

plt.figure()
plt.title('Droplet Activation', fontsize=16)
for i in range(1,N+1):
    if i-1 in Data['Act'][4]:
        plt.plot(Data['Drop'][:,0], Data['Drop'][:,i], c='blue', alpha=0.7)
    else:
        plt.plot(Data['Drop'][:,0], Data['Drop'][:,i], c='gray', alpha=0.7)
        
plt.legend(LegLines, ['Activated', 'Unactivated'], prop={'size':12})
plt.yscale('log')
plt.xlabel('Time [s]', fontsize=15)
plt.ylabel('Wet Radius [$\\mathrm{\mu~m}$]', fontsize=15)
plt.tight_layout()
    

Fig = plt.figure(figsize=(7,5))
ax = Fig.subplots()
Fig.suptitle('Droplet Activation', fontsize=16)
ax.plot(Data['Env']['Time (s)'], Data['Act'][1], linewidth=4)
ax2 = ax.twinx()
ax2.plot(Data['Env']['Time (s)'], 100*Data['Act'][1]/N, alpha=0.0)

ax.set_xlabel('Time [s]', fontsize=15)
ax.set_ylabel('Activated Droplet Count', fontsize=15, labelpad=15)
ax2.set_ylabel('Activated Droplet Percentage [%]', fontsize=15, labelpad=15)
plt.tight_layout()
    

# =============================================================================
# Solute Temporal Plots
# =============================================================================

Fig = plt.figure(figsize=(10,8))
ax = Fig.subplots(nrows=3, ncols=2)
ax = ax.flatten()

Fig.suptitle('Solute Properties', fontsize=16)

for c in range(6):
    
    if c in [2,4]:
        ax[c].plot(Data['Sol'][c,:,0], 1e9*Data['Sol'][c,:,1:].mean(axis=1),
                   c='k', linewidth=3, zorder=6)
        for i in range(1, N+1):
            ax[c].plot(Data['Sol'][c,:,0], 1e9*Data['Sol'][c,:,i], alpha=0.4, c='gray')
        
        ax[c].set_yscale('log')
    
    else:
        ax[c].plot(Data['Sol'][c,:,0], Data['Sol'][c,:,1:].mean(axis=1),
                   c='k', linewidth=3, zorder=6, label='Mean')
        for i in range(1, N+1):
            ax[c].plot(Data['Sol'][c,:,0], Data['Sol'][c,:,i], alpha=0.4, c='gray')

ax[0].legend()

ax[0].set_ylabel('Density [$\\mathrm{kg~m^{-3}}$]')
ax[1].set_ylabel('Kappa [-]')
ax[2].set_ylabel('Mass [$\\mathrm{\mu g}$]')
ax[3].set_ylabel('Molar Mass [$\\mathrm{g~mol^{-1}}$]')
ax[4].set_ylabel('Radius [nm]')
ax[5].set_ylabel('Surface Tension [$\\mathrm{J~m^{-2}}$]')

Fig.supxlabel('Time [s]', fontsize=16)
plt.tight_layout()


plt.figure()
plt.title('Solute Radius Change', fontsize=16)
plt.hist(Data['Sol'][4,0,1:]*1e9, bins=50, facecolor='none', 
         edgecolor='black', linewidth=1.5, label='Start')
plt.hist(Data['Sol'][4,-1,1:]*1e9, bins=50, facecolor='none', 
         edgecolor='red', linewidth=1.5, label='End')
plt.xlabel('Solute Radius [nm]', fontsize=15)
plt.ylabel('Particle Count', fontsize=15)
plt.legend(prop={'size':12})
plt.tight_layout()


# =============================================================================
# Organics Plots
# =============================================================================

Names = list(Data['Org'].keys())
EndName = Names[-1]
Names = [name.split('_')[0] for name in Names[:-1]]

plt.figure()
plt.title('Organic Gas Concentration', fontsize=16)
for i in range(Data['Org'][EndName].shape[1]-1):
    plt.plot(Data['Org'][EndName][:,0], Data['Org'][EndName][:,i+1],
             label=Names[i], linewidth=2)
plt.yscale('log')
plt.legend()
plt.xlabel('Time [s]', fontsize=15)
plt.ylabel('[$\\mathrm{molec~cm^{-3}}$]', fontsize=15)
plt.tight_layout()


# =============================================================================
# Mole Fraction Plots
# =============================================================================

Na = 6.0221408e23
Mw = 18.01528e-3
rhoL = 997.0

# Moles of each component (water, inorganic, organic)
InorgMoles = (Data['Inorg']['Mass (kg)']/Data['Inorg']['Molar Mass (kg/mol)']).values
WaterMass = (4/3)*np.pi*((Data['Drop'][:,1:]*1e-6)**3 - (Data['Sol'][4,:,1:])**3)*rhoL
WaterMoles = WaterMass/Mw
OrgMoles = np.array([Data['Org'][key][:,1:]/Na for key in list(Data['Org'].keys())[:-1]])

# Total moles over time for each droplet
TotalMoles = InorgMoles[np.newaxis,:] + WaterMoles + OrgMoles.sum(axis=0)

# Mole fractions
InorgMF = InorgMoles[np.newaxis,:]/TotalMoles
WaterMF = WaterMoles/TotalMoles
OrgnsMF = np.array([OrgMoles[i]/TotalMoles for i in range(OrgMoles.shape[0])])



Fig = plt.figure(figsize=(13,4.5))
ax = Fig.subplots(nrows=1, ncols=3)
ax = ax.flatten()

Fig.suptitle('Particle Component Mole Fractions', fontsize=16)

for i in range(TotalMoles.shape[1]):
    if i in Data['Act'][4]:
        ax[0].plot(OrgnsMF[:,:,i].sum(axis=0), WaterMF[:,i], c='blue', alpha=0.5)
        ax[1].plot(InorgMF[:,i], OrgnsMF[:,:,i].sum(axis=0), c='blue', alpha=0.5)
        ax[2].plot(WaterMF[:,i], InorgMF[:,i], c='blue', alpha=0.5)
    else:
        ax[0].plot(OrgnsMF[:,:,i].sum(axis=0), WaterMF[:,i], c='gray', alpha=0.5)
        ax[1].plot(InorgMF[:,i], OrgnsMF[:,:,i].sum(axis=0), c='gray', alpha=0.5)
        ax[2].plot(WaterMF[:,i], InorgMF[:,i], c='gray', alpha=0.5)

ax[0].legend(LegLines, ['Activated', 'Unactivated'], prop={'size':12})

ax[0].set_xlabel('Organic', fontsize=13)
ax[1].set_xlabel('Inorganic', fontsize=13)
ax[2].set_xlabel('Water', fontsize=13)

ax[0].set_ylabel('Water', fontsize=13)
ax[1].set_ylabel('Organic', fontsize=13)
ax[2].set_ylabel('Inorganic', fontsize=13)

plt.tight_layout()


plt.figure()
plt.title('Average Organic Mole Fractions', fontsize=15)
for i in range(OrgnsMF.shape[0]):
    plt.plot(Time, OrgnsMF[i].mean(axis=1), label=Names[i], linewidth=3)
plt.yscale('log')
plt.ylim([1e-6, 1e0])
plt.legend()
plt.xlabel('Time [s]', fontsize=15)
plt.ylabel('Mole Fraction', fontsize=15)
plt.tight_layout()


