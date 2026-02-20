# The PARCELY cloud parcel model
PARCELY is a particle-resolved cloud parcel model that includes co-condensation of organic compounds and the formation of organic films around droplets.

# Installation
## Linux
To install PARCELY, either download or clone the contents of "PARECLY_V2_Fortran" into your desired directory. In a terminal navigate to the "PARCELY_V2" subfolder and run the command "make".
Afterwards, run the command "make clean", which will remove the temporary files (.mod and .o) created during the compiling. You should be left with an executable "parcely".
The makefile for PARCELY is quite simple and should be easy to change if there are issues with the compilation (for instance, changing to a different Fortran compiler).

## Windows
No installation guide has been made yet for Windows, however in the meantime one can download the source files and create a Visual Studio project (just as an example).

# Running PARCELY
To run PARCELY, one only needs to launch the executable. All options and flags are read from the editable .txt files. The output files are stored in the subfolder "PARCELY_Output", also as .txt files.

## Environment Input File
The file "parcely_inorganic_input_file.txt" contains all the input variables and options that specify information on the parcel evironment, run-time options, aerosol size distribution, and organic options.

1. **DistSeed** - A seed for the pseudo-random generator from which the aerosol size distribution is sampled. It is not guaranteed that the same seed will work across different machines owing to differences in
compilers or other software, but they will recreate the same distribution on the same computer. If the seed is set to 0, the seed itself will be randomly generated (useful for performing many runs with different seeds).
The seed used will be printed in the "inorganics_output.txt" file in the header column.

2. **Cubes** - The parcel volume in cubic centimeters. The default "1" therefore is 1 cm<sup>3</sup>. If one increases the parcel volume, the number of particles will scale accordingly so that the number concentration remains the same.

3. **Prs, Tmp, RH** - Initial values for pressure, temperature, and relative humidity. These three change with time and are outputted in the "environment_output.txt" file.

4. **W, ac, at** - Constant values for the updraft velocity, mass accommodation coefficient, and thermal accommodation coefficient. The latter pair are used in the KineticEffects subroutine (droplet_functions.f90), and the former in the equations
for temperature and relative humidity (env_equations.f90).

5. **RunTime, dt** - The total model runtime in seconds and the timestep size used for each call of the ODE.

6. **ntsteps** - The number of equally spaced timesteps to save and output to the files. This should be equal or larger than dt, the timestep size.

7. **DistType** - The distribution type used for all the aerosol size distribution. (1) is monotonic distributions (all particles have the same dry radius), (2) is normal distributions, and (3) is log-normal distributions.
Currently PARCELY does not have the option to mix multiple kinds of distributions, so all distributions must be of the same kind (e.g., you cannot have one population of monotonic and one of log-normal).

8. **npops** -  The number of individual size distributions. For example, if one is simulating a total particle population that is defined by a mix of three distinct distributions, npops is 3. **This must match the inorganic input file entries** (see next section).

9. **DistConcs, DistRads, DistStds** - The number concentration, mean radius (in micrometers), and standard deviation of each size distribution, separated by commas.

10. **Xi** - The mixing state of the inorganics. A value of 0 indicates that every particle is purely one compound, and a value of 1 indicates that every particle is an equal mix of all inorganic compounds used.

11. **INCLUDE_ORGANICS** - The flag to set whether organic compounds are used or not. If .FALSE., PARCELY ignores the organic input file entirely and solves only for the water growth of droplets and environmental variables.

12. **INCLUDE_COCONDENSE** - The flag to set whether the mass growth equations for organics are included in the solver. This is useful if one wants to test purely the effects of organic film formation without co-condensation.

13. **Msft** - The mode of surface tension used in the simulation. (1) is a constant value (the value set in sftc, usually set to that of pure water), (2) is the temperature-dependent surface tension of pure water, (3) is a mole-fraction weighted average of water and organics, and (4) is the organic film method.
**Options (3) and (4) cannot work if INCLUDE_ORGANICS is set to .FALSE.!**

14. **sftc** - The constant surface tension value to use if option 1 is selected in Msft.

## Inorganic Input File
The file "parcely_inorganic_input_file.txt" is the input file listing all of the "inorganic" compounds. It should be noted that while this file and the inputs are called "inorganics", they do not have to be limited to inorganics. One can use a surrogate organic compound as the inorganic - the only difference though is in the compound properties used.

1. **Name** - The name of the inorganic. This is mostly irrelevant in the code itself and primarily for user input clarity. It must only be 6 characters long.

2. **PopFractions** - The fraction of each size distribution attributed to this inorganic, separated by spaces. **This must match npops in the environment file.** As an example, if two size distributions of 100 particles have been set (npops = 2) and you have two inorganic compounds where one has 75 particles in the first distribution and 25 particles in the second distribution, and the second is identical but reversed, the input would be:
INORG1 0.75 0.25 ...
INORG2 0.25 0.75 ...

3. **MolarMass, Density, Kappa** - The molar mass, liquid density, and kappa/hygroscopicity values associated with that inorganic.

## Organic Input File
The file "parcely_organic_input_file.txt" is the input file listing all of the "organic" compounds.

1. **Name** - The name associated with the organic. This will be used to name the output file that contains the condensed-phase mass per particle for that organic. It must also be limited to 6 characters.

2. **Total** - The total concentration of that organic across both phases. The sum of all the condensed mass of this organic and the gas phase concentration should be equal to this value over the runtime of the simulation.

3. **C0, MolarMass, Density, Kappa, Sigma** - The pure-component values for the saturation concentration/volatility, molar mass, liquid density, kappa/hygroscopicity, and surface tension of the organic.

4. **InitFrac** - The initial condensed **mass fraction** of the organic for each size distribution, separated by spaces. PARCELY will adjust the inorganic moles accordingly such that the radius of each particle is maintained after incorporating the organics. The sum of all the initial condensed mass fractions for all the organics must be less than 1 (to keep inorganic mass). If the total calculated condensed mass of the organic is larger than the user-inputted total concentration, the gas-phase concentration is set to essentially zero.
