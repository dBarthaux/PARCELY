# The PARCELY cloud parcel model
PARCELY is a particle-resolved cloud parcel model that includes co-condensation of organic compounds and the formation of organic films around droplets.

# Installation
## Linux
To install PARCELY, either download or clone the contents of "PARECLY_V2_Fortran" into your desired directory. In a terminal navigate to the "PARCELY_V2" subfolder and run the command "make".
Afterwards, run the command "make clean", which will remove the temporary files (.mod and .o) created during the compiling. You should be left with an executable "parcely".
The makefile for PARCELY is quite simple and should be easy to change if there are issues with the compilation (for instance, a different Fortran compiler).

## Windows
No installation guide has been made yet for Windows, however in the meantime one can download the source files and create a Visual Studio project (just as an example).

# Running PARCELY
To run PARCELY, one only needs to launch the executable. All options and flags are read from the editable .txt files. The output files are stored in the subfolder "PARCELY_Output", also as .txt files.

## Environment Input File
The file "parcely_inorganic_input_file.txt" contains all the input variables and options that specify information on the parcel evironment, run-time options, aerosol size distribution, and organic options.

1. DistSeed - A seed for the pseudo-random generator from which the aerosol size distribution is sampled. It is not guaranteed that the same seed will work across different machines owing to differences in
compilers or other software, but they will recreate the same distribution on the same computer. If the seed is set to 0, the seed itself will be randomly generated (useful for performing many runs with different seeds).
The seed used will be printed in the "inorganics_output.txt" file in the header column.

2. Cubes - The parcel volume in cubic centimeters. The default "1" therefore is 1 cm<sup>3</sup>. If one increases the parcel volume, the number of particles will scale accordingly so that the number concentration remains the same.

3. Prs, Tmp, RH - Initial values for pressure, temperature, and relative humidity. These three change with time and are outputted in the "environment_output.txt" file.

4. W, ac, at - Constant values for the updraft velocity, mass accommodation coefficient, and thermal accommodation coefficient. The latter pair are used in the KineticEffects subroutine (droplet_functions.f90), and the former in the equations
for temperature and relative humidity (env_equations.f90).

5. RunTime, dt - The total model runtime in seconds and the timestep size used for each call of the ODE.

6. ntsteps - The number of equally spaced timesteps to save and output to the files. This should be equal or larger than dt, the timestep size.

7. DistType - The distribution type used for all the aerosol size distribution. (1) is monotonic distributions (all particles have the same dry radius), (2) is normal distributions, and (3) is log-normal distributions.
Currently PARCELY does not have the option to mix multiple kinds of distributions, so all distributions must be of the same kind (e.g., you cannot have one population of monotonic and one of log-normal).

8. npops -  The number of individual size distributions. For example, if one is simulating a total particle population that is defined by a mix of three distinct distributions, npops is 3. **This must match the inorganic input file entries** (see next section).

9. DistConcs, DistRads, DistStds - The number concentration, mean radius (in micrometers), and standard deviation of each size distribution, separated by commas.

10. Xi - The mixing state of the inorganics. A value of 0 indicates that every particle is purely one compound, and a value of 1 indicates that every particle is an equal mix of all inorganic compounds used.

11. INCLUDE_ORGANICS - The flag to set whether organic compounds are used or not. If .FALSE., PARCELY ignores the organic input file entirely and solves only for the water growth of droplets and environmental variables.

12. INCLUDE_COCONDENSE - The flag to set whether the mass growth equations for organics are included in the solver. This is useful if one wants to test purely the effects of organic film formation without co-condensation.

13. Msft - The mode of surface tension used in the simulation. (1) is a constant value (the value set in sftc, usually set to that of pure water), (2) is the temperature-dependent surface tension of pure water, (3) is a mole-fraction weighted average of water and organics, and (4) is the organic film method.
**Options (3) and (4) cannot work if INCLUDE_ORGANICS is set to .FALSE.!**

14. sftc - The constant surface tension value to use if option 1 is selected in Msft.
