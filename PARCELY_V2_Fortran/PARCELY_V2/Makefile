FC=gfortran
FFLAGS=-O3 -ffree-line-length-none
SRC=dvode_kinds_module.f90 dvode_blas_module.f90 dvode_linpack_module.f90 dvode_module.f90 fmin_module.f90 root_module.f90 \
env_constants.f90 water_functions.f90 air_functions.f90 distribution_functions.f90 model_parameters.f90 env_equations.f90 \
input_read.f90 droplet_functions.f90 solute_functions.f90 dydt_main.f90 parcely_program.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

parcely: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod