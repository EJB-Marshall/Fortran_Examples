FC = gfortran
FFLAGS = -O3 # This is an optimisation flag
SRC = real_type_mod.o finite_difference_mod.o \
	boundary_conditions_mod.o advection_mod.o \
	time_stepping_mod.o evolve_advection.o

.PHONY: clean

main: $(SRC)
	$(FC) $(FFLAGS) -o $@ $(SRC) 
%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
clean:
	@rm *.o *.mod main