FC = gfortran
FFLAGS = -O3 -Wall # This is an optimisation flag
SRC = real_type_mod.o read_parameters_mod.o \
	diff_op_mod.o boundary_conditions_mod.o \
	burgers_mod.o time_stepping_mod.o \
	evolve_burgers.o

.PHONY: clean

main: $(SRC)
	$(FC) $(FFLAGS) -o $@ $(SRC) 
%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
clean:
	@rm *.o *.mod main