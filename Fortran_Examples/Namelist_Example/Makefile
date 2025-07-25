FC = gfortran
FFLAGS = -Wall -Wextra -O3 # This is an optimisation flag
SRC = real_type_mod.o read_parameters_mod.o namelist_read.o

.PHONY: clean

main: $(SRC)
	$(FC) $(FFLAGS) -o $@ $(SRC) 
%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
clean:
	@rm *.o *.mod main