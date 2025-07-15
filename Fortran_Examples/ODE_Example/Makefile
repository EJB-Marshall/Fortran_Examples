FC = gfortran
FFLAGS = -O3 # This is an optimisation flag
SRC = real_type_mod.o rk2.o ode_solve.o

ode_solve: $(SRC)
	$(FC) $(FFLAGS) -o $@ $(SRC) 
%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
clean:
	rm *.o *.mod