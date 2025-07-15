FC = gfortran
FFLAGS = -O3 # This is an optimisation flag
SRC = real_type_mod.o arrays.o main.o

main: $(SRC)
	$(FC) $(FFLAGS) -o main $(SRC) 
%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
clean:
	rm *.o *.mod