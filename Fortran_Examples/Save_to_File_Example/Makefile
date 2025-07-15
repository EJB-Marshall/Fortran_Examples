FC = gfortran
FFLAGS = -O3 # This is an optimisation flag
SRC = real_type_mod.o io.o

io: $(SRC)
	$(FC) $(FFLAGS) -o io $(SRC) 
%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
clean:
	rm *.o *.mod