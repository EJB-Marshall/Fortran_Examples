FC = gfortran
FFLAGS = -O3 # This is an optimisation flag
SRC =  maths.o area.o 

maths: $(SRC)
	$(FC) $(FFLAGS) -o maths $(SRC) 
%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
maths.o: area.o
clean:
	rm *.o *.mod



