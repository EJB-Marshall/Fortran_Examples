import subprocess
import numpy as np


### Define Grid Parameters

nml_1 = "&Grid_Params"
x_min = 0.0
x_max = 2.0*np.pi
Nx = 2000
Ngz = 2
Nvar = 1


nml_2 = "&Sim_Params"
t0 = 0.0
tend = 20.0


### Create namelist file
f = open("params.nml","w")

f.write(nml_1+"\n")
f.write("x_min = " + str(x_min) + "\n")
f.write("x_max = " + str(x_max) + "\n")
f.write("Nx = " + str(Nx) + "\n")
f.write("Ngz = " + str(Ngz) + "\n")
f.write("Nvar = " + str(Nvar) + "\n")
f.write("/\n")
f.write("\n")
f.write(nml_2+"\n")
f.write("t0 = " + str(t0) + "\n")
f.write("tend = " + str(tend) + "\n")
f.write("/\n")
f.close() # Need to close namelist file before running Fortran!


### Run our Fortran script
subprocess.run("./main")