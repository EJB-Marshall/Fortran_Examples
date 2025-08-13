import numpy as np
from matplotlib import pyplot as plt

fortran_soln  = np.fromfile("Burgers_soln.dat", dtype=np.float64)
fortran_soln  = fortran_soln.reshape(-1, 2000)
fortran_times = np.fromfile("Burgers_soln_times.dat", dtype=np.float64)
grid = np.linspace(0,2*np.pi,np.shape(fortran_soln)[1]-4) # Subtract the number of ghost points from array size


### Plotting parameters I use when saving a single figure 
### for a paper

# plt.rcParams.update({"text.usetex": True,
#     "font.family": "serif",
#     "font.serif": "Computer Modern",
#     "savefig.bbox": "tight",
#     "savefig.format": "pdf"})
# plt.rc('font', size=16)


### Plot the solution
for i in range(np.shape(fortran_times)[0]):
    # plt.plot(grid,fortran_soln[i,2:-2],'x',linestyle=None)
    plt.plot(grid,fortran_soln[i,2:-2])
    plt.xlabel('x')
    plt.title("t = " +str(fortran_times[i]))
    plt.ylim([-1.5,10])
    plt.draw()
    plt.pause(0.01)
    plt.cla()
