import numpy as np
from matplotlib import pyplot as plt

fortran_soln = np.loadtxt('advection_soln.dat')
fortran_times = np.loadtxt('advection_soln_times.dat')
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
    plt.plot(grid,fortran_soln[i,2:-2],'x',linestyle=None)
    plt.xlabel('x')
    plt.title("t = " +str(fortran_times[i]))
    plt.ylim([-1.5,1.5])
    plt.draw()
    plt.pause(0.01)
    plt.cla()
