import numpy as np
from matplotlib import pyplot as plt

python_soln = np.loadtxt('python_soln.txt')
python_times = np.loadtxt('python_soln_times.txt')

fortran_soln = np.loadtxt('fortran_soln.txt')
fortran_times = np.loadtxt('fortran_soln_times.txt')

plt.rcParams.update({"text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern",
    "savefig.bbox": "tight",
    "savefig.format": "pdf"})
plt.rc('font', size=16)

### Plot outputs against the exact solution
plt.plot(python_times,python_soln, 'x',label='Python Solution')
plt.plot(fortran_times,fortran_soln, '.',label='Fortran Solution')
plt.plot(fortran_times, 0.01*np.exp(2*fortran_times), label='Exact Solution')
plt.xlabel("Time")
plt.ylabel("Solution")
plt.legend()
plt.show()