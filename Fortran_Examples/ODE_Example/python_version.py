import numpy as np
import time
t_start = time.time()
### Heun's Method
def rk2(y0,t0,dt,rhs):
    k1 = rhs(t0,y0)
    k2 = rhs(t0 + dt, y0 + dt*k1)
    y0 = y0 + (dt/2)*(k1 + k2)
    return y0

### Exponential Growth/Decay ODE
def ODE(t,y):
    dydt = 2*y
    return dydt

### Define initial parameters
t0 = 0.0
tend = 5.0
dt = 0.01
y0 = 0.01

### Create list to store solution
y_output = []
y_output.append(y0)
t_output = []
t_output.append(t0)

### Main integration loop
while t0 < tend:
    y0 = rk2(y0,t0,dt,ODE)
    t0 += dt
    y_output.append(y0)
    t_output.append(t0)

### Save solution as a text file
y_soln = np.array(y_output) # Convert list to numpy arrays
t_soln = np.array(t_output)
np.savetxt('python_soln.txt', y_soln)
np.savetxt('python_soln_times.txt', t_soln)

t_end = time.time()
print("Elapsed time is",t_end-t_start) #Prints simulation length of time 