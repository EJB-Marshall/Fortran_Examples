#!/usr/bin/env python
# encoding: utf-8

"""
Euler_DS.py

Code to evolve the Relativistic Euler Equations
on a exponentially expanding de Sitter background.

Re-coded version of old scheme used in 2022 paper. 
Created to make new plots for thesis.
"""


### Import Python Libraries
import sys
import numpy as np
import h5py
import argparse
import time
import matplotlib.pyplot as plt
import time
import math


t_start = time.time()


#######################################
#Parser Settings
#######################################

# Initialise parser
parser = argparse.ArgumentParser(description=\
"""This program numerically solves an IVP for the
relativistic Euler Equations in T2 Symmetry.""")

# Parse files
parser.add_argument('-d','-display', default = False, 
    action='store_true', help=\
"""A flag to indicate if visual display is required.""")
parser.add_argument('-N', type=int,help=\
"""The number of grid points.""")
parser.add_argument('-f','-file', help=\
"""The name of the hdf file to be produced.""")
args = parser.parse_args()

#Output Settings
display_output = args.d
store_output = args.f is not None
if store_output and args.f is None:
    print("Euler_DS.py: error: argument -f/-file is required")
    sys.exit(1)

#Check Inputs
if args.N is None:
    print("Euler_DS.py: error: argument -N is required.")
    sys.exit(1)

    

##############################################################
# Finite Difference Derivative Operators
##############################################################  
# @njit(cache=True) 
def Dx(f,r):
    """ This function calculates the first derivative
        using a 2nd order central finite difference stencil"""
    h = r[1]
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    #2nd Order Central Finite Difference
    # df[1:n-1] = (-f[0:n-2]+f[2:n])/(2*h)

    # #Periodic Boundary Condition
    # df[0] = (-f[n-1]+f[1])/(2*h)
    # df[n-1] = (-f[n-2]+f[0])/(2*h)

    #4th Order Central Finite Difference
    df[2:n-2] = (f[0:n-4]-8*f[1:n-3]+8*f[3:n-1]-f[4:n])/(12*h)

    #Periodic Boundary Condition
    df[1] = (f[n-1]-8*f[0]+8*f[2]-f[3])/(12*h)
    df[0] = (f[n-2]-8*f[n-1]+8*f[1]-f[2])/(12*h)
    df[n-2] = (f[n-4]-8*f[n-3]+8*f[n-1]-f[0])/(12*h)
    df[n-1] = (f[n-3]-8*f[n-2]+8*f[0]-f[1])/(12*h)


    return df


##################################################################
# SOLVE PDES
##################################################################


def advec(f,r,t):

    """This function computes the evolution equations.
    f - The functions we are evolving
    r - This is the grid variable (not the ratio of slopes!)
    t - Time
    K - Sound Speed"""

    z = f[0]
    

    ######################################################################
    # Collect Evolved Quantities
    ######################################################################
    

    return np.array([-Dx(z,r)])


###########################################################################
# Runge-Kutta Scheme
##########################################################################
def rk2(f,tspan,yinit,dt):
    Nt = np.shape(tspan)[0]
    Nvar = np.shape(yinit)[0] #Number of variables to solve for
    Nx = np.shape(yinit)[1]
    y = np.zeros((Nt,Nvar,Nx))
    t = tspan
    y[0,:,:] = yinit[:]
    t0 = tspan[0]
    y0 = yinit[:]
    for i in range(1,Nt):
        deltaT = tspan[i]-tspan[i-1]; #Time steps to output
        M = math.ceil(abs(deltaT/dt)); #Number of time steps in between
        h = deltaT/M; #Length of time step
        for l in range(0,M):
            k1 = f(t0,y0)
            k2 = f(t0+h,y0+h*k1)
            y1 = y0 + h*(0.5*k1 + 0.5*k2)
            y0 = y1
            t0 = t0 + h
        y[i,:,:] = y1
    return [t,y]

def rk(f,t0,tf,yinit,dt,x_step):
    """Strong Stability preserving 3rd order Runge-Kutta"""
    y = []
    t = []
    y.append(yinit[:])
    t.append(t0)
    y0 = yinit[:]
    h = dt
    i = 0
    while t0<tf:
            i += 1 

            ### 'Classic' RK4 Method
            k1 = f(t0,y0)
            k2 = f(t0+0.5*h,y0+0.5*h*k1)
            k3 = f(t0 +0.5*h,y0+0.5*h*k2)
            k4 = f(t0+h,y0+h*k3)
            y0 = y0 + h*(1/6*k1 + 1/3*k2 +1/3*k3 +1/6*k4)

            t0 = t0 + h
            h = 0.5*(x_step) ### Adjust timestep based on characteristic speeds
            if h > 0.5*x_step:
                h = 0.5*x_step

            if i == 100: ### Timesteps to output
                i = 0
                y.append(y0)
                t.append(t0)
    return [np.array(t),np.array(y)]




    
############################################################################
#Creating grid and initial data
############################################################################
tspan = np.linspace(0,20,num=500) #300 time steps to output
# tspan = np.flip(tspan) 
R = 2*np.pi #Length of grid
N = args.N #Number of Grid Points 
r = np.delete(np.linspace(0,2*np.pi,num=N+1),-1) #Uniform Grid
h = r[1] #Grid Spacing

##########################
# Initial Data
##########################

#Parameters for A0, A etc.
a = 0.1

# z = 0.3*np.exp(-10*(r-np.pi)**2) 
# z[np.abs(z)<1e-15] = 0
z = np.sin(r)


y0 = np.array([z])


######## Evolution Equation ############
rhs = lambda t,y: advec(y,r,t)
# [t,y] = rk2(rhs, tspan, y0, 0.1*(r[1])**(4/3))
[t,y] = rk(rhs,tspan[0],20.0,y0,0.5*h,h)

t_stop = time.time()
print("Elapsed time is " + str(t_stop-t_start))


# for i in range(np.shape(t)[0]):
#     plt.plot(y[i,0,:])
#     plt.draw()
#     plt.pause(0.01)
#     plt.cla()









            
            



