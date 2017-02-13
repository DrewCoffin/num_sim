# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 22:12:25 2017

Plotting Laplace Equation
Drew Coffin
February 1, 2017
Python 3.5
"""
#Necessary modules
import matplotlib.pyplot as plt
import numpy as np
import math as m
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import Axes3D

#Calculate the leading coefficient for the nth term
def coeff(n):
    c = 4*To / (n*np.pi*m.sinh(n*np.pi))
    return c

#What resolution do you want on the individual hyperbolic curves?
def xres(n):
    dx = 1 #Width of the square
    nsteps = n #number of steps to cross the square
    return dx/nsteps
    
#Calculates the nth term in the Fourier series as an array of Z-values
def indivterm(X,Y,n): 
    c = coeff(n)
    Z = c * np.sin(n*m.pi*X) * np.sinh(n*m.pi*(1-Y))
    return Z
    
'''How many terms in the infinite series do you want to include?'''
nterms = 100 
To = 10 #Setting some value for To
A_n = coeff(nterms)
'''How many grid points do you want in each dimension?'''
step = xres(100) 

# Make data.
X = np.arange(0, 1, step) #Initialize gridpoints
Y = np.arange(0, 1, step)
X, Y = np.meshgrid(X, Y) #Make meshgrid
Z = indivterm(X,Y,1) #First term in infinite series
for i in range(3, nterms, 2):#Successive terms
    Z += indivterm(X,Y,i)

    
#Set up the plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.view_init(azim=30) #Viewing angle in degrees (change to rotate)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
                       
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
