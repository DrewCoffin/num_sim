"""
Solving Dedimensionalized equation via Euler Method
Drew Coffin
January 25, 2017
Python 3.5
"""
#Necessary modules
import matplotlib.pyplot as plt
import numpy as np
import math as m

#Various constants  
def initial():
    start = 1 #initial value
    return start
    
#Is this a Modified or Improved Euler?
def coeffs(flag):
    if flag.strip() == 'mod':
        alpha = 0.5
        gamma1 = 0 #Clear any old value
        gamma2 = 1 
    elif flag.strip() == 'imp':
        alpha = 1
        gamma1 = 0.5
        gamma2 = 0.5
    else: #Variable coefficients
        alpha = 0
        gamma1 = 1
        gamma2 = 0
    return alpha, gamma1, gamma2
    
def Eulerstep(start, a, deltime):
    nextval = start*(1 - a*deltime)
    return nextval

def SecondEuler(select, y0, nsteps):
    a,g1,g2 = coeffs(select)
    deltime = 1/nsteps #length of a time step
    newtimes = np.arange(0,1,deltime) #Generates n time points
    newvals = []
    newvals.append(y0) #Starting point
    for i in range(nsteps-1):
        partstep = Eulerstep(newvals[i], a, deltime) #Partial time step
        nextval = newvals[i] - deltime*(g1*newvals[i] + g2*partstep)
        #Note the negative sign is due to F(x) = -f(x)
        newvals.append(nextval) #Add new value to Euler array
    return newtimes, newvals

#Main program    
y0 = initial()
nsteps = 200
times, vals = SecondEuler('mod', y0, nsteps) #The Euler method
plt.plot(times, vals) #Make the plot
plt.xlabel('delta x')
plt.ylabel('y value')
plt.show()   
print('For ' + str(nsteps) + ' time steps, value is ' + str((vals[-1]/m.exp(-1) - 1)*100) + '% off')

'''
For 100 time steps, value is 1.0066958545162041% off
For 200 time steps, value is 0.5016703138565948% off

This is first order, so something is wrong.
'''