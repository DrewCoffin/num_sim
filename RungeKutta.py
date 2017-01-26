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
    if flag.strip() == 'imp':
        alpha = 1
        gamma1 = 0.5
        gamma2 = 0.5
    else: #Variable coefficients
        alpha = 0.7
        gamma1 = 0.7
        gamma2 = 0.7
    return alpha, gamma1, gamma2

def ModEuler(y0, nsteps):
    a,g1,g2 = coeffs('mod')
    deltime = 1/nsteps #length of a time step
    newtimes = np.arange(0,1,deltime) #Generates n time points
    newvals = []
    newvals.append(y0) #Starting point
    for i in range(nsteps-1):
        nextval = newvals[i]*(1 - deltime) #Note factoring
        newvals.append(nextval) #Add new temperature to Euler array
#    print(len(newtimes), len(newTemps))
    return newtimes, newvals

#Main program    
y0 = initial()
times, vals = ModEuler(y0, 100) #The Euler method
plt.plot(times, vals) #Make the plot
plt.xlabel('delta x')
plt.ylabel('y value')
plt.show()   
print((vals[-1]/m.exp(-1) - 1)*100) #How accurate am I?

'''
Returned accuracy is percent off
'''