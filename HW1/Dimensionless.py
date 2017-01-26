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

def Euler(y0, nsteps):
    deltime = 1/nsteps #length of a time step
    newtimes = np.arange(0,1,deltime) #Generates n time points
    newvals = []
    newvals.append(y0) #Starting point
    for i in range(nsteps-1):
        nextval = newvals[i]*(1 - deltime) #Note factoring
        newvals.append(nextval) #Add new temperature to Euler array
    return newtimes, newvals

#Main program    
y0 = initial()
nsteps = 2501
times, vals = Euler(y0, nsteps) #The Euler method
plt.plot(times, vals) #Make the plot
plt.xlabel('x')
plt.ylabel('y value')
plt.show()   
print('For ' + str(nsteps) + ' time steps, value is ' + str((vals[-1]/m.exp(-1) - 1)*100) + '% off')

'''
For 101 time steps, Value is 0.4979270251807666% off
For 1001 time steps, Value is 0.0499791770809388% off
For 2501 time steps, value is 0.019996667330324236% off

Note that our error is linear with delta t, thus Euler is first order in time.


Problem 3

This can be seen by noting that the Euler method involves two terms.
First, we take the value of the function at that moment. 
Then we subtract the product of the rate of change at that moment, and the time step to our next value.
This is equivalent to the first two terms of the Taylor series:

f(x) = f(x0) + f'(x0)*(x-x0)
'''