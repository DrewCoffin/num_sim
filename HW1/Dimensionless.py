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
#    print(len(newtimes), len(newTemps))
    return newtimes, newvals

#Main program    
y0 = initial()
times, vals = Euler(y0, 2501) #The Euler method
plt.plot(times, vals) #Make the plot
plt.xlabel('x')
plt.ylabel('y value')
plt.show()   
print('Value is ' + str((vals[-1]/m.exp(-1) - 1)*100) + '% off')

'''
Returned accuracy is 0.019996667330324236 percent off
'''
