"""
Solving Newton Cooling via Euler Method
Drew Coffin
January 25, 2017
Python 3.5
"""
#Necessary modules
import matplotlib.pyplot as plt
import numpy as np

#The data set we wish to match
def GetArrays():
    times = list(range(16)) #0,1,2...15 array. Note Python starts at 0
    Temps = [83, 77.7, 75.1, 73, 71.1, 69.4, 67.8, 66.4, 64.7, 63.4, 62.1, 61, 59.9, 58.7, 57.8, 56.6]
    return times, Temps
  
#Various constants  
def constants(times, Temps):
    T0 = Temps[0] #initial temperature
    tf = times[-1] #Gets final time value
    Tair = 25 #ambient temperature
    R = 0.042 #rate of cooling
    return T0, tf, Tair, R

def Euler(tf, T0, R, nsteps):
    deltime = tf/nsteps #length of a time step
    newtimes = np.arange(0,15,deltime) #Generates n time points
    newTemps = []
    newTemps.append(T0) #Same start point
    for i in range(nsteps-1):
        nextTemp = newTemps[i] - R*(newTemps[i] -Tair)*deltime
        newTemps.append(nextTemp) #Add new temperature to Euler array
#    print(len(newtimes), len(newTemps))
    return newtimes, newTemps

#Main program
times, Temps = GetArrays()
T0, tf, Tair, R = constants(times, Temps)
Eulertimes, EulerTemps = Euler(tf, T0, R, 100) #The Euler method
plt.plot(times, Temps, Eulertimes, EulerTemps) #Make the plot
plt.xlabel('time (min)')
plt.ylabel('Temp (C)')
plt.show()    

'''
Part C:
Estimate M = 0.300 kg (about 10 ounces).
Estimate A = 25 square centimeters (.0025 square meters)

In[56]: 0.0237*.0075/(.042*.3*4.19*10**3)
Out[56]: 3.3668598704398217e-06 m

Part D
'''