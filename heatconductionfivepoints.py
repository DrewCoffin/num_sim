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
import pylab
    
'''The sum inside our exact expression'''
def sumcalc(xn,t,alpha):
    '''How many terms in the infinite series do you want to include?'''
    nterms = 1000   
    value = 0 #Start
    for j in range(nterms):
        constant = 400/((2*j-1)*m.pi)
        sinterm = m.sin((2*j-1)*m.pi*xn)
        expterm = m.exp(-alpha*(2*j-1)**2*m.pi**2*t)
        value += constant*sinterm*expterm
    return value

'''Returning the exact solution'''    
def exact(x,t,alpha):
    exactarr = []
    for i in range(len(x)):
        exactarr.append(100-sumcalc(x[i],t,alpha))
    return exactarr

'''Making the next time step'''
def nexttime(s,oldarr):
    nexttime = [100] #Endcap
    firstnum = 11*oldarr[0]/12 - 5*oldarr[1]/3 + 0.5*oldarr[2] + oldarr[3]/3 - oldarr[4]/12
    nexttime.append(oldarr[1]+s*firstnum) #Note time term, oldarr[i]
    for i in range(2, len(oldarr)-2):
        '''Five point symmetric calculation'''
        nextnum = -oldarr[i-2]/12 + 4*oldarr[i-1]/3 - 5*oldarr[i]/2 + 4*oldarr[i+1]/3 - oldarr[i+2]/12
        '''Note time term, oldarr[i]'''        
        nexttime.append(oldarr[i]+s*nextnum) #Note time term, oldarr[i]
    lastnum = 11*oldarr[99]/12 - 5*oldarr[98]/3 + 0.5*oldarr[97] + oldarr[96]/3 - oldarr[95]/12
    nexttime.append(oldarr[99]+s*lastnum)
    nexttime.append(100) #Other end cap
    return nexttime

'''What value of s?'''
s = 0.1
alpha = 10**-5
delx = 0.01
deltime = s*delx**2/alpha #Explicitly calculate timestep
nsteps = int(1/delx)
finaltime = 10000 #User-set final time
ntime = finaltime/deltime

# Make data.
xarr = np.arange(0,1,delx)
exactarr = exact(xarr,finaltime*1.618,alpha) #Mysterious golden ratio on the final time
initarr = [100] + (nsteps-2)*[0] + [100]
newarr = nexttime(s,initarr) #First time step
i = 1
while i < ntime:
    newarr = nexttime(s,newarr)
    if i/1000 - int(i/1000) == 0: #Plot every 1000th time step
        print(i) #Tracking on script progress
        plt.plot(xarr,newarr, 'r--', xarr, exactarr, 'b--')
    i += 1
print('Blue is exact solution at final time step')
print('Red is approximate solution at every 1000th time step')

