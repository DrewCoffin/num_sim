# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 16:49:41 2017

@author: Drew Coffin
"""
import numpy as np
import matplotlib.pyplot as plt
import random as r

'''Make random numbers'''
Rsarr = []
Rtharr = []
nvals = 1000
i = 0
while i < nvals:
    Rsarr.append(r.random())
    Rtharr.append(r.random())
    i += 1
    
col = np.arange(nvals) #Make color values for the points
    
vth = 1
vx = [vth*(-2*np.log(Rsarr[i]))*np.cos(2*np.pi*Rtharr[i]) for i in range(len(Rsarr))]
vy = [vth*(-2*np.log(Rsarr[i]))*np.sin(2*np.pi*Rtharr[i]) for i in range(len(Rsarr))]

plt.scatter(vx, vy, c=col)
plt.axis('equal')
plt.title('Maxwellian distribution')
plt.show()

'''Plotting 10 random particles'''

def vplus(vold,Bold,delT):
    Bmag = np.dot(Bold,Bold)**(0.5)
    alpha = (1 - Bmag*delT**2/4)*vold
    beta = delT*np.cross(vold,Bold)
    gamma = 0.5*delT**2*np.dot(vold,Bold)*Bold
    denom = 1 + Bmag*delT**2/4
    return np.array((alpha+beta+gamma)/denom)
    
def vnew(vplus, delT, Eold):
    return np.array(vplus + delT*Eold)

B0 = np.array([0,0,1])
E0 = np.array([0,0.1,0])
Tgyro = 2*np.pi/B0[2]
Tfin = 10*Tgyro
delT = Tgyro*0.01
ntsteps = int(Tfin/delT)

pos = np.array([0,1.0,0])

i = 0
while i < 10:
    xarr = [pos[0]]
    yarr = [pos[1]]
    vel = np.array([vx[100*i], vy[100*i], 0])
    k = 0
    while k < ntsteps:
        vp = vplus(vel,B0,delT)
        vn = vnew(vp,delT,E0)
        pos += delT*vn
        vel = vn #Update the velocity
        xarr.append(pos[0])
        yarr.append(pos[1])
        k += 1
    plt.plot(xarr,yarr)
    plt.title('Trajectories (cumulative)')
    print("Calculating the ", i, "th particle's trajectory") #Tracker
    i += 1