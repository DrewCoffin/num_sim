# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 21:16:45 2017

@author: Drew Coffin
"""

'''Note the electric field is q*E/m, which is units of acceleration.
The magnetic field is q*B/m, which is units of (gyro)frequency.
'''

import numpy as np
import matplotlib.pyplot as plt

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
E0 = 0
v0 = np.array([1,0,0])
Tgyro = 2*np.pi/B0[2]
Tfin = 10*Tgyro
delT = Tgyro*0.01
ntsteps = int(Tfin/delT)

pos0 = np.array([0,1,0])

i = 0
pos = [float(pos0[j]) for j in range(len(pos0))] #Make sure initial position is in floats
xarr = [pos[0]]
yarr = [pos[1]]
Bnew = np.zeros(3)
vel = v0
while i < ntsteps:
    Bnew[2] = B0[2]*pos[0]
    vp = vplus(vel,Bnew,delT)
    vn = vnew(vp,delT,E0)
    pos += delT*vn
    vel = vn #Update the velocity
    xarr.append(pos[0])
    yarr.append(pos[1])
    #if i/1000 - int(i/1000) == 0:
     #   print(i) #tracker
    i += 1

plt.plot(xarr, yarr, 'b')
plt.axis('equal')
plt.show()