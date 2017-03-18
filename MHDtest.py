# -*- coding: utf-8 -*-
"""
Plotting MHD Simulation
Drew Coffin
March 21, 2017
Python 3.5
"""
#Necessary modules
import matplotlib.pyplot as plt
import numpy as np
    
'''What was our initial condition?''' 
def initial(z,delx):
    z0 = 0.05
    u0 = rho0 = 1
    vinit = rho0 = np.zeros(shape=(3,len(z)))
    vinit[0,:] = u0*np.exp(-(z[:]/z0)**2) #Note x is 0th row
    return vinit, rho0
    
'''Leapfrog one liner'''
def leapfrog(C,t0, t1):
    t2 = np.zeros(shape=(3,len(t0[2])))
    t2[2,1:-1] = t0[1:-1] - C * (t1[2:] - t1[:-2]) 
    return t2
    
'''Lax-Wendroff one liner'''
def lax(C,dens,vel,old):
    newstar = new = np.zeros(len(old))
    flux = massflux(dens,vel)
    newstar[1:-1] = 0.5*(old[1:-1]+old[2:]) - 0.5*C*(flux[2:] + -1*flux[1:-1])
    new[1:-1] = old[1:-1] - C * (newstar[1:-1]*vel[1:-1] - newstar[:-2]*vel[:-2])
    return newstar
    
'''Quick transposing'''
def star(a):
    return list(map(list,zip(*a)))    
    
'''MHD equations'''
def massflux(dens,vel):
    return np.array([dens*vel[0][i] for i in range(len(vel[0]))])
    
def momflux(dens,vel,press,B0):
    Bmag = np.linalg.norm(B0[:,0]) #Magnitude of B
    vtrans = star([vel[0]])
    termone = dens*np.dot(vtrans,vel[0])
    termtwo = 0.5*(press+Bmag)*np.eye(len(vel[0]))
    termthree = np.dot(star(B0[2]),B0[2])
    final = termone + termtwo - termthree
    return final
    
def magflux(B0,vel):
    Btrans = star([B0[2]])
    vtrans = star([vel[0]])
    return np.dot(Btrans,vel[0]) - np.dot(vtrans,B0[2])

# Grid resolution
delt = 0.008
zstep = 0.01
nzsteps = 1/zstep
ntsteps = 1/delt
c = delt/zstep

'''Generate initial conditions'''
zarr = np.arange(0,1+zstep,zstep)
initarr = initial(zarr,zstep)
dens0 = np.ones(len(zarr))
B0 = [np.zeros(len(zarr)), np.zeros(len(zarr)), np.ones(len(zarr))] #B only in z direction
#print(initarr[2], len(initarr[2]))

'''Do the initial steps for each method to overwrite initarr'''
laxarr = lax(c, dens0[0], initarr[0], dens0)
'''lf1 = leapfrog(c, initarr, laxarr) #First two steps for leapfrog
lf2 = leapfrog(c, laxarr, lf1)'''
    
'''Iterative loop
i = 1
while i <= ntsteps: 
    leapfrogarr = leapfrog(c, lf1, lf2)
    [lf1, lf2] = [lf2, leapfrogarr] #Update multiple time steps
    laxarr = lax(c,laxarr)
    if int(i%25) == 0:
        plt.plot(zarr, initarr[2], 'r--', zarr, laxarr, 'b--')
    i += 1
 
#plt.plot(zarr, initarr, 'r--', zarr, laxarr, 'b--')
                               
plt.show()'''
