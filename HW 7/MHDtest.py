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
    
'''What was our initial velocity?''' 
def initial(z):
    z0 = 0.05
    u0 = 1
    vinit = np.zeros(shape=(3,len(z)))
    vinit[0,:] = u0*np.exp(-(z[:]/z0)**2) #Note x is 0th row
    return vinit
    
'''Leapfrog'''
def leapfrog(C,t0, t1):
    t2 = np.zeros(shape=(3,len(t0[2])))
    t2[2,1:-1] = t0[1:-1] - C * (t1[2:] - t1[:-2]) 
    return t2
    
'''Lax-Wendroff'''
def lax(C,dens,vel,old,flux):
    half = new = np.ones(len(old))
    half[1:-1] = 0.5*(old[1:-1]+old[2:]) - 0.5*C*(flux[2:] + -1*flux[1:-1])
    newflux = np.array([half[i]*vel[0,i] for i in range(len(half))])
    #print(half[90:110],vel[0,90:110])
    new[1:-1] = old[1:-1] - C * (newflux[2:] + -1*newflux[1:-1])
    return new, newflux
    
'''Quick transposing'''
def star(a):
    return list(map(list,zip(*a)))    
    
'''MHD equations'''
def Fmass(dens,vel):
    return np.array([dens*vel[0][i] for i in range(len(vel[0]))])
    
def Fmom(dens,vel,press,B0):
    print(B0)
    Bmag = np.linalg.norm(B0[:,0]) #Magnitude of B
    vtrans = star([vel[0]])
    print(vtrans, vel[0])
    termone = dens*np.dot(vtrans,vel[0])
    termtwo = 0.5*(press[0]+Bmag)*np.eye(len(vel[0]))
    termthree = np.dot(star(B0[2]),B0[2])
    final = termone + termtwo - termthree
    return final
    
def Fmag(B0,vel):
    Btrans = star([B0[2]])
    vtrans = star([vel[0]])
    return np.dot(vtrans,B0[2]) - np.dot(Btrans,vel[0]) 

# Grid resolution
delt = 0.08
zstep = 0.01
nzsteps = 1/zstep
ntsteps = 1/delt
c = delt/zstep

'''Generate initial conditions'''
zarr = np.arange(-1,1+zstep,zstep)
v0 = initial(zarr)
dens0 = p0 = np.ones(len(zarr))
#B only in z direction
B0 = np.zeros(shape=(3,len(zarr)))
B0[2] = np.ones(len(zarr))
mom0 = [dens0[i]*v0[0,i] for i in range(len(zarr))]

'''Do the initial step for the continuity equation'''
massflux = Fmass(dens0[0],v0)
laxdens,newmassflux = lax(c, dens0[0], v0, dens0,massflux)
newvel = np.zeros(shape=(3,len(zarr)))
newvel[0] = [newmassflux[i]/laxdens[i] for i in range(len(laxdens))]

'''Initial step for momentum equation'''
momflux = Fmom(dens0[0],v0,p0,B0)
laxmom,newmomflux = lax(c, dens0[0], v0, mom0, momflux)

'''Initial step for Faraday's Law'''
magflux = Fmag(B0,v0)
laxmag,newmagflux = lax(c, dens0[0], v0, B0[2],magflux)

'''lf1 = leapfrog(c, initarr, laxarr) #First two steps for leapfrog
lf2 = leapfrog(c, laxarr, lf1)'''
    
'''Iterative loop'''
i = 1
while i <= ntsteps: 
    #leapfrogarr = leapfrog(c, lf1, lf2)
    #[lf1, lf2] = [lf2, leapfrogarr] #Update multiple time steps
    laxdens, newmassflux = lax(c,laxdens[0], newvel, laxdens,newmassflux)
    newvel[0] = [newmassflux[i]/laxdens[i] for i in range(len(laxdens))]
    if int(i%10) == 0:
        #print(i, laxdens)
        plt.plot(zarr, v0[0], 'r--', zarr, newvel[0], 'b--')
        #plt.plot(zarr, laxarr, 'b--')
    i += 1
 
#plt.plot(zarr, initarr, 'r--', zarr, laxarr, 'b--')#                               

plt.show()