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
def laxhalf(C,old,flux):
    half = np.zeros(len(old))
    half[:-1] = 0.5*(old[:-1]+old[1:]) - 0.5*C*(flux[1:] + -1*flux[:-1])
    half[-1] = 0.5*(old[-1]+old[0]) - 0.5*C*(flux[0] + -1*flux[-1]) #Periodic boundary
    return half
    
def lax(C,old,half,halfvel):
    new = np.zeros(len(half))
    newflux = np.array([half[i]*halfvel[i] for i in range(len(halfvel))])
    #print(newflux)
    new[1:] = old[1:] - C * (newflux[1:] + -1*newflux[:-1])
    new[0] = old[0] - C * (newflux[0] + -1*newflux[-1]) #Periodic boundary
    return new
    
'''Quick transposing'''
def star(a):
    return np.array(map(list,zip(*a)))    
    
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
delt = 0.008
zstep = 0.01
nzsteps = 1/zstep
ntsteps = 1/delt
c = delt/zstep

'''Generate initial conditions'''
zarr = np.arange(-1,1+zstep,zstep)
v0 = initial(zarr)
dens0 = np.ones(len(zarr))
p0 = np.ones(len(zarr))
#B only in z direction
B0 = np.zeros(shape=(3,len(zarr)))
B0[2] = np.ones(len(zarr))
mom0 = np.zeros(shape=(3,len(zarr)))
mom0[0] = [dens0[i]*v0[0,i] for i in range(len(zarr))]
#print(mom0)

'''Do the initial iteration'''
massflux = Fmass(dens0[0],v0)
momflux = Fmom(dens0[0],v0,p0[0],B0)
magflux = Fmag(B0,v0)
rhohalf = laxhalf(c,dens0,massflux)
momhalf = laxhalf(c,mom0[0],momflux)
maghalf = laxhalf(c,B0[2],magflux)
halfvel = np.zeros(shape=(3,len(zarr)))
halfvel[0] = [momhalf[i]/rhohalf[i] for i in range(len(momhalf))]
newrho = lax(c,dens0[0],rhohalf,halfvel[0])
newmom = lax(c,mom0[0],momhalf,halfvel[0])
newmag = lax(c,B0[0],maghalf,halfvel[0])
newvel = [newmom[i]/newrho[i] for i in range(len(newmom))]

'''Mag equation to update B0, then lax full step for new values'''

'''lf1 = leapfrog(c, initarr, laxarr) #First two steps for leapfrog
lf2 = leapfrog(c, laxarr, lf1)'''
    
'''Iterative loop'''
i = 1
while i <= ntsteps: 
    massflux = Fmass(newrho,newvel)
    momflux = Fmom(newrho,newvel,p0[0],newmag)
    magflux = Fmag(newmag,newvel)
    rhohalf = laxhalf(c,newrho,massflux)
    momhalf = laxhalf(c,newmom,momflux)
    maghalf = laxhalf(c,newmag,magflux)
    halfvel = np.zeros(shape=(3,len(zarr)))
    halfvel[0] = [momhalf[i]/rhohalf[i] for i in range(len(momhalf))]
    newrho = lax(c,newrho,rhohalf,halfvel[0])
    newmom = lax(c,newmom,momhalf,halfvel[0])
    newmag = lax(c,newmag,maghalf,halfvel[0])
    newvel = [newmom[i]/newrho[i] for i in range(len(newmom))]
    if int(i%5) == 0:
        #print(i, laxdens)
        plt.plot(zarr, v0[0], 'r--', zarr, newvel[0], 'b--')
        #plt.plot(zarr, laxarr, 'b--')
    i += 1
 
#plt.plot(zarr, initarr, 'r--', zarr, laxarr, 'b--')#                               

plt.show()