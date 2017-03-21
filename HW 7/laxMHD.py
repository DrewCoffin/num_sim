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
    z0 = 0.25
    u0 = 1
    vinit = np.zeros(shape=(3,len(z)))
    vinit[0,:] = u0*np.exp(-(z[:]/z0)**2) #Note x is 0th row
    return vinit

'''Lax-Wendroff.

Note I am doing two steps. 
First, laxhalf calculates the "half-step" value.
Then, combined with the calculated "half-step" flux, the full step value is obtained.'''
def laxhalf(C,i,old,flux):
    half = np.zeros(len(old))
    half[i] = 0.5*(old[i-1]+old[i]) - 0.5*C*(flux[i] + -1*flux[i-1])
    if i == len(old) - 1:
        half[i] = 0.5*(old[-1]+old[0]) - 0.5*C*(flux[0] + -1*flux[-1]) #Periodic boundary
    return half
    
def lax(C,i,old,half,halfflux):
    new = np.zeros(len(half))
    #print(newflux)
    new[i] = old[i] - C * (halfflux[i] + -1*halfflux[i-1])
    if i == len(half) - 1:
        new[i] = old[0] - C * (halfflux[0] + -1*halfflux[-1]) #Periodic boundary
    return new
    
'''Quick transposing of an array, a'''
def star(a):
    return np.array(list(map(list,zip(*a))))
    
'''Energy expression'''
def nrg(i,rho, vel, B0, p):
    gamma = 0.5
    velmag = np.linalg.norm(vel[:,i])
    Bmag = np.linalg.norm(B0[:,i]) #Magnitude of B
    termone = 0.5*rho[i]*velmag**2
    termtwo = 0.5*Bmag**2
    termthree = p/(2*gamma-2)
    return termone + termtwo + termthree
    
'''Here, the fluxes are calculated from each expression.'''

def Fmass(i,dens,vel):
    massflux = np.array(dens[i]*vel[:,i])
    return massflux
    
def Fmom(i, dens,vel,press,B0):
    Bmag = np.linalg.norm(B0[:,i]) #Magnitude of B
    vtrans = star([vel[:,i]])
    termone = dens*np.dot(vtrans,[vel[:,i]])
    termtwo = 0.5*(press[0]+Bmag)*np.eye(3)
    termthree = np.dot(star([B0[:,i]]),[B0[:,i]])
    final = termone + termtwo - termthree
    return final
    
def Fmag(i,B0,vel):
    Btrans = star([B0[:,i]])
    vtrans = star([vel[:,i]])
    magflux = np.dot(vtrans,[B0[:,i]]) - np.dot(Btrans,[vel[:,i]]) 
    return magflux
    
def Fnrg(i,rho,vel,B0,p):
    Bmag = np.linalg.norm(B0[:,i]) #Magnitude of B
    w = nrg(i,rho, vel, B0, p)
    vtrans = star([vel[:,i]])
    termone = w + 0.5*(p[i]+Bmag**2)    
    termtwo = np.dot(vtrans,B0[:,1])*B0[:,1]
    return termone*vel[:,1] - termtwo
    
def nrg2press(nrg,dens,vel,B0):
    gamma = 0.5
    Bmag = term = p = np.zeros(len(dens))
    for i in range(len(B0)):
        Bmag[i] = np.linalg.norm(B0[:,i]) #Magnitude of B
        term[i] = nrg[i] - 0.5*(dens[i]*vel[0,i]**2-Bmag[i]**2)
        p[i] = 2*(gamma-1)*term[i]
    return p
    
# Grid resolution
delt = 0.05
zstep = 0.05
nzsteps = 1/zstep
ntsteps = 1/delt
c = delt/zstep

'''Generate initial conditions'''
zarr = np.arange(-1,1+zstep,zstep)
vel = initial(zarr)
dens = np.ones(len(zarr))
press = np.ones(len(zarr))
#B only in z direction
mag = np.zeros(shape=(3,len(zarr)))
mag[2] = np.ones(len(zarr))
mom = np.zeros(shape=(3,len(zarr)))
mom[0] = [dens[i]*vel[0,i] for i in range(len(zarr))]
nrgarr = np.ones(len(zarr))
rhohalf = nrghalf = np.zeros(len(zarr))
momhalf = maghalf = halfvel = np.zeros(shape=(3,len(zarr)))

'''The iterative loop to calculate all components.'''
for i in range(len(zarr)):
    '''Initial fluxes, turned into column vectors or tensors per point'''
    massflux = star([Fmass(i, dens,vel)]) #vector
    momflux = Fmom(i,dens[i],vel,press,mag) #tensor
    Zmomflux = star([momflux[2]]) #Del is z-derivative wrt mag field, so z dot tensor
    magflux = Fmag(i, mag,vel)
    Zmagflux = star([magflux[2]])
    nrgflux = star([Fnrg(i, dens,vel,mag,press)])
    '''Half-step values'''
    rhohalf[i] = laxhalf(c,i,dens,massflux)
    momhalf[:,i] = laxhalf(c,i,mom[0],Zmomflux)
    maghalf[:,i] = laxhalf(c,i,mag[2],magflux)
    nrghalf[i] = laxhalf(c,i,nrgarr,nrgflux)
    halfvel[:,i] = [momhalf[k]/rhohalf[k] for k in range(len(momhalf))]
    '''Half-step fluxes'''
    halfFmass = Fmass(i,rhohalf,halfvel)
    halfFmom = Fmom(i,rhohalf,halfvel,press,maghalf)
    halfFmag = Fmag(i,maghalf,halfvel)
    halfFnrg = Fnrg(i,rhohalf,halfvel,maghalf,press)
    '''New values'''
    newrho = lax(c,i,dens,rhohalf,halfFmass)
    newmom = lax(c,i,mom[0],momhalf,halfFmom)
    newmag = np.zeros(shape=(3,len(zarr)))
    newmag[2] = lax(c,i,mag[2],maghalf[2],halfFmag)
    newnrg = lax(c,i,nrgarr,nrghalf,halfFnrg)
    newvel = np.zeros(shape=(3,len(zarr)))
    newvel[0] = [newmom[k]/newrho[k] for k in range(len(newmom))]
    newpress = nrg2press(newnrg,newrho,newvel,newmag)
    '''Make plots for each parameter'''
    if int(i%5) == 0: #Plot every fifth step
        plt.figure(1)
        plt.subplot(411)
        plt.plot(zarr,vel[0], 'r--', zarr, newvel[0], 'b--')
        plt.title('FTCS scheme (unstable)')

        plt.subplot(412)
        plt.plot(zarr,mag[2], 'r--', zarr, newmag[2], 'b--')
        plt.title('Upwind scheme')

        plt.subplot(413)
        plt.plot(zarr,dens, 'r--', zarr, newrho, 'b--')
        plt.title('Leapfrog scheme')

        plt.subplot(414)
        plt.plot(zarr,press, 'r--', zarr, newpress, 'b--')
        plt.title('Lax-Wendroff scheme')

        plt.subplots_adjust(top = 3) #Make titles spacing behave
        plt.show()
    i += 1
    

plt.show()