# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 14:27:09 2020

@author: Baihua Wu
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt



"""
g=np.empty([2,2],dtype=np.float64)
#g=[[1,2],[4,3]]
g[0,0]=1
g[0,1]=2
g[1,0]=4
g[1,1]=3
e,v=np.linalg.eig(g)
print(e)
print(v)

"""



def potential( itype , x ):
    if itype==1:
        U=0.5*(x-1)**2
    elif itype==2:
        U=0.25*x**4
    #elif itype==3:
    return U

def dpotential( itype , x ):
    if itype==1:
        U=x-1
    elif itype==2:
        U=x**3
    #elif itype==3:
    return U

itype=1
x_0=-30.0
x_e=30.0

mass=1.0
dx=0.1

Ngrid=int((x_e-x_0)/dx)
#print(Ngrid)

#psi=np.empty([Ngrid,Ngrid],dtype=np.float64)
T=np.zeros([Ngrid,Ngrid],dtype=np.float64)
V=np.zeros([Ngrid,Ngrid],dtype=np.float64)
x=np.zeros([Ngrid],dtype=np.float64)
psi=np.zeros([Ngrid,Ngrid],dtype=np.float64)
e=np.zeros([Ngrid],dtype=np.float64)

for i in range(0,Ngrid):
    x[i]=x_0+i*dx
    V[i,i]=potential( itype , x[i] )
    T[i,i]=2/(2*mass*dx**2)
    if i>0:
        T[i,i-1]=-1/(2*mass*dx**2)
    
        T[i-1,i]=-1/(2*mass*dx**2)

#print(T[0])
#print(V[0])

H=T+V
#print(H)
eig,eigp=np.linalg.eig(H)


index=np.argsort(eig)
for i in range(Ngrid):
    e[i]=eig[index[i]]
    psi[i]=eigp[:,index[i]]
    
"""    
print(e)
print(psi[1])

plt.figure(1)
for i in range(0,4):
    plt.plot(x,psi[i])
#plt.plot(x,potential( itype , x ))
plt.show()
"""

phi0=np.zeros([Ngrid],dtype=np.float64)
alpha=0.1
x0=4
for i in range(Ngrid):
    phi0[i]=(alpha/np.pi)**0.25*np.exp(-0.5*alpha*(x[i]-x0)**2)

coeff=np.zeros([Ngrid],dtype=np.float64)
for i in range(Ngrid):
    for j in range(Ngrid):
        coeff[i]=coeff[i]+phi0[j]*psi[i,j]*dx
t=np.pi/2      
phit=np.zeros([Ngrid,2],dtype=np.float64)   
for i in range(Ngrid):
    for j in range(Ngrid):
        phit[i,0]=phit[i,0]+np.cos(e[j]*t)*coeff[j]*psi[j,i]
        phit[i,1]=phit[i,1]-np.sin(e[j]*t)*coeff[j]*psi[j,i]
    
"""
plt.figure(1)

plt.plot(x,phit[:,0]**2+phit[:,1]**2)
#plt.plot(x,potential( itype , x ))
plt.show()
"""

A=np.zeros([Ngrid],dtype=np.float64)
A=np.sqrt(phit[:,0]**2+phit[:,1]**2)
d2x=np.zeros([Ngrid,Ngrid],dtype=np.float64)
for i in range(0,Ngrid):
    d2x[i,i]=potential( itype , x[i] )
    
    if i>0:
        d2x[i,i-1]=1
        d2x[i-1,i]=1


Qpot=np.zeros([Ngrid],dtype=np.float64)
Qpot=-0.5/(dx**2)*np.matmul(d2x,A)/A

plt.figure(1)
plt.plot(x,Qpot)
plt.show()

"""
d1x=np.zeros([Ngrid,Ngrid],dtype=np.float64)
for i in range(Ngrid):
    
    d1x[i,i]=-1
    if i>0:
        d1x[i-1,i]=1
        
    
dQpot=np.zeros([Ngrid],dtype=np.float64)        
dQpot=1/(dx)*np.matmul(d1x,Qpot)
"""
plt.figure(1)
plt.plot(x,dQpot)
plt.show()


t_start=0
t_end=4*np.pi
dt=0.1
q=4.0
p=0.0

nstep=int((t_end-t_start)/dt)
x_traj=np.zeros([nstep],dtype=np.float64)
time=np.zeros([nstep],dtype=np.float64)
for i in range(nstep):
    x_traj[i]=q
    time[i]=i*dt
    for j in range(Ngrid):
        if x[j]<q and x[j+1]>q:
            p=p-(dpotential(itype,q)+dQpot[j])*dt/2
            break

    q=q+p*dt
    
    for j in range(Ngrid):
        if x[j]<q and x[j+1]>q:
            p=p-(dpotential(itype,q)+dQpot[j])*dt/2
            break

plt.figure(1)
plt.plot(time,x_traj)
plt.show()
