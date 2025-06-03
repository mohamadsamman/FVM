#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 09:20:51 2024

@author: masssonr

resolution par la methode VF de l'equation elliptique 1D

- u''(x) = f(x) sur (0,L)
u(0) = uD
-u'(L) = g

"""

import numpy as np
import matplotlib.pyplot as plt

# longueur du domaine 
L=1.0 

def u(x):
    return np.exp(np.sin(np.pi*x))  

def up(x): 
    return np.pi*np.cos(np.pi*x)*np.exp(np.sin(np.pi*x))  

def f(x):
    s =  np.pi**2*np.exp(np.sin(np.pi*x))
    s = s*(np.sin(np.pi*x) - (np.cos(np.pi*x))**2 )  
    return s


uD = u(0)
g = - up(L)



def VF(f,uD,g,L,N):

# pas du maillage constant 
    h = L/N

    X = np.linspace(h/2,L-h/2,N)

# matrice A tridiagonale 
    dvec = np.ones(N)*2/h
    dvec[0] = 3/h
    dvec[N-1] = 1/h
    hdvec = - np.ones(N-1)*1/h 

    Ah = np.diag(dvec) + np.diag(hdvec,1) + np.diag(hdvec,-1) 

    Sh = h*f(X)
    Sh[0] = Sh[0] + uD*2/h 
    Sh[N-1] = Sh[N-1] - g  

# solution Uh 
    Uh = np.linalg.solve(Ah,Sh)
    
    return X,Uh


# nombre de mailles
N= 40

X,Uh = VF(f,uD,g,L,N)

#tracer les courbes
plt.figure(1)
plt.clf()
Xfine = np.linspace(0,L,200)
plt.plot(Xfine,u(Xfine), label="solution_exacte",marker='*') 
plt.plot(X,Uh, label="solution_approchee")
plt.legend(loc="upper left")
plt.ylim(0.8, 3.5)
plt.show()


# etude de la convergence du schema fct de h = L/N 
Nmesh = 8
sizeh = np.zeros(Nmesh)
erreurl2 = np.zeros(Nmesh)
erreurh10 = np.zeros(Nmesh)

for imesh in range(Nmesh):

    N = 10*2**imesh

    X,Uh = VF(f,uD,g,L,N)

    h = L/N 
    
    sizeh[imesh] = h

# erreur l2 discrete 
    erl2 = 0
    for i in range(N):
        erl2 = erl2 + h*(Uh[i]-u(X[i]))**2 # quadrature du pt milieu
    
    erl2 = np.sqrt(erl2)
    print("erreur l2 ",erl2)    
    erreurl2[imesh] = erl2

# erreur h10 discrete 
    erh10 = 0
    for i in range(N-1): 
        eri = u(X[i]) - Uh[i]
        erip1 = u(X[i+1]) - Uh[i+1] 
        erh10 = erh10 + (eri-erip1)**2/h 

    er0 = u(X[0]) - Uh[0]
    erh10 = erh10 + er0**2*2/h     
    erh10 = np.sqrt(erh10)
    print("erreur h10 ",erh10)   
    erreurh10[imesh] = erh10
    
    
plt.figure(2)
plt.ylabel(" erreur l2 et h10 ")
plt.xlabel(" pas du maillage h ")
plt.loglog(sizeh,erreurl2,"-xb", label="erreur l2")    
plt.loglog(sizeh,erreurh10,"-xr", label="erreur h10")    
plt.legend(loc="upper left")

droite=np.polyfit(np.log(sizeh),np.log(erreurl2),1)
print("ordre de convergence l2",droite[0])

droite=np.polyfit(np.log(sizeh),np.log(erreurh10),1)
print("ordre de convergence h10",droite[0])








