import numpy as np
import matplotlib.pyplot as plt
import sys
#jeu de donnees

L=1000 # longueur du reservoir 

tf=3600*24*30*12*10  # temps final de simulation 

#rapport des viscosite de l'eau et du gaz muw/mug > 1
rmu = 10/(2*2*2*2*2*2*2)

#Vitesse en m/s  
VT = 1

#nombre d'inconnues 
N=100

#pas du maillage 
h=L/N
#fonction flux scalaire f(s) = s**2/(s**2 + (1-s)**2/rmu)  
def f(s):
    z = VT*s**2/( s**2 + (1-s)**2/rmu )
    return z 

#derivee de la fonction flux scalaire f'(s)  
def fp(s):
    z = VT*2/rmu*s*(1-s)/(s**2 + (1-s)**2/rmu )**2   
    return z

#max de la derivee de la fonction flux scalaire f'(s)  sur [0,1]
def maxfp():  
    zz = np.linspace(0,1,1000)
    return np.max(fp(zz))     

print(maxfp())

# Calcul du pas de temps Δt
max_fp = maxfp()
delta_t = h / (VT * max_fp)

print(f"Le pas de temps Δt correspondant est : {delta_t}")