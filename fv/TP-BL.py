import numpy as np
import matplotlib.pyplot as plt
import sys
#jeu de donnees

L=1000 # longueur du reservoir 

tf=3600*24*30*12*10  # temps final de simulation 

#rapport des viscosite de l'eau et du gaz muw/mug > 1
rmu = 10

#Vitesse en m/s  
VT = 1E-6

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


#fonction F(S) vectorielle avec decentrage 
def F(S):
    V = np.zeros(N)
    
    for i in range(1,N):
        V[i] = f(S[i-1]) - f(S[i])
    
    V[0] = f(1)-f(S[0])
    V = V/h 
    return V

X = np.linspace(h/2,L-h/2,N)

#schema Euler Explicite  
#y = x + dt*f(x)
def EulerExplicite(x,dt):
    y = x + dt*F(x)
    return y

dt0 = h/maxfp()
t = 0
S = np.zeros(N)
while (t<tf):
    dt = min(dt0,tf-t)
    t = t + dt 
    
    S = EulerExplicite(S,dt)
    
    plt.figure(1)
    plt.plot(X,S,'-b')
    
plt.figure(3)
plt.plot(X,S,'-b')
   

# derivee A = F'(S)
def DF(S):    
    A = np.zeros([N,N])
    
    for i in range(1,N):
        A[i,i-1] = fp(S[i-1])
        A[i,i]   = - fp(S[i])
        
    A[0,0] = -fp(S[0]) 
    A = A/h
    
    return A 


def EulerImplicite(x,dt):
#schema Euler Implicite 
#on resoud G(y) = 0 avec G(y) = y-x - dt*F(y)    

    eps = 1.e-6
    kmax = 100
    dsobj = 0.1
    y = x
    r = y-x - dt*F(y)
    nr = np.linalg.norm(r)
    nr0 = nr
    k=0
    while (nr/nr0>eps)and(k<kmax)and(nr>eps):
        A = np.eye(N) - dt*DF(y)
        dy = - np.linalg.solve(A,r)
        dymax = np.max(np.abs(dy))
        alpha = np.min([1,dsobj/dymax])
        #alpha = 1 # sans relaxation 
        y = y + dy*alpha
        r = y-x - dt*F(y)
        nr = np.linalg.norm(r)  
        k = k+1
    if (nr/nr0> 2*eps):
        print(" newton not converged",nr/nr0)
        sys.exit()
    return y


dt0 = tf/40
t = 0
S = np.zeros(N)
while (t<tf):
    dt = min(dt0,tf-t)
    t = t + dt 
    
    S = EulerImplicite(S,dt)
    
    plt.figure(2)
    plt.plot(X,S,'-r')
    
plt.figure(3)
plt.plot(X,S,'-r')



