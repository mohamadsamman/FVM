clear

// Equation elliptique 1D sur le domaine (0,L)
//
// -u''(x) = f(x) sur (0,L)
//  u(0) = uD
//  -u'(L) = g 
//

//jeu de donnees

L=1; // longueur du domaine (0,L) 

// solution exacte u 
function  s = u(x)
s = exp(sin(%pi.*x))  
endfunction

// derivee de la solution exacte u' 
function  s = up(x)
s = %pi*cos(%pi.*x).*exp(sin(%pi.*x))  
endfunction

// second membre f = -u'' 
function  s = f(x)
s = %pi^2.*exp(sin(%pi.*x)).*(sin(%pi.*x) - (cos(%pi.*x))^2 )  
endfunction

uD = u(0)

g = -up(L)



//N = 64

for imesh=1:8
    N = 4*2**(imesh-1)
//////////////// discretization FV ////////

//pas du maillage uniforme 
h=L/N;

//coordonnees des centres de maille 

X = linspace(h/2,L-h/2,N)
Xe = linspace(0,L,1000)

dvec = ones(N,1)*2/h
dvec(1) = 3/h 
dvec(N) = 1/h 

hdvec = - ones(N-1,1)*1/h

Ah = diag(dvec) + diag(hdvec,1) + diag(hdvec,-1)

Sh = h.*f(X)'
Sh(1) = Sh(1) + uD*2/h
Sh(N) = Sh(N) - g

Uh = Ah\Sh

scf(1);
xset("font size",5)
xlabel("X","fontsize",5) ;
ylabel("U","fontsize",5) ;
plot (Xe,u(Xe),'-r')  
plot (X,Uh,'-b')  

erl2 = 0
for i=1:N
    erl2 = erl2 + h*( Uh(i)-u(X(i)) )**2
end
erl2 = sqrt(erl2)
disp('erl2',erl2)

erh10 = 0
for i=1:N-1
    erh10 = erh10 + ( u(X(i))-Uh(i) - (u(X(i+1)) - Uh(i+1)) )**2/h
end
erh10 = erh10 + (u(X(1)) - Uh(1))**2*2/h
erh10 = sqrt(erh10)
disp('erh10',erh10)

sizeh(imesh) = h
erreurl2(imesh) = erl2
erreurh10(imesh) = erh10

end


scf(32);
xset("font size",5)
xlabel("log2 h","fontsize",5) ;
ylabel("log2 erreur","fontsize",5) ;
plot(log(sizeh)/log(2),log(erreurl2)/log(2),'-xr')  
plot(log(sizeh)/log(2),log(erreurh10)/log(2),'-xb')   




