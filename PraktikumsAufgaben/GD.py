import numpy as np
import matplotlib.pyplot as plt
import uncertainties as unc
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const
import sympy
import os

#aufgabe1
m1=np.array([2,3,4,5,6])
x=np.array([1.6,2.7,3.2,3.5,4.0])

m,n = np.polyfit(m1, x, deg=1, cov=True)
err = np.sqrt(np.diag(n))
print('Steigung und startpunkt der Ausgleichsgerade')
print (m)
print('Fehler der ausgleichsgerade')
print(err)

k=const.g/m[0]

print('k=')
print(k)

plt.rcParams['figure.figsize'] = (10, 5)
plt.rcParams['font.size'] = 12
plt.rcParams['lines.linewidth'] = 0.75

plt.plot(m1,x,'+r', label="Werte von T")
plt.plot(m1, m[0]*m1+m[1], '-k',label='Regression von T')
plt.tight_layout()
plt.legend()
plt.grid()
plt.show()

#aufgabe 2
#   a)

g=np.array([60,80,100,110,120,125])
b=np.array([285,142,117,85,86,82])
n=6

def poly(g,b):
    return 1/((1/g)+(1/b))
Brennweitef = poly(g,b)
print('Brennweite f')
print (Brennweitef)

def mittelf(Brennweitef,n):
    return np.sum(Brennweitef)/n
def stabwf(Brennweitef,n):
    return np.sqrt(1/((n-1))*np.sum((mittelf(Brennweitef,n)-(Brennweitef))**2))
def fehlermittel(Brennweitef,n):
    return np.sqrt(1/(n*(n-1))*np.sum((mittelf(Brennweitef,n)-Brennweitef)**2))

print('Mittelwert von f')
print(mittelf(Brennweitef,n))
print('StandardAbweichung von f')
print(stabwf(Brennweitef,n))
print('Fehler des Mittelwerte')
print(fehlermittel(Brennweitef,n))

#   b)
G=1/g
B=1/b

params , _ = np.polyfit(B,G, deg=1, cov=True)
err = np.sqrt( np.diag( _ ) )

m=ufloat(params[0],err[0])
a=ufloat(params[1],err[1])

print('steigung der ausgleichsgeraden')
print(m)
print('anfang der Ausgleichsgeraden')
print(a)

plt.plot(B,G,'.',label=r'Messwerte')
plt.plot(B,(m.n*B+a.n),label=r'lineare Regression')
plt.legend()
plt.legend(loc='best')
plt.xlabel(r'$\frac{1}{b}$ / $\frac{1}{mm}$' )
plt.ylabel(r'$\frac{1}{g}$ / $\frac{1}{mm}$' )
plt.show()

#Aufgabe 3

d=np.array([0.1,0.2,0.3,0.4,0.5,1.0,1.2,1.5,2.0,3.0,4.0,5.0])
N=np.array([7565,6907,6214,5531,4942,2652,2166,1466,970,333,127,48])

def exp(x,n,m):
    return n*np.exp(- x * m)
params, covariance_matrix = curve_fit(exp,d,N)

uncertainties = np.sqrt(np.diag(covariance_matrix))

n=ufloat(params[0],uncertainties[0])
m=ufloat(params[1],uncertainties[1])
print('m')
print(m)
print('n')
print(n)
fig, (ax1,ax2)=plt.subplots(2,1)

ax1.set_title(r"Nicht Logarithmisch")
ax1.plot(d,np.sqrt(N),'.',label=r'Messunsicherheit')
ax1.plot(d,exp(d,n.n,m.n),label=r'Fit')
ax1.errorbar(d, N, yerr=np.sqrt(N)*5 , fmt='.' ,label=r'Messwerte')
ax1.legend()
ax1.legend(loc='best')
ax1.set_xlabel('$d/cm$')
ax1.set_ylabel('$N/60s^-1$')

ax2.set_title(r"$Logarithmisch$")
ax2.errorbar(d, N, yerr=np.sqrt(N)*5, fmt='.',label=r'Messwerte')
ax2.plot(d,exp(d,n.n,m.n),label=r'Fit')
ax2.plot(d,np.sqrt(N),'.',label=r'Messunsicherheit')
ax2.set_yscale('log')
ax2.legend()
ax2.legend(loc='best')
ax2.set_xlabel('$d/cm$')
ax2.set_ylabel('$N/60s^-1$')

fig.tight_layout()

plt.show()