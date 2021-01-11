import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import uncertainties.unumpy as unp
import uncertainties as unc

PressureP1 = np.array([30,50,108,153,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000])#,1040,1050])
TemperatureT1 = np.array([22+273.15,26+273.15,40+273.15,51+273.15,58+273.15,64+273.15,69+273.15,73+273.15,76+273.15,79+273.15,83+273.15,86+273.15,88+273.15,91+273.15,94+273.15,97+273.15,98+273.15,100+273.15,102+273.15,104+273.15,108+273.15])#,108+273.15,109+273.15])

#umrechnung des Drucks von bar zu Pascal
PressureP1 = PressureP1*10**2


print('Druck in Pascal')
print (PressureP1)
print('Temperatur in Kelvin')
print (TemperatureT1)

#Logarythmus vom Druck ausrechnen
def poly(PressureP1):
    return np.log(PressureP1/10**5)
Logarytmus = poly(PressureP1)
print('logarythmus von dem Druck')
print (Logarytmus)

#Kehrwert von der Temperatur ausrechnen
def poly2(TemperatureT1):
    return 1/TemperatureT1
Kehrewert = poly2(TemperatureT1)
print('Kehrwert der Temperatur')
print (Kehrewert)

#ausgleichsgerade ausrechnen
m,n = np.polyfit(Kehrewert, Logarytmus, deg=1, cov=True)
err = np.sqrt(np.diag(n))
print('Steigung und startpunkt der Ausgleichsgerade')
print (m)
print('Fehler der ausgleichsgerade')
print(err)

#ausrechnen der molaren Verdampfungswärme
um = unp.uarray(m,err)
def poly3(um):
    return (-um[0])*const.gas_constant
uL = poly3(um)
print('molare Verdampfungswärme mit unsicherheit')
print(uL)

#ausrechnen der Äußere Verdampfungswärme
La = 373*const.gas_constant
uLa=unc.ufloat(La,0)
print('Äußere verdampungswärme mit unsicherheit')
print(uLa)

#ausrechnen der gesammten Verdampfungswärme
Li = uL-uLa
print('verdampfungswärme mit unsicherheit')
print(Li)

#umrechnung zu elekrtonenvolt
Li = Li/const.N_A
Li = Li/const.e
print('Li in elektronenvolt')
print(Li)

plt.rcParams['figure.figsize'] = (10, 5)
plt.rcParams['font.size'] = 12
plt.rcParams['lines.linewidth'] = 0.75

plt.plot(Kehrewert,Logarytmus,'+r', label="Werte von T")
plt.plot(Kehrewert, m[0]*Kehrewert+m[1], '-k',label='Regression von T')
plt.ylabel('[ln(p/p0)] in Pas')
plt.xlabel('[1/T] in K')
plt.tight_layout()
plt.legend()
plt.grid()

plt.show()

PressureP2 = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
TemperatureT2 = np.array([118+273.15,129.5+273.15,141+273.15,150+273.15,157+273.15,163.5+273.15,168+273.15,172.5+273.15,177+273.15,181.5+273.15,185.5+273.15,189+273.15,192+273.15,195+273.15,198.5+273.15])
LitTemp = np.array([0+273.15,25+273.15,40+273.15,60+273.15,80+273.15,100+273.15,120+273.15,140+273.15,160+273.15,180+273.15,200+273.15,220+273.15,240+273.15,260+273.15,280+273.15,300+273.15,320+273.15,340+273.15,360+273.15,374+273.15])
LitVerd = np.array([45054,43990,43350,42482,41585,40657,39684,38643,37518,36304,34962,33968,31809,29930,27795,25300,22297,18502,12966,2067])
#umrechnen des Drucks von Bar zu Pascal
PressureP2 = PressureP2*1e5
x = TemperatureT2

print('Druck in Pascal')
print(PressureP2)
print('Temperatur in Kelvin')
print(TemperatureT2)

#ausgleichsgerade und Fehler berechnen
g,b = np.polyfit(TemperatureT2, PressureP2, deg=3, cov=True)
err2 = np.sqrt(np.diag(b))
print('Werte der Ausgleichsfunktion')
print(g)
print('Fehler der Ausgleichsfunktion')
print(err2)

a = g[0]
b = g[1]
c = g[2]
d = g[3]

R=const.gas_constant
C=0.9

#ausrechnen von L(T) mit positiver und negativer Wurzel
def L1(x,a,b,c,d):
    return (((R*x/2)+np.sqrt((R*x/2)**2-0.9*(a*x**3 + b*x**2+c*x+d)))*((3*a*x**3+2*b*x**2+c*x)/(a*x**3+b*x**2+c*x+d)))
L1(x,a,b,c,d)
print('Werte für L(T) bei positiver Wurzel')
print(L1(x,a,b,c,d))

def L2(x,a,b,c,d):
    return (((R*x/2)-np.sqrt((R*x/2)**2-0.9*(a*x**3 + b*x**2+c*x+d)))*((3*a*x**3+2*b*x**2+c*x)/(a*x**3+b*x**2+c*x+d)))
L2(x,a,b,c,d)
print('werte für L(T) bei negativer Wurzel')
print(L2(x,a,b,c,d))

plt.plot(TemperatureT2,PressureP2, '+r', label="Messwerte des Drucks")
plt.plot(TemperatureT2, g[0]*TemperatureT2**3+g[1]*TemperatureT2**2+g[2]*TemperatureT2+g[3], '-k', label="Ausgleichskurve")
plt.grid()
plt.xlabel("T [K]")
plt.ylabel("P [P]")
plt.tight_layout()
plt.legend()
plt.show()

plt.plot(TemperatureT2,L1(x,a,b,c,d), '+k', label="Berechnete L Werte")
plt.grid()
plt.xlabel("T [K]")
plt.ylabel("L [J/mol]")
plt.tight_layout()
plt.legend()
plt.show()

plt.plot(TemperatureT2,L2(x,a,b,c,d), '+k', label="Berechnete L Werte")
plt.grid()
plt.xlabel("T [K]")
plt.ylabel("L [J/mol]")
plt.tight_layout()
plt.legend()
plt.show()

plt.plot(LitTemp,LitVerd, '+k', label="Literaturwerte für L")
plt.grid()
plt.xlabel("T [K]")
plt.ylabel("L [J/mol]")
plt.tight_layout()
plt.legend()
plt.show()