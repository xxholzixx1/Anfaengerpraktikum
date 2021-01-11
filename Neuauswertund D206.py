import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import uncertainties.unumpy as unp
import uncertainties as unc

Time = np.array([0, 60,120,180,240,300,360,420,480,540,600,660,720,780,840,900,960,1020,1080,1140,1200,1260,1320,1380,1440,1500,1560,1620,1680,1740,1800,1860,1920,1980,2040,2100])
TemperatureT1 = np.array([21.7,23.0,24.3,25.3,26.4,27.5,28.8,29.7,30.9,31.9,32.9,33.9,34.8,35.7,36.7,37.6,38.4,39.2,40.0,40.7,41.4,42.2,42.9,43.6,44.3,44.9,45.5,46.1,46.7,47.3,47.8,48.4,48.9,49.4,49.9,50.3])
TemperatureT2 = np.array([21.7,21.7,21.6,21.5,20.8,20.1,19.2,18.5,17.7,16.9,16.2,15.5,14.9,14.2,13.6,13.0,12.4,11.7,11.3,10.9,10.4,9.9,9.5,9.1,8.7,8.3,8.0,7.7,7.4,7.1,6.8,5.6,4.3,3.4,3.0,2.9])
Time = Time*60
Time2 = ([480,960,1440,1920])
TemperatureT1 = TemperatureT1+273.15
TemperatureT2 = TemperatureT2+273.15
TemperatureT11 = np.array([304.5,311.55,317.45,322.05])
TemperatureT22 = np.array([209.85,285.55,281.85,277.45])
N = np.array([120,120,120,122])
pb = np.array([4.1,3.2,3.4,3.5,3.5,3.4,3.3,3.2,3.2,3,3,2.9,2.8,2.8,2.7,2.6,2.6,2.6,2.5,2.5,2.4,2.4,2.4,2.4,2.4,2.4,2.3,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2])
pa = np.array([4.0,5.0,5.5,6.0,6.0,6.0,6.5,6.5,7.0,7.0,7.0,7.5,7.5,8.0,8.0,8.0,8.5,8.5,9.0,9.0,9.0,9.0,9.5,9.5,10.0,10.0,10.0,10.0,10.5,10.5,10.75,11.0,11.0,11.0,11.0,11.0])
pa1 = np.array([7,8.5,10,11])
pb1 = np.array([3.2,2.6,2.4,2.2])
pb = pb+1
pa = pa+1
pa1 = pa1+1
pb1 = pb1+1
q = np.array([41.397,50.0719,58.739,65.095])

m,n = np.polyfit(Time, TemperatureT1, deg=2, cov=True)
err = np.sqrt(np.diag(n))
print('Steigung und startpunkt der Ausgleichsgerade')
print (m)
print('Fehler der ausgleichsgerade')
print(err)
um = unp.uarray(m,err)

o,p = np.polyfit(Time, TemperatureT2, deg=2, cov=True)
err2 = np.sqrt(np.diag(p))
print('Steigung und startpunkt der Ausgleichsgerade 2')
print (o)
print('Fehler der ausgleichsgerade 2')
print(err2)
um2 = unp.uarray(o,err2)

def poly(TemperatureT11):
    return um[0]*2*TemperatureT11+um[1]
dT1 = poly(TemperatureT11)
print('Ableitung Temperatur 1 ( mit fehler)')
print (dT1)

def poly2(TemperatureT22):
    return um[0]*2*TemperatureT22+um[1]
dT2 = poly(TemperatureT22)
print('Ableitung Temperatur 2 (mit fehler)')
print (poly(TemperatureT22))

m1c1 = 16732
mkck = 750 
def poly3(dT3, N):
    return (m1c1+mkck)*dT3*(1/N)
print('Güteziffer empirisch')
print(poly3(dT1, N))

def poly4(TemperatureT11, TemperatureT22):
    return TemperatureT11/(TemperatureT11-TemperatureT22)
print ('Güteziffer ideal')
print (poly4(TemperatureT11, TemperatureT22))

def poly5(pb):
    return np.log(pb)
log = poly5(pb)
print('logarythmus von dem Druck')
print (log)

def poly6(TemperatureT2):
    return 1/TemperatureT2
KW = poly6(TemperatureT2)
print('Kehrwert der Temperatur')
print (KW)

q,r = np.polyfit(KW, log, deg=1, cov=True)
err3 = np.sqrt(np.diag(r))
print('Steigung und startpunkt der Ausgleichsgerade')
print (q)
print('Fehler der ausgleichsgerade')
print(err3)
um3 = unp.uarray(q,err3)

L = -um3[0]*8.31451
print('Verdampfungswärme')
print(L)

L = L/120.91
print('Verdampfungswärme')
print(L)

def poly7(dT3):
    return (m1c1+mkck)*dT3*(1/L)
MD = poly7(dT1)
print('Massendurchsatz')
print (MD)

K = 1.14
T0 = 273.15
q0 = 5.51
p0 = 1

def N_mech(k,pa2,pb2,p0,T0,q0,mdt):
    return ((1/(K-1)) * ((pb2 *  (pa2/pb2)**(1/K)) -(pa2) ) * ((TemperatureT22*p0)/(q0*T0*pa2)) * mdt)*(-1)

print('Kompressorleistung in bar*l/s')
print(N_mech(K,pa1,pb1,p0,T0,q0,MD))
print ('kopressorleistung in kg*m/s')
print(N_mech(K,pa1,pb1,p0,T0,q0,MD)*100)

plt.rcParams['figure.figsize'] = (15, 10)
plt.rcParams['font.size'] = 12
plt.rcParams['lines.linewidth'] = 0.75

plt.plot(KW, log, '+k')
plt.plot(KW, (q[0]*KW)+q[1], '-k')
plt.grid()
plt.xlabel('1/T2 in [1/K]')
plt.ylabel("ln(Pb/Bar)")
plt.tight_layout()
plt.legend()
plt.savefig('Plot1D206.png')
plt.show()

plt.plot(Time, TemperatureT1, '+b')
plt.plot(Time, (m[0]*Time**2)+(m[1]*Time)+m[2], '-b')
plt.plot(Time, TemperatureT2, '+r')
plt.plot(Time,(o[0]*Time**2)+(o[1]*Time)+o[2], '-r')
plt.grid()
plt.xlabel('Zeit in [s]')
plt.ylabel("Temperatur in [K]")
plt.tight_layout()
plt.legend()
plt.savefig('Plot2D206.png')
plt.show()