import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import uncertainties.unumpy as unp
import uncertainties as unc

tk = np.array([15.7,15.8,15.6,15.6,15.6,15.6,15.6,15.6,15.7,15.7])
tg = np.array([132.3,131.9,132.1,131.8,131.8,131.3,131.5,131.0,131.2,131.1,130.6])
t1 = np.array([131.2,109.8,93.3,78.3,70.5,62.4,55,49.2,44.0,42.5])
t2 = np.array([131.1,109.1,92.7,80.1,70.5,62.9,55.1,49,45.2,42.8])
t3 = np.array([130.6,109.6,92.4,79.6,70.6,62.5,55,48.9,44.5,42])
dw = np.array([0.998103,0.996236,0.99379,0.99083,0.9874,0.98361,0.794,0.97484,0.96996,0.96605])
Temp = np.array([20.5,28,35.7,43.5,51.4,59.2,67.1,75,82.9,88.9])
#def poly(a,b,c,d,e,f,g,h,i,j):
#    return (a+b+c+d+e+f+g+h+i+j)/10
#print('mittelwert von tk')
#print(poly(tk[0],tk[1],tk[2],tk[3],tk[4],tk[5],tk[6],tk[7],tk[8],tk[9]))
#print('mittelwert von tg')
#print(poly(tg[0],tg[1],tg[2],tg[3],tg[4],tg[5],tg[6],tg[7],tg[8],tg[9]))

#mittelwerte und Standardabweichungen
mtk = np.mean(tk, axis=0)
stk = np.std(tk) 
umtk = unp.uarray(mtk,stk)
print('mittelwert von tk')
print(umtk)
mtg = np.mean(tg, axis=0)
stg = np.std(tg)
umtg = unp.uarray(mtg,stg)
print('mittelwert von tg')
print(umtg)

upk = unc.ufloat(2.233,0.009)
upg = unc.ufloat(2.229,0.008)
upfl = unc.ufloat(0.99821,0)
uKk = unc.ufloat(5.08,0.08)
uKk = uKk*10**(-8)
#Viskosität berechnen
v = uKk*(upk-upfl)*umtk
print('Viskosität')
print(v)
#Kg berechnen
uKg = v/(umtg*(upg-upfl))
print('K der großen Kugel')
print(uKg)

v1 = uKg*(upg-dw)*t1
print('viskosität vont t1')
print(v1)
v2 = uKg*(upg-dw)*t2
print('viskosität von t2')
print(v2)
v3 = uKg*(upg-dw)*t3
print('viskosität von t3')
print(v3)

#for a,b,c in zip(v1,v2,v3):
#    print(np.mean([a,b,c]))
#mV = np.mean([a,b,c])

mV = [np.mean([*x]) for x in zip(v1,v2,v3)]
print('Mittelwert der Viskosität')
print(mV)

def poly3(P1):
    return np.log(unp.nominal_values(P1))
print('logarythmus von mV')
print(poly3(mV))



def poly2(TemperatureT1):
    return 1/TemperatureT1
print('Kehrwert der Temperatur')
print(poly2(Temp))

q,r = np.polyfit(poly2(Temp), poly3(mV), deg=2, cov=True)

plt.rcParams['figure.figsize'] = (15, 10)
plt.rcParams['font.size'] = 12
plt.rcParams['lines.linewidth'] = 0.75

plt.plot(poly2(Temp), poly3(mV), '+k')
plt.plot(poly2(Temp), (q[0]*poly2(Temp)**2)+q[1]*poly2(Temp)+q[2], '-r')
plt.grid()
plt.xscale('log')
plt.xlabel('1/T in [1/K]')
plt.ylabel("ln(Pb/Pas)")
plt.tight_layout()
plt.legend()
plt.savefig('Plot1V106.png')
plt.show()

plt.plot(Temp, dw, '+k')
plt.plot(poly2(Temp), (q[0]*poly2(Temp)**2)+q[1]*poly2(Temp)+q[2], '-k')
plt.grid()
plt.xlabel('T in [K]')
plt.ylabel("D in [g/cm^3")
plt.tight_layout()
plt.legend()
plt.savefig('Plot2V106.png')
plt.show()
