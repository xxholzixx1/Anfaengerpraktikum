import numpy as np
import matplotlib.pyplot as plt
import uncertainties as unc
import uncertainties.unumpy as unp
import sympy

C=100               #willkürlicher wert für die summe

n0=(C/100)+1    #mit c willkürlich und n0=2 willkürlich für einfachere rechnung berechnung von n.
                #=> n zusätzliche rechnungen nötig

n1=(C/9)+1
n2=(C/0.25)+1
print('\nAufgabe 2)\n\n\tVorher beliebige Wahl von von der Summe als C=100, damit man für den Fehler 10 n0 als 2(kleinstmögliches n) setzen kann')
print('\tZusätzliche Durchführungen für Standardabweichung 3:\n\t',n1-n0,'\n\tZusätzliche Durchführungen für Standardabweichung 0.5:\n\t',n2-n0)

Rin=unc.ufloat(10,1)
Rau=unc.ufloat(15,1)
h=unc.ufloat(20,1)

V=np.pi*(Rau**2-Rin**2)*h

print('\n\n\nAufgabe 3)\n\n\tVergleichswert des Volumen:',V,'\n')

r1, r2, H = sympy.var('r1 r2 H')

err=0.1

f=H* np.pi*((r2**2)-(r1**2))

errf=sympy.sqrt((f.diff(r1)**2*err**2+f.diff(r2)**2*err**2+f.diff(H)**2*err**2))
print(f'\tFunktion: {f}\n\tFehlerfunktion: {errf}\n')
print(f'\tWert der Funtion: {f.evalf(subs={r1:10,r2:15,H:20})}')
print(f'\tFehlerwert: {errf.evalf(subs={r1:10,r2:15,H:20})}\n')
