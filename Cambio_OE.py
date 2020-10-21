# -*- coding: utf-8 -*-
"""
@author: Elias Obreque
"""

import matplotlib.pyplot as plt
import math
import numpy as np


NombreTle = 'GomX-4B'
#NombreTle = 'sat42788'


m_earth     = 5.9722e24
G           = 6.674e-11
mu          = 3.986044418e5
R_earth     = 6378.0000


def julianepoch(y):
    m   = 1
    J0  = 367.0*y - math.floor(7.0*(y + math.floor((m + 9.0)/12.0))/4.0) + math.floor(275.0*m/9.0) + 1721013.5
    return J0


def kepler_ite(mm, ecc):
    e_   = 0
    err  = 1.0e-10
    e0   = mm
    flag = True
    itt  = 0
    while flag:
        e_  = mm + ecc*np.sin(e0)
        if abs(e_ - e0) < err:
            flag = False
        e0  = e_
        itt = itt + 1
        e_  = e0
    return e_


def getOE(Data):
    leng = len(Data)
    i    = [ ]
    e    = [ ]
    n    = [ ]
    a    = [ ]
    raan = [ ]
    ap   = [ ]
    t    = [ ]
    f    = [ ]
    r    = [ ]
    beta = [ ]
    con  = 0
    for k in range(leng-1):
        if Data[k][0] == '2':
            i.append(float(Data[k][9:16]))
            raan.append(float(Data[k][17:25]))
            ap.append(float(Data[k][34:42]))
            e.append(float(Data[k][26:33])/10000000.0)
            n.append(float(Data[k][52:63]))
            nsat = 2.0*math.pi*float(Data[k][52:63])/86400.0
            a.append(((mu**(1.0/3.0))/(nsat**(2.0/3.0))))
            MA = float(Data[k][43:51])
            h = np.sqrt(mu*a[con]*(1 - e[con]**2))
            E = kepler_ite(MA, e[con])
            f.append(2*np.arctan(np.sqrt((1.0 + e[con])/(1.0 - e[con]))*np.tan(0.5*E)))
            r.append(((h**2)/mu)/(1 + e[con]*np.cos(f[con])))
            con = con + 1
        elif Data[k][0] == '1':
            if float(Data[k][18:20])<20:
                year = float('20'+Data[k][18:20])
            else:
                year = float('19'+Data[k][18:20])    
            J0epoch = julianepoch(year)
            epoch   = float(Data[k][20:32])
            Jepoch  = J0epoch + epoch
            beta.append(float(Data[k][53:59])*10**float(Data[k][59:61]))
            t.append(Jepoch)    
    return i, a, n, e,  raan, ap, t, r, f, beta


def Trap(fx, x):
    area = 0
    for i in range(len(x) - 1):
        area = area + 0.5*(x[i + 1] - x[i])*(fx[i + 1] + fx[i])
    return area


TLE_open    = open(NombreTle+'.txt', 'r')
TLE_read    = TLE_open.read()
Data        = TLE_read.split('\n')
TLE_open.close()

inc, a, n, e, RAAN, AP, tJD, r, f, beta = getOE(Data)

t = (np.array(tJD) - tJD[0])

NumRev = Trap(n, t)

KilTrav = 2*math.pi*(np.mean(a) + R_earth)*NumRev

print("Number of revolutions:", NumRev)
print("Kilometers traveled:", KilTrav)

plt.figure()
plt.title('Inclination')
plt.ylabel('$i$ [°]')
plt.xlabel('Time [day]')
plt.plot(t, inc)
plt.grid()

plt.figure()
plt.title('Alt')
plt.ylabel('Alt [km]')
plt.xlabel('Time [day]')
plt.plot(t, np.array(r) - R_earth)
plt.grid()

plt.figure()
plt.title('RAAN')
plt.ylabel('$RAAN$ [°]')
plt.xlabel('Time [day]')
plt.plot(t, RAAN)
plt.grid()

plt.figure()
plt.title('Arg. Per')
plt.ylabel('$\omega$ [°]')
plt.xlabel('Time [day]')
plt.plot(t, AP)
plt.grid()

plt.figure()
plt.title('Mayor Semi-axe')
plt.ylabel('$Apogee$ [km]')
plt.xlabel('Time [day]')
plt.plot(t, a)
plt.grid()

plt.figure()
plt.title('Revolutions by day')
plt.ylabel('$n$ [rev/día]')
plt.xlabel('Time [day]')
plt.plot(t, n)
plt.grid()

plt.figure()
plt.title('Eccentricity')
plt.ylabel('$e$ [-]')
plt.xlabel('Time [day]')
plt.plot(t, e)
plt.grid()

plt.figure()
plt.title('Ballistic coefficient')
plt.ylabel('beta [-]')
plt.xlabel('Time [day]')
plt.plot(t, beta)
plt.grid()
plt.show()




