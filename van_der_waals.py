import numpy as np
import matplotlib.pyplot as plt
#
# Computes the isotherms of the Van der Waals equation
# given input constants [a,b]
#
# (P + a(n/V)^2) ( V/n - b) = RT
#

def van_der_waals_P(V,T,a,b,n):

    nonV = n/V

    P = (R*T - a*nonV*nonV)/(1.0/nonV - b)
    return P

def van_der_waals_T(P,V,a,b,n):

    nonV = n/V

    T = (P + a*nonV*nonV)*(1.0/nonV - b)/R
    return T

def van_der_waals_V(P,T,a,b,n):

    coeff = [P-b, a, 0.0, -R*T-a*b]

    return np.roots(coeff)



blank = [1.0]


# Ideal Gas
idealgaslabel = 'Ideal Gas'
aideal=0
bideal=0


# Constants for water
waterlabel = r'Water (H$_2$O)'
awater = 553.6
bwater = 0.03049

# Number of moles
n=100.0
# Ideal gas constant (SI)
R = 8.314

npoints = 100

Pcrit = awater/(27.0*bwater*bwater)
Tcrit = 8.0*awater/(27*R*bwater)

Vcrit = 3.0*bwater*n

print ("Water Critical Pressure, Temperature, Volume: ",Pcrit/1.0e3, Tcrit,Vcrit)

# Pressure (kiloPascals)
P = np.linspace(0.1,5.0e1, npoints)

# Volume (Litres)
V = np.linspace(0.01,50, npoints)

Tideal = np.zeros((npoints,npoints))
Twater = np.zeros((npoints,npoints))


for i in range(npoints):
    for j in range(npoints):

        Tideal[i,j] = van_der_waals_T(1000.0*P[i],V[j],aideal,bideal,n)
        Twater[i,j] = van_der_waals_T(1000.0*P[i],V[j],awater,bwater,n)

fig1 = plt.figure(figsize=(8,6))
fig1.suptitle("Isotherms for "+str(n)+" moles of "+idealgaslabel, fontsize=22)
ax1 = fig1.add_subplot(111)
ax1.set_xlabel('Volume (L)',fontsize=22)
ax1.set_ylabel('Pressure (kPa)',fontsize=22)
CS1 = ax1.contour(V,P,Tideal, colors='black',levels=list(range(0,1000,50)))
ax1.plot(blank,blank, color='black', label=r'T(K)') # Dummy plot to add T to legend
plt.clabel(CS1, inline=1, fmt='%i')

ax1.legend(loc='upper right', fontsize=22)

fig1.savefig('PV_idealgas.png')

fig2 = plt.figure(figsize=(8,6))
fig2.suptitle("Isotherms for "+str(n)+" moles of "+waterlabel, fontsize=22)
ax2 = fig2.add_subplot(111)
ax2.set_xlabel('Volume (L)',fontsize=22)
ax2.set_ylabel('Pressure (kPa)',fontsize=22)
CS2 = ax2.contour(V,P,Twater, colors='black',levels=list(range(0,1000,50)))
CS3 = ax2.contour(V,P,Twater,colors='red', linewidth=2, levels=[Tcrit])
ax2.plot(blank,blank, color='black', label=r'T(K)') # Dummy plot to add T to legend
ax2.plot(blank,blank, color='red',linewidth=2, label=r'$T_{\rm crit}$') # Dummy plot to add Tcrit to legend
ax2.axhline(Pcrit/1.0e3, linestyle='dashed', color = 'green', label = r'$P_{\rm crit}$')
ax2.axvline(Vcrit, linestyle='dashed', color='blue', label = r'$V_{\rm crit}$')

ax2.legend(loc='upper right',fontsize=22)
plt.clabel(CS2, inline=1, fmt='%i')

fig2.savefig('PV_water.png')

plt.show()








