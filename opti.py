# Circuit WeeMet-all
# Optimisation du bloc inductance et de l'oscillateur
# Groupe 11.61

from numpy import *
from matplotlib import pyplot as plt

Vcc = 6 # [V] alimentation
C_G = 470e-12 # [F] capacité du bloc oscillateur
R_L = 330 # [ohm] résistance du bloc inductance
    
L_n = 0.457e-3 # [H] inductance sans métal
tau_n = L_n/R_L # [s]
L_m = 0.427e-3 # [H] inductance avec métal
tau_m = L_m/R_L # [s] 

R = linspace(10, 35000, 100000) # [ohm] résistances du bloc oscillateur
T = 2*log(2)*C_G*R # [s] périodes du bloc oscillateur
f = 1/T # [Hz] fréquences du bloc oscillateur

D = lambda T : Vcc*((1/(1+exp(-T/(2*tau_n)))) - (1/(1+exp(-T/(2*tau_m)))))
Diff = abs(D(T)) # Différence entre le max de VL avec et sans métal

print("Différence maximale : {} [V]".format(max(Diff)))
print("Résistance optimale : {} [ohm]".format(R[argmax(Diff)]))
print("Période optimale : {} [s]".format(T[argmax(Diff)]))
print("Fréquence optimale : {} [Hz]".format(f[argmax(Diff)]))

print(T[argmax(Diff)]/(2*log(2)))



plt.figure("Optimisation oscillateur")

#plt.plot(R, Diff)
#plt.xlabel("Résistance $R_G$ [$\Omega $]")

plt.plot(T*1e6, Diff)
plt.xlabel("Période [$\mu$s]")

#plt.plot(f, Diff)
#plt.xlabel("Fréquence [Hz]")

plt.ylabel("$\Delta V$ [V]")
plt.show()