# Simulation circuit WeeMet-all
# Groupe 11.61

from numpy import *
import matplotlib.pyplot as plt

###############################
#### Paramètres simulation ####
###############################


t_begin = 0 # 0.4e-3 # [s]
t_end = 2e-5 # [s]
n = 10000 # nombre de points dans t
t = linspace(t_begin, t_end, n)


#######################
##### Oscillateur #####
#######################


def simulation_oscillateur(R_G, C_G):
    
    VS_haut = alimentation # [V]
    VS_bas = 0 # [V]
    VA_haut = (2/3)*VS_haut # [V]
    VA_bas = (1/3)*VS_haut # [V]

    tau = R_G*C_G # [s]

    VS = zeros(n) ; VS[0] = VS_bas
    VA = zeros(n) ; VA[0] = VA_bas
    VG = zeros(n) ; VG[0] = VA_haut

    def charge_oscillateur(t, t_0):
        return VS_haut*(1-(2/3)*exp(-(t-t_0)/tau))

    def décharge_oscillateur(t, t_0):
        return (2/3)*VS_haut*exp(-(t-t_0)/tau)

    t_0 = 0
    chargement = False

    for i in range(n-1):
        
        if VG[i] > VA[i]:       # décharge
            
            if chargement:
                t_0 = t[i]
                chargement = False
                
            VG[i+1] = décharge_oscillateur(t[i+1], t_0)
            VA[i+1] = VA_bas
            VS[i+1] = VS_bas
            
        else:                   # charge
            
            if not chargement:
                t_0 = t[i]
                chargement = True
            
            VG[i+1] = charge_oscillateur(t[i+1], t_0)
            VA[i+1] = VA_haut
            VS[i+1] = VS_haut
            
    return VG, VA, VS

def graphe_oscillateur():
    
    plt.figure("simulation oscillateur")
    #plt.title("Tensions du bloc oscillateur en fonction du temps")
    plt.plot(t, VS, label="$V_S$", color="r")
    plt.plot(t, VA, label="$V_A$", color="y")
    plt.plot(t, VG, label="$V_G$", color="b")
    plt.xlabel("Temps [s]")
    plt.ylabel("Tension [V]")
    plt.legend()
    plt.show()


######################
##### Inductance #####
######################


def simulation_inductance_double(VS, R_L, L_n, L_m): # prend comme paramètre le signal carré du bloc oscillateur
    
    VS_haut = max(VS) # [V]
    
    tau_n = L_n/R_L # [s]
    tau_m = L_m/R_L # [s]
    
    Vn = zeros(n)
    Vm = zeros(n)
    
    def charge_inductance(t, t_0, tau, V_0):
        return (V_0-VS_haut)*exp(-(t-t_0)/tau) + VS_haut
    
    def décharge_inductance(t, t_0, tau, V_0):
        return V_0*exp(-(t-t_0)/tau)
    
    t_0 = 0
    Vn_0 = alimentation/(1+exp(-T/(2*tau_n)))
    Vm_0 = alimentation/(1+exp(-T/(2*tau_m)))
    chargement = False
    
    for i in range(n):
        
        if VS[i] == VS_haut: # charge
            
            if not chargement:
                t_0 = t[i]
                Vn_0 = Vn[i-1]
                Vm_0 = Vm[i-1]
                chargement = True
                
            Vn[i] = charge_inductance(t[i], t_0, tau_n, Vn_0)
            Vm[i] = charge_inductance(t[i], t_0, tau_m, Vm_0)
        
        else:               # décharge
            
            if chargement:
                t_0 = t[i]
                Vn_0 = Vn[i-1]
                Vn_max = Vn_0 # garde en mémoire le max auquel Vn est arrivé
                Vm_0 = Vm[i-1]
                Vm_max = Vm_0 # garde en mémoire le max auquel Vm est arrivé
                chargement = False
            
            Vn[i] = décharge_inductance(t[i], t_0, tau_n, Vn_0)
            Vm[i] = décharge_inductance(t[i], t_0, tau_m, Vm_0)
            
    return Vn, Vm


def simulation_inductance_simple(VS, R_L, L_n, L_m): # prend comme paramètre le signal carré du bloc oscillateur
    
    VS_haut = max(VS) # [V]
    
    L = linspace(L_n, L_n + 2*(L_m-L_n), n)
    L[L>L_m] = L_m
    
    VL = zeros(n)
    
    def charge_inductance(t, t_0, tau, V_0):
        return (V_0-VS_haut)*exp(-(t-t_0)/tau) + VS_haut
    
    def décharge_inductance(t, t_0, tau, V_0):
        return V_0*exp(-(t-t_0)/tau)
    
    t_0 = 0
    VL_0 = alimentation/(1+exp(-T*R_L/(2*L_n)))
    chargement = False
    
    for i in range(n):
        
        if VS[i] == VS_haut: # charge
            
            if not chargement:
                t_0 = t[i]
                VL_0 = VL[i-1]
                chargement = True
            
            tau = L[i]/R_L
            VL[i] = charge_inductance(t[i], t_0, tau, VL_0)
        
        else:               # décharge
            
            if chargement:
                t_0 = t[i]
                VL_0 = VL[i-1]
                VL_max = VL_0 # garde en mémoire le max auquel VL est arrivé
                chargement = False
                
            tau = L[i]/R_L
            VL[i] = décharge_inductance(t[i], t_0, tau, VL_0)
            
    return VL

def graphe_inductance_double():
    
    plt.figure("Bloc détecteur de métal")
    plt.title("Tensions du bloc détecteur de métal en fonction du temps")
    plt.plot(t, VS, label="$V_S$", color="r")
    plt.plot(t, Vn, label="$V_{L,n}$", color="y")
    plt.plot(t, Vm, label="$V_{L,m}$", color="b")
    plt.xlabel("Temps [s]")
    plt.ylabel("Tension [V]")
    plt.legend()
    plt.show()

def graphe_inductance_simple():
    
    plt.figure("Bloc détecteur de métal")
    #plt.title("Tensions du bloc détecteur de métal en fonction du temps")
    plt.plot(t, VS, label="$V_S$", color="b")
    plt.plot(t, VL, label="$V_L$", color="r")
    plt.xlabel("Temps [s]")
    plt.ylabel("Tension [V]")
    plt.legend()
    plt.show()
    

##############################
##### Détecteur de crête #####
##############################


def simulation_crête(VL, C_det, R_det): # prend en paramètre le signal d'entrée du bloc
    
    tau = R_det*C_det # [s]
    
    VC = zeros(len(VL)); VC[0] = VL[0]
    
    def décharge_crête(t, t_0, V_0):
        return V_0*exp(-(t-t_0)/tau)
    
    t_0 = 0
    VC_0 = VC[0]
    chargement = False
    
    for i in range(n-1):
        
        if VC[i] > VL[i+1]: # si VC est plus grand que VL, VC se décharge
            
            if chargement:
                t_0 = t[i]
                VC_0 = VC[i]
                chargement = False
                
            VC[i+1] = décharge_crête(t[i+1], t_0, VC_0)
            
        else : # VC se charge avec VL
            
            if not chargement:
                chargement = True
            
            VC[i+1] = VL[i+1]
    
    return VC, VL, min(VC)

def graphe_crête():
    
    plt.figure("Bloc détecteur de crête")
    #plt.title("Tensions du bloc détecteur de crête en fonction du temps")
    plt.plot(t, VL, label="$V_L$", color="y")
    plt.plot(t, VC, label="$V_C$", color="b")
    plt.xlabel("Temps [s]")
    plt.ylabel("Tension [V]")
    #plt.ylim(4, 4.25)
    #plt.ylim(-0.1, 5.1)
    #plt.xlim(0, 0.001)
    plt.legend()
    plt.show()
    
def graphe_multiple_crête(VL, VC1, VC2, VC3):
    
    plt.figure()
    plt.plot(t[1:]-t_begin, VL[1:], label="$V_L$")
    plt.plot(t[1:]-t_begin, VC1[1:], label="$V_C$ 330$\Omega$")
    plt.plot(t[1:]-t_begin, VC2[1:], label="$V_C$ 990$\Omega$")
    plt.plot(t[1:]-t_begin, VC3[1:], label="$V_C$ 10k$\Omega$")
    plt.xlabel("Temps [s]")
    plt.ylabel("Tension [V]")
    plt.ylim(-0.5, 2.5)
    plt.legend(loc="lower left")
    plt.show()
    
    
def signal_triangle(f, V_max):
    
    T = 1/f # [s]
    V = zeros(n); V[0] = 0
    
    t_0 = 0
    chargement = True
    k = 0 # nombre de périodes effectuées
    
    for i in range(n-1):
        
        if t[i] - k*T < T/2: # charge
            
            if not chargement:
                t_0 = t[i]
                chargement = True
            
            V[i+1] = (2*V_max/T)*(t[i+1]-t_0)
        
        else: # décharge
            
            if chargement:
                t_0 = t[i]
                chargement = False
                
            V[i+1] = -(2*V_max/T)*(t[i+1]-t_0-(T/2))
            
            if t[i] - k*T > T:
                k+=1
            
    return V


######################
#### Soustracteur ####
######################


def simulation_soustracteur(VC, R_L, L_n, L_m):
    
    tau_n = L_n/R_L # [s]
    tau_m = L_m/R_L # [s]
    Vn_max = alimentation/(1+exp(-T/(2*tau_n)))
    Vm_max = alimentation/(1+exp(-T/(2*tau_m)))
    
    rapport = 3/(Vm_max-Vn_max)
    Vsous = (3*Vm_max)/(3+Vm_max-Vn_max)
    #print(Vsous, rapport)
    
    VF = (1+rapport)*Vsous - rapport*VC
    VF[VF>alimentation] = alimentation # toutes les valeurs de VF plus grandes que 5 sont ramenées à 5
    VF[VF<0] = 0 # toutes les valeurs de VF plus petites que 0 sont ramenées à 0
        
    return VF, Vsous*ones(n), rapport

def graphe_soustracteur_double():
    
    fig, (ax1, ax2) = plt.subplots(2)
    fig.suptitle("Tensions du bloc soustracteur de métal en fonction du temps")
    
    #ax1.plot(t, (1+rapport)*Vsous, label="$aV_{sous}$"+", a={:.2f}".format(1+rapport), color="b")
    #ax1.plot(t, rapport*VC, label="$bV_C$"+", b={:.2f}".format(rapport), color="r")
    #ax1.plot(t, Vsous, label="$V_{sous}$", color="b")
    ax1.plot(t, VC, label="$V_C$", color="r")
    ax1.set(ylabel="Tension [V]")
    ax1.legend()
    
    #ax2.plot(t, Vsous, label="$V_{sous}$", color="b")
    #ax2.plot(t, VC, label="$V_C$", color="r")
    ax2.plot(t, VF, label="$V_F = aV_{sous} - bV_C$", color="y")
    ax2.set(xlabel="Temps [s]", ylabel="Tension [V]")
    ax2.legend()
    
    plt.show()
    
def graphe_soustracteur_simple():
    
    plt.figure()
    plt.plot(t, VC, label="$V_C$", color="r")
    plt.plot(t, VF, label="$V_F$", color="b")
    plt.plot(t, (1+rapport)*Vsous, "--", label="$a V_{sous}$", color="y")
    plt.xlabel("Temps [s]")
    plt.ylabel("Tension [V]")
    plt.legend()
    plt.show()


###################
#### Signaleur ####
###################


def simulation_signaleur(VF, Vref):
    
    Vout = zeros(n)
    Vout[VF<=Vref] = alimentation
    
    return Vout, Vref*ones(n)

def graphe_signaleur():
    
    plt.figure("Bloc signaleur")
    #plt.title("Tensions du bloc signaleur en fonction du temps")
    plt.plot(t, VF, label="$V_F$", color="r")
    plt.plot(t, Vout, label="$V_{out}$", color="b")
    plt.plot(t, Vref, "--", label="$V_{ref}$", color="y")
    plt.xlabel("Temps [s]")
    plt.ylabel("Tension [V]")
    plt.legend()
    plt.show()


########################
#### Tous les blocs ####
########################


def graphe_blocs(Vcc, VS, VL, VC, VF, Vout):
    
    fig, axs = plt.subplots(3, 2)
    fig.suptitle("Tensions dans le circuit WeeMet-all")

    axs[0, 0].plot(t, Vcc, label="$V_{cc}$")
    axs[0, 0].legend()
    axs[0, 0].set(ylabel="Tension [V]")
    #axs[0, 0].set_title("Entrée bloc oscillateur")

    axs[0, 1].plot(t, VS, "tab:purple", label="$V_S$")
    axs[0, 1].legend()
    #axs[0, 1].set_title("Sortie bloc oscillateur")

    axs[1, 0].plot(t, VL, "tab:pink", label="$V_L$")
    axs[1, 0].legend()
    axs[1, 0].set(ylabel="Tension [V]")
    #axs[1, 0].set_title("Sortie bloc détecteur de métal")

    axs[1, 1].plot(t, VC, "tab:red", label="$V_C$")
    axs[1, 1].legend()
    axs[1, 1].set_ylim(0.7, 5.4) #4.4
    #axs[1, 1].set_title("Sortie détecteur de crête")

    axs[2, 0].plot(t, VF, "tab:orange", label="$V_F$")
    axs[2, 0].legend()
    axs[2, 0].set_ylim(-0.1, 4)
    axs[2, 0].set(xlabel="Temps [s]", ylabel="Tension [V]")
    #axs[2, 0].set_title("Sortie bloc soustracteur")

    axs[2, 1].plot(t, Vout, "tab:green", label="$V_{out}$")
    axs[2, 1].legend()
    axs[2, 1].set(xlabel="Temps [s]")
    #axs[2, 1].set_title("Sortie bloc signaleur")

    plt.show()

def graphe_schéma():
    
    plt.figure()
    plt.title("Avec métal")
    #plt.plot(t[2500:7500], Vcc[2500:7500], label="$V_{cc}$", color = "tab:blue")
    #plt.plot(t[2500:7500], VS[2500:7500], label="$V_S$", color = "tab:purple")
    #plt.plot(t[2500:7500], VL[2500:7500], label="$V_L$", color="tab:pink")
    #plt.plot(t[2500:7500], VC[2500:7500], label="$V_C$", color="tab:red")
    #plt.plot(t[2500:7500], VF[2500:7500], label="$V_F$", color="tab:orange")
    plt.plot(t[2500:7500], Vout[2500:7500], label="$V_{out}$", color="tab:green")
    plt.xlabel("Temps")
    plt.ylabel("Tension")
    plt.ylim(-0.15, 6)
    plt.xticks([])
    plt.yticks([])
    #plt.yticks([2.5], ["$V_{ref}$"])
    #plt.yticks([2.5, 5], ["$V_{ref}$", "$V_{cc}$"])
    plt.legend()
    plt.show()


############################
#### Paramètres circuit ####
############################


alimentation = 6 # [V]
Vcc = alimentation*ones(n)

R_G = 6340 # [ohm]
C_G = 0.47e-9 # [F]
T = 2*log(2)*R_G*C_G # [s]

R_L = 330 # [ohm]
L_n = 0.457e-3 # [H] inductance sans métal
L_m = 0.427e-3 # [H] inductance avec métal

C_det = 1e-6 # [F]
R_det = 100000 # [ohm]

Vref = 2.5 # [V]


######################
#### Instructions ####
######################


"""Oscillateur"""
VG, VA, VS = simulation_oscillateur(R_G, C_G)
#graphe_oscillateur()

"""Inductance"""
VL, Vm = simulation_inductance_double(VS, R_L, L_n, L_m)
#print(max(Vn), max(Vm), max(Vn)-max(Vm))
#graphe_inductance_double()

#VL = simulation_inductance_simple(VS, R_L, L_n, L_m)
#print(max(VL))
graphe_inductance_simple()

"""Détecteur de crête"""
V_tri = signal_triangle(1000, 2) # signal triangulaire 1 [kHz] entre 0 et 2 [V]
VC1, VL, VC_min = simulation_crête(V_tri, C_det, 330)
VC2, VL, VC_min = simulation_crête(V_tri, C_det, 990)
VC3, VL, VC_min = simulation_crête(V_tri, C_det, 10000)
#graphe_multiple_crête(VL, VC1, VC2, VC3)

VC, VL, VC_min = simulation_crête(Vm, C_det, R_det)
#print(max(VC))
#print("Tension minimale de VC en régime = {} [V]".format(VC_min))
#print("VC_max = {} [V]".format(max(VC)))
#graphe_crête()

"""Soustracteur"""
#VC = signal_triangle(1000, 1)
VF, Vsous, rapport = simulation_soustracteur(VC, R_L, L_n, L_m)
#graphe_soustracteur_double()
#graphe_soustracteur_simple()

"""Signaleur"""
#VF = signal_triangle(1000, 2)
Vout, Vref = simulation_signaleur(VF, Vref)
#graphe_signaleur()

"""Tous les blocs"""
#graphe_blocs(Vcc, VS, VL, VC, VF, Vout)

#graphe_schéma()