# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 22-23
# Problème 5
#
# Script de test
#  Vincent Legat
#  Nathan Coppin
#  Nicolas Roisin
#
# Largement inspiré d'un code préliminaire de Nicolas Roisin :-)
# Ou les méthodes numériques pour obtenir la solution du projet P2 !
#
# -------------------------------------------------------------------------
#
# Solution d'Alexandre Mahy
# Modifié pour correspondre aux données du groupe 11.61
# dans le cadre du projet P2
#
# -------------------------------------------------------------------------


from numpy import *
from scipy.special import roots_legendre
import matplotlib.pyplot as plt


# ============================================================
# FONCTIONS A MODIFIER [begin]
#
#
# ------------------------------------------------------------------------------------
#
# Calcul des composantes radiales et verticales
# du champ magnétique généré
# par une série de boucles circulaires de courants !
#
#    X,Z : tableaux numpy des coordonnées x,z des points en [m]
#    Rsource : liste des rayons des boucles en [m]
#    Zsource : liste de hauteurs des boucles en [m]
#    Isource : liste des courants des boucles en [A]
#    data : structure contenant les paramètres matériels
#   
#  La fonction renvoie une liste [Bx,Bz] où Bx,Bz sont des tableaux/variables 
#  du même type que X,Z.   Le résultat est exprimé en [Tesla] ou [kg /s2 A]
#

def inductanceMegneticField(X,Z,Rsource,Zsource,Isource,data):
 
# 
# A COMPLETER / MODIFIER
#

    Bx  = zeros(shape(X))
    Bz  = zeros(shape(X))
  
    for j in range(len(Rsource)): # pour chaque boucle
      
        R_boucle = Rsource[j] # nombre
        Z_boucle = Zsource[j] # nombre
        I_boucle = Isource[j] # nombre
      
        Xi = R_boucle*data.Xcircle # liste
        Yi = R_boucle*data.Ycircle # liste
      
        dTheta = (2*pi)/data.nTheta # nombre
        dxi = -dTheta*Yi # liste
        dyi = dTheta*Xi # liste
      
        for i in range(data.nTheta): # pour chaque point de la boucle discrétisée
            Ri_cube = ((X-Xi[i])**2 +Yi[i]**2 +(Z-Z_boucle)**2)**(3/2) # tableau
          
            Bx += I_boucle * (dyi[i]*(Z-Z_boucle))/Ri_cube
            Bz += I_boucle * (-dxi[i]*Yi[i]-dyi[i]*(X-Xi[i]))/Ri_cube
      
    Bx *= 10**(-7)
    Bz *= 10**(-7)
  
# 
# A COMPLETER / MODIFIER
# 
 
    return [Bx,Bz]

# ------------------------------------------------------------------------------------
#
# Calcul des points et points d'intégration
# pour une intégration double de Gauss-Legendre
# afin de calculer les flux du champs magnétique
#
#    X0,Xf : intervalle radial d'intégration
#    Z0,Zf = intervalle vertical de moyenne du flux
#    nX : nombre de sous-intervalles en x
#    nZ : nombre de sous-intervalles en z
#    nGL : nombre de points de Gauss-Legendre par intervalle (en x et en z)
#   
#  La fonction renvoie une liste [X,Z,W] contenant les abscisses et les poids.
#  Ce seront des tableaux unidimensionnels de taille n*nGL*m*nGL
#  Si n ou m sont nuls, on fera une intégration unidimensionnelle en utilisant 
#  X0 ou Z0 et les tableaux auront une taille m*nGL ou n*nGL respectivement
#


def inductanceGaussLegendre(X0,Xf,Z0,Zf,n,m,nGaussLegendre):
 
# 
# A COMPLETER / MODIFIER
#
 
    xi,we = roots_legendre(nGaussLegendre)
    
    if n == 0:
        X_uni = array([X0])
        WX = 1.0
    else:
        hx = (Xf-X0)/n
        X_uni = zeros(n*nGaussLegendre)
        
        for i in range(n):
            X1 = X0 + i*hx ; X2 = X0 + (i+1)*hx
            X_uni[i*nGaussLegendre : (i+1)*nGaussLegendre] = ((X2-X1)/2)*xi + (X1+X2)/2
            
        WX = array(list(we)*n)
        WX = WX * X_uni * pi * hx
    
    
    if m == 0:
        Z_uni = array([Z0])
        WZ = 1.0
    else:
        hz = (Zf-Z0)/m
        Z_uni = zeros(m*nGaussLegendre)
        
        for j in range(m):
            Z1 = Z0 + j*hz ; Z2 = Z0 + (j+1)*hz
            Z_uni[j*nGaussLegendre : (j+1)*nGaussLegendre] = ((Z2-Z1)/2)*xi + (Z1+Z2)/2
        
        WZ = array(list(we)*m)/(2*m)
    

    W = outer(WZ, WX); W = W.flatten()
    
    X, Z = meshgrid(X_uni, Z_uni, indexing="xy"); X = X.flatten(); Z = Z.flatten()

# 
# A COMPLETER / MODIFIER
#
  
    return [X,Z,W]


#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
#
# ------------------------------------------------------------------------------------
#
# Classe du projet P2 
# Un détecteur de métal...
#


class MetalDetectorProject:
  mu0     = 4e-7*pi          # permeabilité du vide en [H/m] 
  muR     = 1-6.4e-6 #1.00022          # permeabilité relative de l'aluminium

  Rplate  = 0.238e-2 #2.425e-2/2           # rayon de la plaque [m]
  Hplate  = 2.425e-2/2 #0.238e-2          # épaisseur de la plaque [m]
  Zplate  = Hplate/2 #0.5e-2           # position verticale de la plaque [m]
                         
  Rcoil   = 1.7e-2/2          # rayon de la bobine [m]
  Hcoil   = 1.4e-2           # épaisseur de la bobine [m]
  Zcoil   = -0.5e-2          # position verticale de la bobine [m]
  nSpires = 200              # nombre de spires
  I       = 0.1              # courant dans la bobine [A]
  f       = 200000             # fréquence du courant [kHz]
    
  nZcoil  = 20               # discrétisation vertical de la bobine
  nTheta  = 8                # discrétisation du cercle
  Theta   = linspace(0,2*pi,nTheta+1)[:-1] + pi / nTheta
  Xcircle = cos(Theta)
  Ycircle = sin(Theta)
  
  def streamPlot(self,X,Z,Bx,Bz,colorStream,colorIplate):  
    plt.streamplot(X,Z,Bx,Bz,density=1.0,linewidth=None,color=colorStream)
    x = array([-self.Rcoil,self.Rcoil])
    z = repeat([linspace(0,-self.Hcoil,self.nZcoil)],2,axis=0)
    plt.plot(x,z,"o-r",linewidth=2)
    r = self.Rplate; h = self.Hplate/2
    x = array([-r, r, r,-r,-r])
    z = array([-h,-h, h, h,-h]) + self.Zplate
    plt.fill(x,z,facecolor='gray')
    plt.xlim((-0.03,0.03)); plt.ylim((-0.03,0.03))
    x = linspace(-r+h/2,r-h/2,30)
    z = self.Zplate * ones_like(x)
    plt.plot(x,z,'or',color=colorIplate)

def main():
  
# ------------------------------------------------------------------------------------ 
#
# Script de test
#
#
# -0- Initialisation du projet
#
# ------------------------------------------------------------------------------------

  p = MetalDetectorProject() 

# ------------------------------------------------------------------------------------
#
# -1- Calcul de l'inductance, inductance mutuelle, inductance induite
#     et inductance effective pour Mr. Oestges !
#
# ------------------------------------------------------------------------------------
  
  Zcoil = -linspace(0,p.Hcoil,p.nZcoil)
  Rcoil = p.Rcoil * ones_like(Zcoil)
  Icoil = p.I     * ones_like(Zcoil) * p.nSpires / p.nZcoil  

   
# 
# -1.1- Flux dans la bobine du champs induit par la bobine
#

  nX = 10; nZ = 5; nGL = 2
  X,Z,W = inductanceGaussLegendre(0,p.Rcoil,0,-p.Hcoil,nX,nZ,nGL)
  Bz    = inductanceMegneticField(X,Z,Rcoil,Zcoil,Icoil,p)[1]
  flux  = W @ Bz
  L     = p.nSpires * (flux/p.I)
  print("Flux across bobine         = %14.7e [Tesla m2]" % flux)
  print("Inductance                 = %14.7e [Henry]" % L)
  
# 
# -1.2- Flux dans la plaque du champs induit par la bobine
#
 
  nX = 20; nZ = 0; nGL = 1
  X,Z,W = inductanceGaussLegendre(0,p.Rplate,p.Zplate,0,nX,nZ,nGL)
  Bz    = inductanceMegneticField(X,Z,Rcoil,Zcoil,Icoil,p)[1]
  flux  = W @ Bz
  M     = p.nSpires * (flux/p.I)
  print("Flux across plate          = %14.7e [Tesla m2]" % flux)
  print("Mutual inductance          = %14.7e [Henry]" % M)

# 
# -1.3- Courants de Foucault 
#
   
  Iplate = -(p.Hplate) / (pi * p.mu0 * p.muR * X[1:]**2)
  for i in range(nX-1):
    Iplate[i] *= (W[0:i+2] @ Bz[0:i+2])  
    
  Rplate = X[1:]
  Zplate = p.Zplate * ones_like(Iplate)
   
# 
# -1.4- Flux dans la bobine du champs induit par les courants de la plaque métallique 
#

  nX = 10; nZ = 5; nGL = 2
  X,Z,W = inductanceGaussLegendre(0,p.Rcoil,0,-p.Hcoil,nX,nZ,nGL)
  Bz    = inductanceMegneticField(X,Z,Rplate,Zplate,Iplate,p)[1]
  flux  = W @ Bz
  Linduced = p.nSpires * (flux/p.I)
  print("Induced flux across bobine = %14.7e [Tesla m2]" % flux)
  print("Induced inductance         = %14.7e [Henry]" % Linduced)
  print("Effective inductance       = %14.7e [Henry]" % (L + Linduced))

# ------------------------------------------------------------------------------------
#
# -2- Représentation des champs magnétiques
#
# ------------------------------------------------------------------------------------

  plt.rcParams['toolbar'] = 'None'
   
  n = 20
  X,Z = meshgrid(linspace(-0.03,0.03,n),linspace(-0.03,0.03,n))
  Bcoilx,Bcoilz = inductanceMegneticField(X,Z,Rcoil,Zcoil,Icoil,p)    
  Bplatex,Bplatez = inductanceMegneticField(X,Z,Rplate,Zplate,Iplate,p)    
 
  plt.figure("A metal finder :-)",figsize=(10, 10))
  p.streamPlot(X,Z,Bcoilx,Bcoilz,'blue','gray')
  #plt.savefig("champ bobine.pdf")
  #plt.title('Electromagnetic field generated by the coil')
  
  plt.figure("Magnetic field from eddy currents..",figsize=(10, 10))
  p.streamPlot(X,Z,Bplatex,Bplatez,'red','blue')
  #plt.savefig("champ plaque.pdf")
  #plt.title('Electromagnetic field induced by the eddy currents of the plate')
 
  plt.figure("Both fields...",figsize=(10, 10))
  p.streamPlot(X,Z,Bcoilx+Bplatex,Bcoilz+Bplatez,'blue','red')
  #plt.savefig("champs superposés.pdf")
  #plt.title('The two superposed fields')
 
# ------------------------------------------------------------------------------------
#
# -3- Bz dans la bobine et la plaque
#
# ------------------------------------------------------------------------------------

  plt.figure("Champ magnétique : Bz dans la bobine et la plaque") 
  
  n = 20

  X = linspace(0,p.Rplate,n) 
  Z = p.Zplate * ones_like(X)
  Bz = inductanceMegneticField(X,Z,Rcoil,Zcoil,Icoil,p)[1]
  plt.plot(X,Bz,'-b')
  
  X = linspace(0,p.Rcoil,n) 
  Z = p.Zcoil * ones_like(X)
  Bz = inductanceMegneticField(X,Z,Rcoil,Zcoil,Icoil,p)[1]
  plt.plot(X,Bz,'-b')
  Bz = inductanceMegneticField(X,Z,Rplate,Zplate,Iplate,p)[1]
  plt.plot(X,Bz,'-r')
  
  Z = (p.Zcoil - p.Hcoil/2) * ones_like(X)
  Bz = inductanceMegneticField(X,Z,Rcoil,Zcoil,Icoil,p)[1]
  plt.plot(X,Bz,'-b')
  Bz = inductanceMegneticField(X,Z,Rplate,Zplate,Iplate,p)[1]
  plt.plot(X,Bz,'-r')

  Z = (p.Zcoil + p.Hcoil/2) * ones_like(X)
  Bz = inductanceMegneticField(X,Z,Rcoil,Zcoil,Icoil,p)[1]
  plt.plot(X,Bz,'-b')
  Bz = inductanceMegneticField(X,Z,Rplate,Zplate,Iplate,p)[1]
  plt.plot(X,Bz,'-r')
  plt.grid()
  
  
# ------------------------------------------------------------------------------------
#
# -4- Points d'intégration !
#
# ------------------------------------------------------------------------------------

  
  plt.figure("Integration nodes")
  
  nX = 5; nZ = 4; nGL = 3
  
  X,Z,W = inductanceGaussLegendre(0,p.Rcoil,0,-p.Hcoil,nX,nZ,nGL) 
  plt.plot(X,Z,'ob',markersize=5)
  
  x = array([0,p.Rcoil])
  z = repeat([linspace(0,-p.Hcoil,nZ+1)],2,axis=0)
  plt.plot(x,z,'-k')
  x = array([-0.003,p.Rcoil])
  z = repeat([linspace(0,-p.Hcoil,nZ+1)],2,axis=0)
  plt.plot(x,z,'--k')  
  z = array([0,-p.Hcoil])
  x = repeat([linspace(0,p.Rcoil,nX+1)],2,axis=0)
  plt.plot(x,z,'-k')
  z = array([0.003,-p.Hcoil])
  x = repeat([linspace(0,p.Rcoil,nX+1)],2,axis=0)
  plt.plot(x,z,'--k')

  X,Z,W = inductanceGaussLegendre(0,p.Rcoil,0.002,0.002,nX,0,nGL)
  plt.plot(X,Z,'or',markersize=5)
  X,Z,W = inductanceGaussLegendre(-0.002,-0.002,0,-p.Hcoil,0,nZ,nGL)
  plt.plot(X,Z,'or',markersize=5)
   
  plt.axis('equal')
  plt.axis('off')
  plt.show()


main()