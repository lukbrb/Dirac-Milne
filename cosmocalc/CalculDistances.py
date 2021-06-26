import numpy as np
import scipy.integrate as itg

class CalculDistances:
  def __init__(self, H_0, w_m, w_r, w_l):
    self.H_0 = H_0
    self.w_m = w_m
    self.w_r = w_r
    self.w_l = w_l
    self.w_k = 1 - (self.w_m + self.w_l + self.w_r)
    
    # -------- Unités-----------------------
    self.c = 299_792.458 #Vitesse de la lumière en km/s
    self.Gyr = 3600 * 24 * 365 * 1e9 # Un milliard d'années en secondes
    self.pars = 3.0856775 * 1e13 # Valeur d'un parsec en km
    #Conversion d'unité, unité de base : Mpc 
    self.unite_conv = {"Gpc": 1e-3, "Mpc": 1, "pc": 1e6, "km": self.pars * 1e6, "m": self.pars * 1e9}
  def integrande_DCMR(self, a):
    a_point = np.sqrt(self.w_k + self.w_m/a + self.w_r/a**2 + self.w_l * a**2)
    return 1/(a * a_point)


  def integrande_DTT(self, a):
    a_point = np.sqrt(self.w_k + self.w_m/a + self.w_r/a**2 + self.w_l * a**2)
    return 1/a_point


  def J(self, x):
    x = abs(x)
    if x > 1e-6:
      return np.sin(np.sqrt(x)/np.sqrt(x))
    else:
      return 1 + x/6 + x**2/120 + x**3/5040  #... x^n/(2n + 1)!
  

  def calcul_Z(self, z):
    Z, _ = itg.quad(self.integrande_DCMR, 1/(1 + z), 1) 
    # Le deuxième élément du tuple donne l'erreur de calcul
    return Z


  def calcul_DCMR(self, z):
    Z = self.calcul_Z(z)
    return Z * (self.c/self.H_0)
  

  def calcul_DTT(self, z):
    sol_DTT, _ = itg.quad(self.integrande_DTT, 1/(1 + z), 1)
    # Le deuxième élément du tuple donne l'erreur de calcul
    return sol_DTT * (self.c/self.H_0)


  def calcul_TVL(self, z):
    DTT = self.calcul_DTT(z)
    return (1/self.Gyr)* (DTT * 1e6 * self.pars)/self.c 


  def calcul_DA(self, z):
    Z = self.calcul_Z(z)
    return self.J(self.w_k * Z**2) * self.calcul_DCMR(z)/(1+z)


  def calcul_DL(self, z):
    DA = self.calcul_DA(z)
    return (1 + z)**2 * DA


  def calcul_COVOL(self, z):
    DCMR = self.calcul_DCMR(z)
    return (4/3) * np.pi * DCMR**3


  def calcul_VP(self, z):
    DA = self.calcul_DA(z)
    return (4/3) * np.pi * DA**3


  def description(self, z):
    print(f"Distance comobile radiale à z = {z}: {self.calcul_DCMR(z)} Mpc")
    print(f"Distance de diamètre angulaire à z = {z}: {self.calcul_DA(z)} Mpc")      
    print(f"Distance de luminosité à z = {z}: {self.calcul_DL(z)} Mpc") 
    print(f"Covolume = {self.calcul_COVOL(z)/1e9} Gpc^3")
    print(f"Volume propre = {self.calcul_VP(z)/1e9} Gpc^3 \n")
    return " "


class DistancesDM:
  
  def __init__(self, H_0):
    self.H_0 = H_0
        # -------- Unités-----------------------
    self.c = 299_792.458 #Vitesse de la lumière en km/s
    self.Gyr = 3600 * 24 * 365 * 1e9 # Un milliard d'années en secondes
    self.pars = 3.0856775 * 1e13 # Valeur d'un parsec en km
    #Conversion d'unité, unité de base : Mpc 
    self.unite_conv = {"Gpc": 1e-3, "Mpc": 1, "pc": 1e6, "km": self.pars * 1e6, "m": self.pars * 1e9}

  def calcul_DL(self, z, unite="Mpc"):
    return (self.c/self.H_0) * (1 + z) * np.sinh(np.log(1 + z)) * self.unite_conv[unite]


  def calcul_DCMR(self, z, unite="Mpc"):
    return (self.c/self.H_0) * np.log(1 + z) * self.unite_conv[unite]

  
  def calcul_DCMT(self, z, unite="Mpc"):
    """ Distance comobile transversale"""
    return (self.c/self.H_0) * np.sinh(np.log(1 + z)) * self.unite_conv[unite]

  def calcul_DA(self, z, unite="Mpc"):
    return self.calcul_DL(z) / (1 + z)**2 * self.unite_conv[unite]


  def calcul_COVOL(self, z, unite="Gpc"):
    f = np.sinh(np.log(1 + z))
    racine = np.sqrt(1 + f**2)
    K = 2 * np.pi * (self.c/self.H_0)**3
    return K * (f * racine - np.log(1 + z)) * self.unite_conv[unite]**3


  def calcul_VP(self, z, unite="Gpc"):
    DA = self.calcul_DA(z, unite)
    return (4/3) * np.pi * DA**3 
  

  def description(self, z, unite="Mpc"):
    print("Données pour l'univers de Dirac-Milne :\n")
    print(f"Distance comobile radiale à z = {z}: {self.calcul_DCMR(z, unite)} {unite}")
    print(f"Distance de diamètre angulaire à z = {z}: {self.calcul_DA(z, unite)} {unite}")      
    print(f"Distance de luminosité à z = {z}: {self.calcul_DL(z, unite)} {unite}") 
    print(f"Covolume = {self.calcul_COVOL(z, 'Gpc')} Gpc^3")
    print(f"Volume propre = {self.calcul_VP(z, 'Gpc')} Gpc^3\n") 
    print("Volumes toujours exprimés en Gpc^3 dans la méthode description")
    return str()
  
