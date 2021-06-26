import numpy as np

H_0 = 70
# -------- Unités-----------------------
c = 299_792.458 #Vitesse de la lumière en km/s
Gyr = 3600 * 24 * 365 * 1e9 # Un milliard d'années en secondes
pars = 3.0856775 * 1e13 # Valeur d'un parsec en km
#Conversion d'unité, unité de base : Mpc 
unite_conv = {"Gpc": 1e-3, "Mpc": 1, "pc": 1e6, "km": pars * 1e6, "m": pars * 1e9}

  
def DL(z, H0=H_0,  unite="Mpc"):
    """Distance de luminsoité"""
    return (c/float(H0)) * (1 + z) * np.sinh(np.log(1 + z)) * unite_conv[unite]


def DCMR(z, H0=H_0, unite="Mpc"):
    """Distance comobile radiale"""
    return (c/float(H0)) * np.log(1 + z) * unite_conv[unite]

  
def DCMT(z, H0=H_0, unite="Mpc"):
    """Distance comobile transversale"""
    return (c/H0) * np.sinh(np.log(1 + z)) * unite_conv[unite]

def DA(z, H0=H_0, unite="Mpc"):
    """Distance de diamètre angulaire"""
    return DL(z, H0, unite) / (1 + z)**2 


def COVOL(z, H0=H_0, unite="Gpc"):
    """Calcul du covolume"""
    f = np.sinh(np.log(1 + z))
    racine = np.sqrt(1 + f**2)
    K = 2 * np.pi * (c/H0)**3
    return K * (f * racine - np.log(1 + z)) * unite_conv[unite]**3


def VP(z, H0=H_0, unite="Gpc"):
    return (4/3) * np.pi * DA(z, H0, unite)**3 
  

def description(z, H0=H_0, unite="Mpc"):
    """ Renvoie les paramètres importants pour un redshift z.\n
    Volumes toujours exprimés en Gpc^3 dans la méthode description"""
    
    print(f"""Données pour l'univers de Dirac-Milne :\n
    Distance comobile radiale à z = {z}: {DCMR(z, H0, unite)} {unite}\n
    Distance de diamètre angulaire à z = {z}: {DA(z, H0, unite)} {unite}\n     
    Distance de luminosité à z = {z}: {DL(z, H0, unite)} {unite}\n 
    Covolume = {COVOL(z, H0, 'Gpc')} Gpc^3\n
    Volume propre = {VP(z, H0, 'Gpc')} Gpc^3\n""")
    return ''
    