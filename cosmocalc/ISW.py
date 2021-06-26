import numpy as np
from scipy.interpolate import interp1d
from .DistancesDM import *


rho_0 = 5.35232 * 1e-27
H_0 = 70
pars = 3.0856775 * 1e13 # Valeur d'un parsec en km
#Conversion d'unité, unité de base : Mpc 
unite_conv = {"Gpc": 1e-3, "Mpc": 1, "pc": 1e6, "km": pars * 1e6, "m": pars * 1e9}
K = 4/3 * np.pi # Constante, pour les volumes
T_H0_MKSA = 977.8 * 1e9 * 3.156E7 
G = 0.667E-10
M_sol = 1.99E30
def rho(z):
    """ Evolution de la densité de matière en fonction du redshift"""
    return rho_0 * (1 + z)**3
    
    
def calcul_L(k, z_data, z_tab, unite="m"):
    """ Fonction qui interpole les valeurs de k et qui renvoie L selon 
        l'unitée choisie (m, km, pc, Mpc, Gpc) 
        --------------------------------------------------------------
        Arguments :
            - Le vecteur d'onde k
            - Les données z_data du redshift
            - Le tableau de valeurs z_tab pour effectuer l'interpolation
        ---------------------------------------------------------------
        Unités :
        
            m, km, pc, Mpc, Gpc"""

    f = interp1d(z_data, k)
    return np.pi/f(z_tab) * unite_conv[unite]


def volume_moyen(L):
    """ Volume moyen d'une structure. 
        Argument:
            - L (en mètres) """
    return 4/3 * np.pi * L**3


def masse_moyenne(L):
    """ Masse moyenne d'une structure, 
        Argument:
            - L (en mètres)"""
    return volume_moyen(L) * rho_0


def nombre_structures(z, L, alpha=0.5, H0=H_0):
    """ Nombre moyen de structures à redshift
        donné. 
        Arguments :
            - Redshift z 
            - Taille de la structure L (en mètres)
            - Alpha : proportion de matière dans un volume
            - H0 : Constante de Hubble actuelle (=70 par défaut)"""

    CV = COVOL(z, H0, "m")
    dCV = np.diff(CV)
    v_moy = volume_moyen(L)
    return alpha * dCV / v_moy[:-1]


def sample_sphere(npoints, r=1, ndim=3):
    """ Projection aléatoire de points à la surface d'une sphère.
    Prend en argument un nombre de points à projeter, et le rayon de 
    la sphère."""
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec * r


def pot_moyen(L, z, unite="m"):
    # Potentiel au bord de la zone de déplétion pour une structure typique de SDSS, corrigée pour le redshift
    return G * masse_moyenne(L) * (1 + z)/(L*unite_conv[unite])/(c * 1e3)**2

# Angle solide, angle en degrés et \ell associé pour cette structure

def angle_sol_moy(L, z, unite="m", H0=H_0):
    return 0.5 * (L/DA(z, H0, unite)/(1 + z))**2


def theta_moy(L, z , unite="m", H0=H_0):
    return np.sqrt(2*angle_sol_moy(L, z, unite, H0))* 180/np.pi


def ell(L, z, unite="m", H0=H_0):
    return  180/theta_moy(L, z, unite, H0)


def angle_sol_tot(L, z, unite="m"):
    """Angle solide total couvert par tranche de redhsift"""
    return angle_sol_moy(L, z, unite)[:-1] * nombre_structures(z, L)


def ISW_moy(z, L, unite="m", H0=H_0):
    """ Variation dmoyenne du potentiel due à l'expansion de l'univers 
        durant le temps de traversée de la structure par le photon """
    cross_time = unite_conv[unite] * L/(1 + z)/(c * 1e3)
    hubble_time = 1./H_0 * T_H0_MKSA/(1 + z)
    return cross_time * pot_moyen(L, z, unite) / hubble_time