import numpy as np
from scipy.interpolate import interp1d
from .DistancesDM import *


rho_0 = 5.35232 * 1e-27
H_0 = 70
pars = 3.0856775 * 1e13 # Valeur d'un parsec en km
#Conversion d'unité, unité de base : Mpc 
unite_conv = {"Gpc": 1e-3, "Mpc": 1, "pc": 1e6, "km": pars * 1e6, "m": pars * 1e9}
    
    
def rho(z):
    """ Evolution de la densité de matière en fonction du redshift"""
    return rho_0 * (1 + z)**3
    
    
def calcul_L(k, z_data, z_tab, unite="m"):
    """ Fonction qui interpole les valeurs de k et qui renvoie L"""
    f = interp1d(z_data, k)
    return np.pi/f(z_tab) * unite_conv[unite]


def volume_moyen(L):
    """ Volume moyen d'une structure. L en mètres"""
    return 4/3 * np.pi * L**3


def masse_moyenne(L):
    """Masse moyenne d'une structure, L en mètres"""
    return volume_moyen(L) * rho_0


def nombre_moyen(z, L):
    """ Nombre moyen de strcutures à redshift
        donné. Prend z et L, en mètres, en argument."""

    CV = COVOL(z, "m")
    dCV = np.diff(CV)
    v_moy = volume_moyen(L)
    return 0.5 * dCV / v_moy[:-1]


def sample_sphere(npoints, r=1, ndim=3):
    """ Projection aléatoire de points à la surface d'une sphère.
    Prend en argument un nombre de points à projeter, et le rayon de 
    la sphère."""
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec * r
