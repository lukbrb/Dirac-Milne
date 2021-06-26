import numpy as np
from CalculDistances import DistancesDM
univ = DistancesDM(70)

with open("donnees/coord_Mpic_Lpic/coord_structures_Mpic_Lpic.txt", "w") as doc:
    for z in np.geomspace(0.05, 20, 40):
        r = univ.calcul_DCMR(z, "Mpc")
        doc.write(f"{str(z)}\t{str(r)}\n")
