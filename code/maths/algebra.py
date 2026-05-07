from .DF import DF
import numpy as np

def eigs(peqLi,params):
    [barra,disco,bulge,halo,parsb] = params #TODO repackage
    DFLk = DF(peqLi[1],barra,disco,bulge,halo,parsb)
    return np.linalg.eig(DFLk)

