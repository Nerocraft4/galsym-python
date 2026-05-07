from matplotlib import pyplot as plt
import matplotlib.animation as animation
from .func2 import func2
from .DF import DF
from models.other import centro_masas_halo, derFdelta
import numpy as np

#TODO També, potencialment, donar l'opció a fer que sigui animat i guardar-ho com a vídeo 
def delta_cont(initial_point, initial_delta, increment_delta, continuation_length, params):
    [barra,disco,bulge,halo,parsb] = params 
    increment_delta = 0.001
    continuation_length = 7 #kpcs
    k = initial_point
    #continuation
    a = k
    continuation = [a]
    delta = initial_delta
    func2val = []
    iterations = int(continuation_length/increment_delta)
    for j in range(iterations):
        DFLk = DF(a,barra,disco,bulge,halo,parsb)
        #all eigenvalues should be purely imaginary
        
        #updating center of masses and all??
        [xcm, ycm, zcm] = centro_masas_halo(delta,[0,0],barra,disco,bulge,halo) #WARN XYDHALO SET TO 0,0 HERE
        #if j%50==0: print(xcm,ycm,zcm)
        despbar = [-xcm,-ycm] 
        barra.xd = despbar[0]
        barra.yd = despbar[1]
        despesf = [delta,0]
        bulge.xd = barra.xd + despesf[0]
        bulge.yd = barra.yd + despesf[1]
        
        Fk_delta = derFdelta(delta=delta,xvec=k[:3],
                             barra=barra,disco=disco,
                             bulge=bulge,halo=halo,parsb=parsb)
        u_dot = np.linalg.solve(DFLk,-Fk_delta) #direction vector
        a = a - u_dot[:3]*increment_delta
        func2val.append(np.linalg.norm(func2(a,barra,disco,bulge,halo,parsb)))
        delta = delta + increment_delta
        continuation.append(a)
    continuation = np.array(continuation).transpose()

    return continuation, func2val
