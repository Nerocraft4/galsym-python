from matplotlib import pyplot as plt
import numpy as np
import matplotlib.animation as animation
from maths.helpers import isopotencial, isodensidad, isodensidad_all

def isopoten(rad,dens,params):
    [barra,disco,bulge,halo,parsb] = params
    curva = isopotencial(rad, dens, barra, disco, bulge, halo, parsb)
    curva = np.array(curva).transpose()
    sc = plt.scatter(curva[0],curva[1],s=3,c=curva[3])
    plt.gca().set_aspect('equal')
    plt.title("isopotencial, zoom")
    plt.colorbar(sc)
    plt.axhline(0, linestyle='--',c="black")
    plt.axvline(0, linestyle='--',c="black")
    plt.show()

def isodensi(rad,dens,params):
    [barra,disco,bulge,halo,parsb] = params
    curva = isodensidad(rad, dens, barra, disco, bulge, halo)
    curva = np.array(curva).transpose()
    sc = plt.scatter(curva[0],curva[1],s=3,c=curva[3])
    plt.gca().set_aspect('equal')
    plt.title("isodensidad")
    plt.colorbar(sc)
    plt.axhline(0, linestyle='--',c="black")
    plt.axvline(0, linestyle='--',c="black")
    plt.show()

def isodensi_parts(rad,dens,params):
    [barra,disco,bulge,halo,parsb] = params
    fig, axs = plt.subplots(2,2)
    tts = ["barra", "disco", "bulge", "halo"]
    curva = isodensidad_all(15, 100, barra, disco, bulge, halo)
    curva = np.array(curva).transpose()
    plt.gca().set_aspect('equal')
    for i in [0,1]:
        for j in [0,1]:
            t = axs[i,j].scatter(curva[0],curva[1],s=3,c=curva[4+2*i+j])
            axs[i,j].set_title(tts[2*i+j])
            axs[i,j].set_aspect('equal')
            plt.colorbar(t,ax=axs[i,j])
    fig.suptitle("densidades de cada parte")
    plt.show()

def param_contin_all(continuations,func2vals):
    for continuation in continuations:
        cs = list(range(len(continuation[0])))
        plt.scatter(continuation[:][0],continuation[:][1],s=3,c=cs)
    plt.show()
    for func2val in func2vals: #TODO currently broken
        plt.plot(func2vals)
    plt.legend(["L1","L2","L3","L4","L5"])
    plt.title("Error evolution when continuating eq points")
    plt.xlabel("Delta Xbulge, $10^3=1kpc$")
    plt.ylabel("Norm of the amended potential")
    plt.show()
