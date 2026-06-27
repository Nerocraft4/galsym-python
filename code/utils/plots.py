from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation
from maths.helpers import isopotencial, isodensidad, isodensidad_all

def isopoten(rad: float, dens: int, galparams: list) -> None:
    curva = isopotencial(rad, dens, galparams)
    curva = np.array(curva).transpose()
    sc = plt.scatter(curva[0],curva[1],s=3,c=curva[3])
    plt.gca().set_aspect('equal')
    plt.title("isopotencial, zoom")
    plt.colorbar(sc)
    plt.axhline(0, linestyle='--',c="black")
    plt.axvline(0, linestyle='--',c="black")
    plt.show()

def isopoten_cont(rad: float, dens: int, nlines: int, galparams: list) -> None:
    curva = isopotencial(rad, dens, galparams)
    curva = np.array(curva).transpose()
    X = np.linspace(-rad,rad,dens,False)
    Y = np.linspace(-rad,rad,dens,False)
    Z = curva[3].reshape(dens,dens).transpose()
    Za = np.geomspace(np.min(Z),np.max(Z),nlines)#TODO passar param
    ylorbr = cm.get_cmap('plasma', nlines)
    colors = ylorbr(np.linspace(0,1,nlines))
    sc = plt.contour(X,Y,Z,Za,colors=colors)
    #plt.contour(X,Y,Z,levels=Za,colors="black") 
    plt.gca().set_aspect('equal')
    plt.title("isopotencial, zoom")
    plt.colorbar(sc)
    plt.axhline(0, linestyle='--',c="black")
    plt.axvline(0, linestyle='--',c="black")
    plt.show()

def isodensi(rad: float, dens: int, galparams: list) -> None:
    curva = isodensidad(rad, dens, galparams)
    curva = np.array(curva).transpose()
    sc = plt.scatter(curva[0],curva[1],s=3,c=curva[3])
    plt.gca().set_aspect('equal')
    plt.title("isodensidad")
    plt.colorbar(sc)
    plt.axhline(0, linestyle='--',c="black")
    plt.axvline(0, linestyle='--',c="black")
    plt.show()

def isodensi_cont(rad: float, dens: int, nlines: int, galparams: list) -> None:
    curva = isodensidad(rad, dens, galparams)
    curva = np.array(curva).transpose()
    X = np.linspace(-rad,rad,dens,False)
    Y = np.linspace(-rad,rad,dens,False)
    Z = curva[3].reshape(dens,dens).transpose()
    Za = np.geomspace(np.min(Z),np.max(Z),nlines)#TODO passar param
    ylorbr = cm.get_cmap('plasma', nlines)
    colors = ylorbr(np.linspace(0,1,nlines))
    sc = plt.contour(X,Y,Z,Za,colors=colors)
    plt.gca().set_aspect('equal')
    plt.title("isodensidad")
    plt.colorbar(sc)
    plt.axhline(0, linestyle='--',c="black")
    plt.axvline(0, linestyle='--',c="black")
    plt.show()

def isodensi_parts(rad: float, dens: int, galparams: list) -> None:
    fig, axs = plt.subplots(2,2)
    tts = ["barra", "disco", "bulge", "halo"]
    curva = isodensidad_all(rad, dens, galparams)
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

def param_contin_all(continuations: list, func2vals: list) -> None:
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

def plot_PO_3D(points,times,var):
    print(points.shape)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    sp = ax.scatter(points[0],points[1],points[2],c=times,s=1)
    #fig.colorbar(sp)
    #ax.set_xlim(-10,10)
    #ax.set_ylim(-10,10)
    #ax.set_zlim(-1,1)
    ax.set_title(str("PO around $L_i$ in the case $\epsilon=$"+str(round(var,3))))
    ax.set_xlabel("X (kpc)")
    ax.set_ylabel("Y (kpc)")
    ax.set_zlabel("Z (kpc)")
    plt.show()
