
#gettime?

#variables globals
global barra, disco, bulge, halo, parsb, ctes_rk78, sect, sectx

import numpy as np
import os
from maths.puntequil import puntequil
from maths.helpers import matriz_rk78, character_of_eq_point
#from maths.param_continuation import delta_cont #TODO DEPRECATED?
from maths.param_continuation_grid import pcg
from maths.aproxlineal import aproxlineal
from maths.DF import DF
from maths.refinamiento import refinamiento

from utils.initializer import initialize
from utils import namer, plots
from utils.io import read_ini_peqs, make_data_dir, extract_galparams

from models.other import setup, CTJAC

from numpy.linalg import eig
from numpy import save, load

#initialize "globals"
base = "." #TODO cuidado amb això, potser automatitzar?
inputfolder = base+"/input"
data = base+"/datos"
fil = "1" #TODO llegir de txt?
col = "3" #TODO idem
arxiu = "_Om"+fil+col
arxi = inputfolder+"/"+"modMiyFer"+arxiu+".dat"
galparams, displacements = initialize(arxi)#{"barra":barra,"disco":disco,"bulge":bulge,"halo":halo,"parsb":parsb}
options = {"verbose":True,"tolerance":1e-8,"maxiter":500} #configure options for aproxlineal,
#TODO maybe set from a .config file
#TODO verbose: False currently breaks puntequil?

galparams = setup(galparams,displacements)

gdgSec1 = lambda xvec : [xvec[1],[0,1,0,0,0,0]]
gdgSec2 = lambda xvec : [xvec[1],[0,1,0,0,0,0]]
gdgSec3 = lambda xvec : [xvec[1],[0,1,0,0,0,0]]
gdgSec4 = lambda xvec : [xvec[0],[1,0,0,0,0,0]]
gdgSec5 = lambda xvec : [xvec[0],[1,0,0,0,0,0]]

make_data_dir(data=data,params=galparams)

ini_peqs = read_ini_peqs(inputdir=inputfolder,params=galparams) #TODO podria estar dins de pcg

if True:
    '''Grid punts equilibri'''
    print("Calculant evolució dels punts d'equilibri")
    pcs = pcg(whichobject= "barra", whichparam= "eps", whichpoints= [4],
          whichsecs= [gdgSec4],
          paramfrom= 0, paramto= 0.1, cjacfrom= 0, cjacto= 1, density= 5,
          ini_peqs= ini_peqs, galparams= galparams, displacements=displacements,
          solveroptions= options, point_evolution=False)

'''Càlcul de punts d'equilibri a partir d'intent inicial'''
galparamslist = extract_galparams(galparams)
peqs = puntequil(ini_peqs,galparamslist,options)

DFpeqs = [DF(peqLi,galparamslist) for peqLi in peqs]
eigens = [eig(DFpeq) for DFpeq in DFpeqs]

'''Caràcter dels punts d'equilibri'''
for i in range(len(eigens)):
    eige = eigens[i]
    print("Point",i+1)
    print(peqs[i][:3])
    print(eige.eigenvalues)
    character_of_eq_point(eig)
    print("CJAC",CTJAC(xvec = peqs[i][:3], pvec = [0,0,0], params = galparamslist))

input()

'''Plots'''
if True:
    galparamslist = extract_galparams(galparams)
    #plots.isopoten_cont(rad=8,dens=100,nlines=100, galparams=galparamslist)
    #plots.isopoten(rad=10,dens=100,galparams=galparamslist)
    #plots.isodensi_cont(rad=15,dens=100,nlines=50, galparams=galparamslist)
    #plots.isodensi(rad=10,dens=100,galparams=galparamslist)
    #plots.isodensi_parts(rad=15,dens=100,galparams=galparamslist)

print("\n\n\n ORB PROP \n")

peqs = np.array(peqs)

if False:
    from maths import orbp
    from matplotlib import pyplot as plt
    import matplotlib.cm as cmap
    from itertools import chain
    X = []
    Y = []
    Z = []
    C = []
    #'''
    #L4, L5
    gdgSec = lambda xvec : [xvec[0],[1,0,0,0,0,0]] #g(xvec)=y, dg(xvec)=[0,1,0,0,0,0]
    cjacs = [-0.136,-0.135,-0.134]#,-0.133, -0.132, -0.131, -0.130, -0.129, -0.128, -0.127, -0.126]
    L = len(cjacs)
    for i in range(L):
        v = cjacs[i]
        temps,pos,cjac = orbp.compute_op(params=galparamslist,
                                         xi=[0,5.467,0e0,1e-3,0,0.e0],
                                         gdgSec=gdgSec,
                                         cjv=v)
        cjac = [cjac]*len(pos[0])
        X.append(pos[0])
        Y.append(pos[1])
        Z.append(pos[2])
        C.append([4]*len(pos[0]))
        temps,pos,cjac = orbp.compute_op(params=galparamslist,
                                         xi=[0,-5.467,0e0,-1e-3,0,0.e0],
                                         gdgSec=gdgSec,
                                         cjv=v)
        cjac = [cjac]*len(pos[0])
        X.append(pos[0])
        Y.append(pos[1])
        Z.append(pos[2])
        C.append([5]*len(pos[0]))
    #'''
    #L1
    gdgSec = lambda xvec : [xvec[1],[0,1,0,0,0,0]] #g(xvec)=y, dg(xvec)=[0,1,0,0,0,0]
    cjacs = [-0.3326,-0.335]
    L = len(cjacs)
    for i in range(L):
        v = cjacs[i]
        temps,pos,cjac = orbp.compute_op(params=galparamslist,
                                         xi=[0.1,0,0,0,1e-4,0.e0],
                                         gdgSec=gdgSec,
                                         cjv=v)
        cjac = [cjac]*len(pos[0])
        X.append(pos[0])
        Y.append(pos[1])
        Z.append(pos[2])
        C.append([1]*len(pos[0]))

    #L2
    gdgSec = lambda xvec : [xvec[1],[0,1,0,0,0,0]] #g(xvec)=y, dg(xvec)=[0,1,0,0,0,0]
    cjacs = [-0.144,-0.143]
    vs = [1e-7] #[0.9e-7,1e-7,1.1e-7] rang observat convergent
    L = len(cjacs)
    for i in range(len(vs)):
        #v = cjacs[i]
        ydot = vs[i]
        print(ydot)
        temps,pos,cjac = orbp.compute_op(params=galparamslist,
                                         xi=[-6.176-0.001*i,0,0,0,ydot,0.e0],
                                         gdgSec=gdgSec,
                                         cjv=0)
        cjac = [cjac]*len(pos[0])
        X.append(pos[0])
        Y.append(pos[1])
        Z.append(pos[2])
        C.append([2]*len(pos[0]))
        temps,pos,cjac = orbp.compute_op(params=galparamslist,
                                         xi=[6.176+0.001*i,0,0,0,-ydot,0.e0],
                                         gdgSec=gdgSec,
                                         cjv=0)
        cjac = [cjac]*len(pos[0])
        X.append(pos[0])
        Y.append(pos[1])
        Z.append(pos[2])
        C.append([3]*len(pos[0]))
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    X = list(chain(*X))
    Y = list(chain(*Y))
    Z = list(chain(*Z))
    C = list(chain(*C))
    sp = ax.scatter(X,Y,Z,c=C,cmap=cmap.summer,s=1)
    peqs = peqs.transpose()
    pts = ax.scatter(peqs[0],peqs[1],peqs[2],s=1,c="red")
    #fig.colorbar(sp)
    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    ax.set_zlim(-1,1)
    ax.set_title("PO around $L_i$ in the case $omega=0$")
    ax.set_xlabel("X (kpc)")
    ax.set_ylabel("Y (kpc)")
    ax.set_zlabel("Z (kpc)")
    plt.show()

#TODO ESTABILIDAD DE CADA PUNTO, però això no calcula la estabilitat, no?
punt = 3 #provant aquest
xkk = 3e-2
paprox,times = aproxlineal(xvec=peqs[punt],xkk=xkk,params=galparamslist) #TODO ver alternativa codi Josep

#plt.scatter(paprox[0],paprox[1],c=times)

#guardar arxiu corresponent al número de punt i a l'xkk que s'ha pres
#guardem temps i punts per separat
arxol = namer.arxol(punt,xkk,arxiu)
path = data + "/" + arxol
save(path+"_ts",times)
save(path+"_ps",paprox)

#TODO currently working on refinamiento
if False:
    refinamiento(pequil=peqs[punt],paprox=paprox,times=times,
                 params=galparamslist,CAMP="",CJAC="",GRADC="",SECCIO="",GRADS="")

'''
#####
DEPRECAT
#####

## Continuació antiga respecte xdbulge
if False:
    continuations = []
    func2vals = []
    for peq in peqs:
            continuation, func2val = delta_cont(initial_point=peq, initial_delta=xdbulge,
                                increment_delta=0.01, continuation_length=3, params=galparams)
            continuations.append(continuation)
            func2vals.append(func2val)

    plots.param_contin_all(continuations,func2vals)
'''
