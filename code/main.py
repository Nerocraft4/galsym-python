
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

from models.other import setup

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

make_data_dir(data=data,params=galparams)

ini_peqs = read_ini_peqs(inputdir=inputfolder,params=galparams) #TODO podria estar dins de pcg

if False:
    '''Grid punts equilibri'''
    print("Calculant evolució dels punts d'equilibri")
    pcs = pcg(whichobject= "halo", whichparam= "yd",
          paramfrom= 0, paramto= -6, cjacfrom= 0, cjacto= 1, density= 50,
          ini_peqs= ini_peqs, galparams= galparams, displacements=displacements,
          solveroptions= options, point_evolution=True)

'''Càlcul de punts d'equilibri a partir d'intent inicial'''
galparamslist = extract_galparams(galparams)
peqs = puntequil(ini_peqs,galparamslist,options)

DFpeqs = [DF(peqLi,galparamslist) for peqLi in peqs]
eigens = [eig(DFpeq) for DFpeq in DFpeqs]

'''Caràcter dels punts d'equilibri'''
for i in range(len(eigens)):
    eig = eigens[i]
    print("Point",i+1)
    print(peqs[i][:3])
    print(eig.eigenvalues)
    character_of_eq_point(eig)

'''Plots'''
if True:
    galparamslist = extract_galparams(galparams)
    #plots.isopoten_cont(rad=8,dens=100,nlines=100, galparams=galparamslist)
    #plots.isopoten(rad=10,dens=100,galparams=galparamslist)
    #plots.isodensi_cont(rad=15,dens=100,nlines=50, galparams=galparamslist)
    #plots.isodensi(rad=10,dens=100,galparams=galparamslist)
    #plots.isodensi_parts(rad=15,dens=100,galparams=galparamslist)

print("\n\n\n ORB PROP \n")

if True:
    from maths import orbp
    from matplotlib import pyplot as plt
    for v in [5e-2,3e-2,9e-3,1e-4]:
        gdgSec = lambda xvec : [xvec[0],[1,0,0,0,0,0]] #g(xvec)=y, dg(xvec)=[0,1,0,0,0,0]
        temps,pos,cjac = orbp.compute_op(params=galparamslist,xi=[0,5.467,0e0,v,0,0.e0],gdgSec=gdgSec)
        print(cjac)
        plt.scatter(pos[0],pos[1],c=temps,s=1)
        plt.xlim(-10,10)
        plt.ylim(-10,10)
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
