
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

'''Grid punts equilibri'''

print("\n\n\n")

if False:
    pcs = pcg(whichobject= "halo", whichparam= "yd",
          paramfrom= 0, paramto= -6, cjacfrom= 0, cjacto= 1, density= 50, 
          ini_peqs= ini_peqs, galparams= galparams, displacements=displacements,
          solveroptions= options, point_evolution=True)   

'''Càlcul de punts d'equilibri a partir d'intent inicial'''

galparamslist = extract_galparams(galparams)
peqs = puntequil(ini_peqs,galparamslist,options)

DFpeqs = [DF(peqLi,galparamslist) for peqLi in peqs]
eigens = [eig(DFpeq) for DFpeq in DFpeqs]

for i in range(len(eigens)):
    eig = eigens[i]
    print(i+1,eig.eigenvalues)
    character_of_eq_point(eig)


if False:
    continuations = []
    func2vals = []
    for peq in peqs:
            continuation, func2val = delta_cont(initial_point=peq, initial_delta=xdbulge, 
                                increment_delta=0.01, continuation_length=3, params=galparams)
            continuations.append(continuation)
            func2vals.append(func2val)
    
    plots.param_contin_all(continuations,func2vals)


if True:
    galparamslist = extract_galparams(galparams)
    plots.isopoten_cont(rad=8,dens=100,nlines=100, galparams=galparamslist)
    #plots.isopoten(rad=10,dens=100,galparams=galparamslist)
    plots.isodensi_cont(rad=15,dens=100,nlines=50, galparams=galparamslist)
    #plots.isodensi(rad=10,dens=100,galparams=galparamslist)
    #plots.isodensi_parts(rad=15,dens=100,galparams=galparamslist)

#TODO ESTABILIDAD DE CADA PUNTO, però això no calcula la estabilitat, no?
punt = 3 #provant aquest
xkk = 3e-2
paprox,times = aproxlineal(xvec=peqs[punt][1],xkk=xkk,params=galparams) #TODO ver alternativa codi Josep

#plt.scatter(paprox[0],paprox[1],c=times)

#guardar arxiu corresponent al número de punt i a l'xkk que s'ha pres
#guardem temps i punts per separat
arxol = namer.arxol(punt,xkk,arxiu)
path = mydatadir + "/" + arxol
save(path+"_ts",times)
save(path+"_ps",paprox)

#TODO currently working on refinamiento
if True:
    refinamiento(pequil=peqs[punt][1],paprox=paprox,times=times,
                 params=params,CAMP="",CJAC="",GRADC="",SECCIO="",GRADS="")




ctes_rk78 = matriz_rk78()



















