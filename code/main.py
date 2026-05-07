
#gettime?

#variables globals
global barra, disco, bulge, halo, parsb, ctes_rk78, sect, sectx

import numpy as np
import os
from maths.puntequil import puntequil
from maths.helpers import matriz_rk78
from maths.algebra import eigs
from maths.param_continuation import delta_cont
#from maths.DF import DF

from utils.initializer import initialize
from utils import namer
#from utils.namers import name_datadir, name_ini_guess_pequi_file #rename to make more compact
from utils.io import read_ini_peqs_file
from models.other import centro_masas_halo, derFdelta

from utils import plots

#initialize "globals"
base = "." #TODO cuidado amb això, potser automatitzar?
inputfolder = base+"/input/"
data = base+"/datos/"
fil = "1" #tbd
col = "3" #tbd
arxiu = "_Om"+fil+col+".dat"
arxi = inputfolder+"modMiyFer"+arxiu
barra, disco, bulge, halo, parsb = initialize(arxi)

#configure options for aproxlineal, 
options = {"verbose":True,"tolerance":1e-8,"maxiter":300} 
#TODO maybe set from a .config file
#TODO verbose: False currently breaks puntequil? 
#llegir valors de la barra i del disc a partir d'un .txt

xd = [0.1] #make sure it is a list
xdhalo = 0 
ydhalo = 0 #max -4
# o desplaçar halo en la direcció contraria a bulge
xydhalo = [xdhalo,ydhalo]

for indxd in range(len(xd)):
    xdbulge = xd[indxd] #REFORMAT THIS S
    
    [xcm, ycm, zcm] = centro_masas_halo(xdbulge,xydhalo,barra,disco,bulge,halo)
    print(xcm,ycm,zcm)

    #PER A TOTA AQUESTA SECCIÓ, comentar i replantejar codi
    # tot això es podria fer amb vectors directament, tipus barra.db = despbar
    despbar = [-xcm,-ycm] 
    #WARN 3Dim, mirar orientacions que siguin correctes
    print(despbar)
    barra.xd = despbar[0]
    barra.yd = despbar[1]

    #idem #TODO make sure the displacement is consistent in 3D
    #TODO density curves will help check if this translations are correct
    despesf = [xdbulge,0]
    bulge.xd = barra.xd + despesf[0]
    bulge.yd = barra.yd + despesf[1]
    
    #idem
    despdisc = [0,0]
    disco.xd = despdisc[0]
    disco.yd = despdisc[1]

    #idem
    desphalo = xydhalo #REW var names and repetition
    halo.xd = desphalo[0]
    halo.yd = desphalo[1]
    
    params = [barra, disco, bulge, halo, parsb]

    mydatadir = data + namer.datadir(xydhalo, xdbulge, halo, bulge)
    
    try:
        print("creating new data directory in",mydatadir)
        os.makedirs(mydatadir)
    except:
        print("directory already exists")
    
    try:
        guess_pequi_file = namer.ini_guess_pequi_file(halo,xdbulge)
        print("reading initial guess for eq points from",guess_pequi_file)
        ini_peqs = read_ini_peqs_file(inputfolder+guess_pequi_file)
    except:
        print("no initial guess for eq. points provided for this case")

    ctes_rk78 = matriz_rk78()

    '''Càlcul de punts d'equilibri a partir d'intent inicial'''
    peqs = puntequil(ini_peqs,barra,disco,bulge,halo,parsb,options)
    eigens = [eigs(peqLi,params) for peqLi in peqs]

    continuations = []
    func2vals = []
    for peq in peqs:
            continuation, func2val = delta_cont(initial_point=peq[1], initial_delta=xdbulge, 
                                increment_delta=0.01, continuation_length=3, params=params)
            continuations.append(continuation)
            func2vals.append(func2val)
    
    plots.param_contin_all(continuations,func2vals)

    plots.isopoten(rad=10,dens=100,params=params)
    plots.isodensi(rad=10,dens=100,params=params)
    plots.isodensi_parts(rad=15,dens=100,params=params)
    
    
        

























