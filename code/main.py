
#gettime?

#variables globals
global barra, disco, bulge, halo, parsb, ctes_rk78, sect, sectx

import numpy as np
import os
from maths.puntequil import puntequil
from maths.helpers import matriz_rk78
from utils.initializer import initializer
from utils.namers import name_datadir, name_ini_guess_pequi_file
from utils.io import read_ini_peqs_file
from models.other import centro_masas_halo

#initialize "globals"
base = "." #TODO cuidado amb això, potser automatitzar?
inputfolder = base+"/input/"
data = base+"/datos/"
fil = "1" #tbd
col = "3" #tbd
arxiu = "_Om"+fil+col+".dat"
arxi = inputfolder+"modMiyFer"+arxiu
barra, disco, bulge, halo, parsb = initializer(arxi)

#configure options for aproxlineal, 
options = {"verbose":True,"tolerance":1e-8,"maxiter":300} 
#TODO maybe set from a .config file
#TODO verbose: False currently breaks puntequil? 
#llegir valors de la barra i del disc a partir d'un .txt

xd = [0.1] #make sure it is a list
xdhalo = 0
ydhalo = 0
xydhalo = [xdhalo,ydhalo]

for indxd in range(len(xd)):
    xdbulge = xd[indxd] #REFORMAT THIS S
    
    [xcm, ycm, zcm] = centro_masas_halo(xdbulge,xydhalo,barra,disco,bulge,halo)
    print(xcm,ycm,zcm)

    #PER A TOTA AQUESTA SECCIÓ, comentar i replantejar codi
    # tot això es podria fer amb vectors directament, tipus barra.db = despbar
    despbar = [-xcm,-ycm] #WARN 3Dim, mirar orientacions que siguin correctes
    print(despbar)
    barra.xd = despbar[0]
    barra.yd = despbar[1]

    #idem
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
    
    datadir = data + name_datadir(xydhalo, xdbulge, halo, bulge)
    
    try:
        print("creating new data directory in",datadir)
        os.makedirs(datadir) #revisar, segurament os.algo
    except:
        print("directory already exists")
    
    try:
        ini_guess_pequi_file = name_ini_guess_pequi_file(halo,xdbulge)
        print("reading initial guess for eq points from",ini_guess_pequi_file)
        ini_peqs = read_ini_peqs_file(inputfolder+ini_guess_pequi_file)
    except:
        print("no initial guess for eq. points provided for this case")

    ctes_rk78 = matriz_rk78()

    '''Càlcul de punts d'equilibri a partir d'intent inicial'''
    [peqL1,peqL2,peqL3,peqL4,peqL5] = puntequil(ini_peqs,barra,disco,bulge,halo,parsb,options)
       
    























