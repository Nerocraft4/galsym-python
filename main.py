
#gettime?

#variables globals
global barra, disco, bulge, halo, parsb, ctes_rk78, sect, sectx

import numpy as np
import os
from .maths.puntequil import puntequil
from .initializer import initializer

#initialize "globals"
base = "pythonbase/"
inputfolder = base+"./input/"
fil = "1" #tbd
col = "3" #tbd
arxiu = "_Om"+fil+col+".dat"
arxi = inputfolder+"modMiyFer"+arxiu
barra, disco, bulge, halo, parsb = initializer(arxi)

#configure options for aproxlineal, #TODO maybe set from a .config file
options = {"verbose":True,"tolerance":1e-8,"maxiter":300} 
#TODO verbose: False currently breaks puntequil? 

def centro_masas_halo(xdbulge, xydhalo):
    md = disco.GM
    centrod = np.array([0,0,0])
    md = 0 #zero perque està centrat i el tenim com a referència
    
    mb = barra.GM
    centrob = np.array([0,0,0])

    mesf = bulge.GM
    centroesf = np.array([xdbulge,0,0])

    mh = halo.GM
    centroh = np.array([xydhalo[0], xydhalo[1],0])

    mt = md+mb+mesf+mh
    
    cm = (1/mt)*(md*centrod +mb* centrob + mesf* centroesf + mh * centroh)

    return cm

def nameit(xydhalo,xdbulge):
    p = r"halomass_"
    p+= str(halo.GM)
    p+= "_bulgemass_"
    p+= str(bulge.GM)
    p+= "_halo_xd_"
    p+= str(xydhalo[0])
    p+= "_yd_"
    p+= str(xydhalo[1])
    p+= r"/"
    p+= "bulgexd_"
    p+= str(xdbulge)
    p+= r"/"
    return p

def matriz_rk78():
    alpha= [0,    2/27,       1/9,        1/6,        5/12,   0.5,    5/6,        1/6,    2/3,    1/3,    1,  0,  1]
    beta = [0,          2/27,       1/36,       1/12,       1/24,   0,      1/8,        5/12,       0,      -25/16,
            25/16,      0.05,       0,          0,          0.25,   0.2,    -25/108,    0,          0,      125/108,
            -65/27,     250/108,    31/300,     0,          0,      0,      61/225,     -2/9,       13/900, 2,
            0,          0,          -53/6,      704/45,     -107/9, 67/90,  3,          -91/108,    0,      0,
            23/108,     -976/135,   311/54,     -19/60,     17/6,   -1/12,  2383/4100,  0,          0,      -341/164,
            4496/1025,  -301/82,    2133/4100,  45/82,      45/164, 18/41,  3/205,      0,          0,      0,
            0,          -6/41,      -3/205,     -3/41,      3/41,   6/41,   0,          -1777/4100, 0,  0
            -341/164,   4496/1025,  -289/82,    2193/4100,  51/82,  33/164, 12/41,      0,          1]
    c = [41/840,    0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 41/840]
    cp= [0,         0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0 , 41/840, 41/840]
    return [alpha, beta, c, cp]
    ctes_rk78.alfa = alpha
    ctes_rk78.beta = beta
    ctes_rk78.c = c
    ctes_rk78.cp = cp


#llegir valors de la barra i del disc a partir d'un .txt

data = "./datos/"

xd = [0.1] #make sure it is a list

xdhalo = 0
ydhalo = 0
xydhalo = [xdhalo,ydhalo]

for indxd in range(len(xd)):
    xdbulge = xd[indxd] #REFORMAT THIS S
    
    [xcm, ycm, zcm] = centro_masas_halo(xdbulge,xydhalo)
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
    
    datadir = data + nameit(xydhalo, xdbulge)
    
    try:
        print("creating new data directory in",datadir)
        os.makedirs(datadir) #revisar, segurament os.algo
    except:
        print("directory already exists")
    
    #READ INITIAL GUESS FOR EQ POINTS, maybe build a separate func for this
    try:
        a = ['peqi_entrada3','_halo_xd_',str(halo.xd),
         '_yd_',str(halo.yd),
         '_bulgexd_',str(xdbulge)]
        ini_guess_pequi = inputfolder+"".join(a).replace(".","_")+".dat"
        print(ini_guess_pequi)
        ini_peqs_data = ""
        with open(file=ini_guess_pequi, mode="r") as f:
            ini_peqs_data = f.readlines() #5 L_i punts de 6D cadascun (p,q)\in R^6
        ini_peqs_data = [x.strip().split(" ") for x in ini_peqs_data]
        ini_peqs = []
        for i in ini_peqs_data:
            peq_i = []
            for j in i:
                if j=="":
                    continue
                peq_i.append(float(j))
            ini_peqs.append(peq_i)
        
    except:
        print("no initial guess for eq. points provided for this case")

    ctes_rk78 = matriz_rk78()

    '''Càlcul de punts d'equilibri a partir d'intent inicial'''
    [peqL1,peqL2,peqL3,peqL4,peqL5] = puntequil(ini_peqs,barra,parsb,options)
       
    























