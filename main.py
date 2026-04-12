
#gettime?

#variables globals
global barra, disco, bulge, halo, parsb, ctes_rk78, sect, sectx

import numpy as np
import os
from maths.puntequil import puntequil

from initializer import initializer
barra, disco, bulge, halo, parsb = initializer("./input/modMiyFer_Om13.dat")

#definir com a objectes?

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

    print(md,mb,mesf,mh)
    mt = md+mb+mesf+mh
    
    print(mt)
    
    cm = (1/mt)*(md*centrod +mb* centrob + mesf* centroesf + mh * centroh) #caldrà tractar com a np arrays

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

dire = "./modMiyFer/"
data = "./datos/"
fil = "1" #tbd
col = "3" #tbd
arxiu = "_Om"+fil+col+".dat"
arxi = dire+"modMiyFer"+arxiu
arxi2 = ""


xd = [2.5] #make sure it is a list

xdhalo = 0
ydhalo = -2
xydhalo = [xdhalo,ydhalo]

for indxd in range(len(xd)):
    xdbulge = xd[indxd] #REFORMAT THIS S
    
    [xcm, ycm, zcm] = centro_masas_halo(xdbulge,xydhalo)
    
    #PER A TOTA AQUESTA SECCIÓ, comentar i replantejar codi
    # tot això es podria fer amb vectors directament, tipus barra.db = despbar
    despbar = [-xcm,-ycm] #WARN 3Dim, mirar orientacions que siguin correctes
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
    #os.path.joing(
    try:
        print("creating new data directory in",datadir)
        os.makedirs(datadir) #revisar, segurament os.algo
    except:
        print("directory already exists")

    ctes_rk78 = matriz_rk78()

    '''Càlcul de punts d'equilibri a partir d'intent inicial'''
    [peqL1,peqL2,peqL3,peqL4,peqL5] = puntequil(arxi2)
       
    























