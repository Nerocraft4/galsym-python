
#gettime?

#variables globals
global barra, disco, bulge, halo, parsb, ctes_rk78, sect, sectx

import numpy as np
import os
from maths.puntequil import puntequil
from maths.helpers import matriz_rk78
from maths.DF import DF
from utils.initializer import initializer
from utils.namers import name_datadir, name_ini_guess_pequi_file
from utils.io import read_ini_peqs_file
from models.other import centro_masas_halo, derFdelta

#TODO temp
from matplotlib import pyplot as plt
import matplotlib.animation as animation

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
    #TODO Check this with isopotential curve / zerovel curves
    #TODO check tol 1e-08
    [peqL1,peqL2,peqL3,peqL4,peqL5] = puntequil(ini_peqs,barra,disco,bulge,halo,parsb,options)
    
    '''
    for k in [peqL1,peqL2,peqL3,peqL4,peqL5]:
        print(k)
        DFLk = DF(k[1],barra,disco,bulge,halo,parsb)
        eigvals_DFLk,eigvecs_DFLk = np.linalg.eig(DFLk)
        print([x for x in eigvals_DFLk])
    print(len(eigvals_DFL1))
    for i in range(len(eigvals_DFL1)):
        eigv = eigvecs_DFL1[i]
        print("eigenvec",i)
        print(eigvals_DFL1[i])
        for j in eigvecs_DFL1[i]:
            print(j)
    '''
    #TODO moure això a un arxiu de plots
    #TODO També, potencialment, donar l'opció a fer que sigui animat i guardar-ho com a vídeo
    h = 0.01
    for k in [peqL1,peqL2,peqL4,peqL5]:
        #continuation
        a = k[1]
        continuation = [a]
        for j in range(400):
            DFLk = DF(a,barra,disco,bulge,halo,parsb)
            eigvals_DFLk,eigvecs_DFLk = np.linalg.eig(DFLk)
            #all eigenvalues should be *almost purely imaginary
            Fk_delta = derFdelta(delta=xdbulge,xvec=k[1][:3],barra=barra,disco=disco,bulge=bulge,halo=halo,parsb=parsb)
            u_dot = np.linalg.solve(DFLk,-Fk_delta) #direction vector
            a = a - u_dot[:3]*h
            continuation.append(a)
        continuation = np.array(continuation).transpose()
        #print()
        plt.scatter(continuation[:][0],continuation[:][1],s=3)
        
        #print(u_dot)
    plt.show()

    
        

























