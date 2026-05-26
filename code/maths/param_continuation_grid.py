from numpy import linspace, array
from maths.puntequil import puntequil
from utils.io import extract_galparams
from models.other import update
import matplotlib.pyplot as plt
import copy

def pcg(whichobject: str, whichparam: str, 
        paramfrom: float, paramto: float, cjacfrom: float, cjacto: float,  
        density: int, ini_peqs: list, galparams: list, displacements: list,
        solveroptions: dict, point_evolution: bool = False):
    '''
    Param Continuation Grid:
    This function calculates each combination of param & jacobi constant to return periodic orbits
    As a first implementation, only the continuation along parameters will be used.
    '''
    #assert paramto>paramfrom, "paramto must be greater than paramfrom"
    #assert cjacto>cjacfrom, "cjacto must be greater than cjacfrom"

    paramlist = linspace(paramfrom,paramto,density)
    cjaclist = linspace(paramfrom,paramto,density)
    pequils = [copy.deepcopy(ini_peqs)]    

    print("Studying",whichobject,"along parameter",whichparam,"from",paramfrom,"to",paramto)
    obj = galparams[whichobject]

    #first compute eq points corresponding to paramlist[1], using ini_peqs as an initial guess
    update(galparams,displacements,whichobject,whichparam,paramlist[1])
    new_peqs = puntequil(ini_peqs,extract_galparams(galparams),solveroptions)
    pequils.append(copy.deepcopy(new_peqs))
    directions = new_peqs-ini_peqs

    for i in range(2,density):
        p = paramlist[i]
        old_pecs = pequils[i-1]
        update(galparams,displacements,whichobject,whichparam,p)
        guess = old_pecs+directions
        new_peqs = puntequil(ini_guess=guess,
                            galparams=extract_galparams(galparams),options=solveroptions)
        pequils.append(copy.deepcopy(new_peqs))

        #print(guess[0],new_peqs[0])
        directions = new_peqs-old_pecs

    pequils = array(pequils)
    if point_evolution:
        pcs = copy.deepcopy(pequils)
        pcs = pcs.transpose([1,2,0])
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for p in pcs:
            ax.scatter(p[0][1:],p[1][1:],p[2][1:],s=1,c=paramlist[1:])
        plt.show()
        for i in range(5):
            p = pcs[i]
            plt.scatter(paramlist[1:],p[0][1:],s=1,c=paramlist[1:])
            plt.xlabel("Paràmetre")
            plt.ylabel("Coordenada X del punt Lagrangià")
            plt.title(str("Evolució de L"+str(i+1)+" versus "+whichobject+" "+whichparam))
            plt.show()
    return pequils


    
