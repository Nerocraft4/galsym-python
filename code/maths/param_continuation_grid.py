from models.other import update, CTJAC
from utils.io import extract_galparams
from utils.plots import plot_PO_3D
from maths.puntequil import puntequil
from maths.DF import DF
from maths import orbp

from numpy import linspace, array, real
from numpy.linalg import eig
import matplotlib.pyplot as plt
import copy

def get_indiv_orbit(peqs,delta,galparamslist,whichsecs,var,prop_cjac):
    DFpeqs = [DF(peqLi,galparamslist) for peqLi in peqs]
    eigens = [eig(DFpeq) for DFpeq in DFpeqs]
    for k in range(len(peqs)):
        peq = peqs[k]
        cjac = CTJAC(peq[:3],[0,0,0],galparamslist)
        eige = eigens[k]
        vals = eige.eigenvalues
        vecs = eige.eigenvectors
        u = vecs[0] #el que tingui tot part real
        xmoved = [peq[0],peq[1],peq[2],0,0,0]
        xmoved += delta*u #array([0.0,1.0,0.0,1.0,0.0,0.0])
        print("xmoved",xmoved)
        print("delta*u",delta*u)
        xmoved = real(xmoved)
        temps,pos,cjac = orbp.compute_op(params=galparamslist,
                                         xi=xmoved,
                                         gdgSec=whichsecs[k],
                                         cjv=prop_cjac)
        #plot_PO_3D(pos,temps,var)
        cjacfinal = CTJAC(xvec = [pos[0][0],pos[1][0],pos[2][0]], 
                         pvec = [pos[3][0],pos[4][0],pos[5][0]],
                         params = galparamslist)
        return {"pos":pos,"temps":temps,"val":var,"cjac":cjacfinal}

def plot_grid_data(data):
    plt.rcParams.update({'font.size': 6})
    M = len(data)
    L = len(data[0])
    fig, axes = plt.subplots(M,L,figsize=(2*L, 2*M))
    fig.suptitle("(Not to scale) Evolution of orbit inclination around L4,5 (XZ Plane projection) wrt $\epsilon$")
    for j in range(M):
        for i in range(L):
            ax = axes[j][i]
            X = data[j][i]["pos"]
            v = data[j][i]["val"]
            cj = data[j][i]["cjac"]
            print("lens",len(X),len(X[0]))
            ax.scatter(X[0],X[2],s=1)
            ax.set_title(str("$\epsilon=$"+str(round(v,2))+", $C_{Jac}$="+str(round(cj,6))))
            ax.set_xlabel("X [kPc]")
            ax.set_ylabel("Z [kPc]")
            ax.set_xlim(-0.75,0.75)
            #ax.set_ylim(-0.25+5.466,0.25+5.466)
            ax.set_ylim(-0.015,0.015)
            #ax.set_ylim(-0.75,0.75)
    fig.tight_layout()
    plt.show()

def pcg(whichobject: str, whichparam: str, whichpoints: list,
        whichsecs: list,
        paramfrom: float, paramto: float, cjacfrom: float, cjacto: float,  
        density: int, ini_peqs: list, galparams: dict, displacements: list,
        solveroptions: dict, point_evolution: bool = False):
    '''
    Param Continuation Grid:
    This function calculates each combination of param & jacobi constant to return periodic orbits
    As a first implementation, only the continuation along parameters will be used.

    whichpoints: list of ints representing the points that will be studied
    '''
    #assert paramto>paramfrom, "paramto must be greater than paramfrom"
    #assert cjacto>cjacfrom, "cjacto must be greater than cjacfrom"

    delta = 1e-5 #multiplier to the eigenvector displacement    

    og_pequils = copy.deepcopy(ini_peqs)
    selected_points = []

    grid_data = []

    print("Studying",whichobject,"along parameter",whichparam,"from",paramfrom,"to",paramto)
    if set(whichpoints) != set([1,2,3,4,5]):
        for v in whichpoints:
            selected_points.append(og_pequils[v-1])
    obj = galparams[whichobject]

    for j in range(5):#TODO hardcode
        print("cjac iter",j)
        line_data = []

        paramlist = linspace(paramfrom,paramto,density)
        update(galparams,displacements,whichobject,whichparam,paramlist[0])
        cjaclist = linspace(paramfrom,paramto,density)#TODO

        curr_prop_cjac = -0.1370+0.00005*j

        galparamslist = extract_galparams(galparams)

        pequils = [copy.deepcopy(selected_points)]
        for peq in selected_points:
            print("iter","0 CTJAC",CTJAC(xvec= peq, pvec= [0,0,0],params= galparamslist))
        
        line_data.append(get_indiv_orbit(selected_points,delta,galparamslist,
                                         whichsecs,paramfrom,curr_prop_cjac))
            

        ##########################################################################################
        #first compute eq points corresponding to paramlist[1], using ini_peqs as an initial guess
        update(galparams,displacements,whichobject,whichparam,paramlist[1])
        new_peqs = puntequil(ini_guessess = selected_points,
                             galparams = galparamslist,
                             options = solveroptions)
        for peq in new_peqs:
            print("iter","1 CTJAC",CTJAC(xvec= peq, pvec= [0,0,0],params= galparamslist))
        pequils.append(copy.deepcopy(new_peqs))
        directions = new_peqs-selected_points

        line_data.append(get_indiv_orbit(new_peqs,delta,galparamslist,
                                         whichsecs,paramlist[1],curr_prop_cjac))

        ##########################################################################################
        for i in range(2,density):
            p = paramlist[i]
            old_pecs = pequils[i-1]
            update(galparams,displacements,whichobject,whichparam,p)
            guess = old_pecs+directions
            new_peqs = puntequil(ini_guessess=guess,
                                galparams=galparamslist,
                                options=solveroptions)
            pequils.append(copy.deepcopy(new_peqs))
            for peq in new_peqs:
                print("iter",i,"CTJAC",CTJAC(xvec= peq, pvec= [0,0,0],params= galparamslist))
            #print(guess[0],new_peqs[0])
            directions = new_peqs-old_pecs

            line_data.append(get_indiv_orbit(new_peqs,delta,galparamslist,
                                             whichsecs,p,curr_prop_cjac))
        grid_data.append(line_data)

    pequils = array(pequils)
    if point_evolution:
        pcs = copy.deepcopy(pequils)
        pcs = pcs.transpose([1,2,0])
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for p in pcs:
            ax.scatter(p[0][1:],p[1][1:],p[2][1:],s=1,c=paramlist[1:])
        plt.show()
        for i in range(len(whichpoints)):
            p = pcs[i-1]
            plt.scatter(paramlist[1:],p[0][1:],s=1,c=paramlist[1:])
            plt.xlabel("Paràmetre")
            plt.ylabel("Coordenada X del punt Lagrangià")
            plt.title(str("Evolució de L"+str(whichpoints[i])+" versus "+whichobject+" "+whichparam))
            plt.show()
    
        

    plot_grid_data(grid_data)
    
    return pequils


    
