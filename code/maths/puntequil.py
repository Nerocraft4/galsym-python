from scipy.optimize import fsolve
from .func2 import func2, func2jac
import numpy as np
from numpy.typing import NDArray

def puntequil(ini_guess: list, galparams: list, options: dict) -> NDArray:
    '''
    This function fins an equilibrium point in the field of the galaxy given galactic
    parameters and an initial guess, through Newton's method (fsolve)
    Input:
        ini_guess:  5x3 vector containing 5 3D initial guesses, one for each Lagr point
                    velocities are assumed to be zero.
        galparams:  list of objects containing galactic params
        options:    dictionary of settings for function solver
            verbose:    boolean, defafult false
            tolerance:  float, default 1e-08
            maxiter:    int, default 300
    Output:
        peqs:   5x3 array containing the equilibrium points of the system
    '''
    [barra,disco,bulge,halo,parsb] = galparams
    omega = barra.omega
    xi=ini_guess
    xacc=1e-14
    OMEGA2 = omega*omega
    pequil = []

    for i in range(5):
        xf=[]
        [xf,infodict,exitflag,msg] = fsolve(func2, xi[i][:3], 
                                        args=galparams, 
                                        fprime=func2jac,
                                        full_output = options["verbose"],
                                        xtol = options["tolerance"],
                                        maxfev = options["maxiter"]); #options?
            
        if exitflag!=1:
            print("Lagrangian Point",i+2,"did not converge:",msg)
        
        for j in range(3):
            if (abs(xf[j])<=xacc):
                xf[j]=0
        pequil.append(xf)
        
    #NO CANVIO ORDRE PERQUE ES TORNARÀ BOIG
    peqL1 = pequil[0] # punto equilibrio derecho
    peqL2 = pequil[1] # punto equilibrio izquierdo
    peqL3 = pequil[2] # punto equilibrio central
    peqL4 = pequil[3] # punto equilibrio superior
    peqL5 = pequil[4] # punto equilibrio inferior
    peqs = np.array([peqL1, peqL2, peqL3,peqL4, peqL5])
    return peqs
