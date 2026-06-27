from scipy.optimize import fsolve
from .func2 import func2, func2jac
import numpy as np
from numpy.typing import NDArray

def puntequil(ini_guessess: list, galparams: list, options: dict) -> NDArray:
    '''
    This function fins an equilibrium point in the field of the galaxy given galactic
    parameters and an initial guess, through Newton's method (fsolve)
    Input:
        ini_guess:  Nx3 vector containing N 3D initial guesses, one for each Lagr point
                    velocities are assumed to be zero. N usually is 5
        galparams:  list of objects containing galactic params
        options:    dictionary of settings for function solver
            verbose:    boolean, defafult false
            tolerance:  float, default 1e-08
            maxiter:    int, default 300
    Output:
        peqs:   5x3 array containing the equilibrium points of the system
    '''
    [barra,disco,bulge,halo,parsb] = galparams
    N_guesses = len(ini_guessess)
    omega = barra.omega
    xacc=1e-14
    OMEGA2 = omega*omega
    pequil = []

    for i in range(N_guesses):
        xf=[]
        [xf,infodict,exitflag,msg] = fsolve(func2, ini_guessess[i][:3], 
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
        
    peqs = np.array(pequil)
    return peqs
