from scipy.optimize import fsolve
from .func2 import func2, func2jac

def puntequil(ini_guess,barra,disco,bulge,halo,parsb,options):
    '''
    ini_guess: 5x6 vector containing 5 6D initial guesses, one for each Lagrangian point
    barra: object of class Barra
    options: dictionary of settings for function solver
        verbose: boolean, defafult false
        tolerance: float, default 1e-08
        maxiter: int, default 300
    '''
    omega = barra.omega
    xi=ini_guess
    xf=[]
    xacc=1e-14

    OMEGA2 = omega*omega
    xa=[]
    xb=[]
    pequil = []

    for i in range(5):
        xa.append(xi[i]) #els dos tenen shape 5x3

        '''
        initially implemented with fsolve in matlab
        fsolve takes a function fun, an initial value x0 and some options,
        and solves F(x)=0 iteratively from F(x0)approx=0. is able to use jacobian matrix

        scipy fsolve works in a similar manner
        ''' 

        print(i,"iniguess",xa[i][:3])
        [xf,infodict,exitflag,msg] = fsolve(func2, xa[i][:3], 
                                        args=(barra,disco,bulge,halo,parsb), 
                                        fprime=func2jac,
                                        full_output = options["verbose"],
                                        xtol = options["tolerance"],
                                        maxfev = options["maxiter"]); #options?
        
        print(i,"final",xf)        
        if exitflag!=1:
            print("Lagrangian Point",i,"did not converge:",msg)
        
        for j in range(3):
            if (abs(xf[j])<=xacc):
                xf[j]=0
        pequil.append(xf)
        
    #ORDRE VE DONAR PER INPUT FILE
    peqL1 = [barra.eps,pequil[1],0,0,0] # punto equilibrio derecho
    peqL2 = [barra.eps,pequil[2],0,0,0] # punto equilibrio izquierdo
    peqL3 = [barra.eps,pequil[0],0,0,0] # punto equilibrio central
    peqL4 = [barra.eps,pequil[3],0,0,0] # punto equilibrio superior
    peqL5 = [barra.eps,pequil[4],0,0,0] # punto equilibrio inferior

    return [peqL1, peqL2, peqL3,peqL4, peqL5]
