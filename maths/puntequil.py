from scipy.optimize import fsolve

def puntequil(barra,options,arxi):
    '''
    barra: object of class Barra
    options: dictionary of settings for function solver
        verbose: boolean, defafult false
        tolerance: float, default 1e-08
        maxiter: int, default 300
    '''
    omega = barra.omega
    xi=load(arxi2) #TODO
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
        and solves F(x)=0 iteratively from F(x0)approx=0

        scipy fsolve works in a similar manner
        ''' 
        #but how do we add a parameter?
        [xf,fval,exitflag,output] = fsolve( func2(omega), xa[i], args=(barra), 
                                            full_output = options["verbose"],
                                            xtol = options["tolerance"],
                                            maxfev = options["maxitere"]); #options?
        if exitflag==0:
            break
        
        for j in range(3):
            if (abs(xf[j])<=xacc):
                xf[j]=0
        pequil.append(xf);

    peqL1 = [barra.eps,pequil[1],0,0,0] # punto equilibrio derecho
    peqL2 = [barra.eps,pequil[2],0,0,0] # punto equilibrio izquierdo
    peqL3 = [barra.eps,pequil[0],0,0,0] # punto equilibrio central
    peqL4 = [barra.eps,pequil[3],0,0,0] # punto equilibrio superior
    peqL5 = [barra.eps,pequil[4],0,0,0] # punto equilibrio inferior

    #TODO perquè canviem ordre si fora es torna a desfer???
    return [peqL1, peqL2, peqL3,peqL4, peqL5]
