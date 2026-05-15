from models.other import CTJAC
from numpy import floor, log10
import copy

def refinamiento(pequil,paprox,times,params,CAMP,CJAC,GRADC,SECCIO,GRADS):
    '''
    Funció que gestiona el refinament de les òrbites de Lyapunov (crec?)
    Input:
        pequil: punt d'equilibri sobre el que s'han centrat les orbites approx
        paprox: llista de punts 6D al voltant de pequil
        times: temps de 0 a 2pi (entenc) corresponents als punts de paprox
    Informació legacy:
    % C**************************************************************************
    % C Rutines pel calcul d'orbites periodiques incloent tir paral.lel
    % C REFOP          Rutina principal del calcul de OP
    % C ITEROP         Auxiliar de REFOP. Fa una iteracio de calcul OP.
    % C OMPLEDF        Rutina auxiliar per omplir la diferencial del sist.
    % C**************************************************************************
    '''
    [barra,disco,bulge,halo,parsb] = params
    dc = 1e-3# valors entre [1e-3, 1e-6]
    # xkk=2e-3, dc=5e-5, ydhalo=-3, xdbulge=2.5 funciona
    n=6
    prec= 1e-10#; % 1e-5; %1e-10;
    prefin = [];
    eps = barra.eps #redundant

    i=0;
    
    xa = [pequil[0], pequil[1], pequil[2], 0,0,0]
    cjacorig = CTJAC(pequil,[0,0,0],params)

    xkkorig = copy.deepcopy(paprox) 
    cjacs = [cjacorig]

    perorig = times[-1] 
    np = i
    mode = 2 #0,1,2 depenent de la informació o equacions que tinguem. 2=fixar energia
    imax = 20#; %50; %20;
    j = 0 #TODO revisar flags
    exito = 0

    while exito == 0 and j<1:
        indn = 0
        xkk = xkkorig
        per = perorig
        cjac = cjacorig+1*dc
        cjacs.append(cjac)
        while j < 1:
            [xkk,times2,indn,per]=refop(xkk,times,cjac,imax,n,np,per,prec,CAMP,CJAC,GRADC,SECCIO,GRADS,mode);
            if indn==1:
                exito = 0
                break #% Salimos del while, para volver a repetir el bucle con otros parametros
            # Guardamos: eps, c, per, xkk(1:6). Es el archivo orbitesL1 de fortran
            prefin.extend(eps, cjacs[1], per, xkk[0:5]) #TODO extrany

            ppinta = pintaorbita(xkk, times2, per, n, np, 100, CAMP) # ppinta contiene: times, puntos orbita
            
            j = j+1 #TODO potencialment redundant?
            xkk[n*np]+=1*dc
            exito = 1
            print('refinamiento: exito ', exito)
        print("dc",dc)
        dc += 10**floor(log10(dc))
        print('refinamiento: dc ',dc)
    
    return c,prefin,ppinta


def refop(xkk,temps,cjac,imax,n,np,per,prec,CAMP,CJAC,GRADC,SECCIO,GRADS,mode):
    '''
    % C***********************************************************************
    % C REFINAMENT D'ORBITA PERIODICA: 
    % C INPUT:
    % C XKK: punts q fem servir com a llavor inicial
    % C      es un vector de dimensio (N x NP)
    % C      XKK=[x_1,x_2,..,x_NP] on xi es el vector de dim N amb p+v.
    % C cjac: energia
    % C TEMPS: TEMPS ASSOCIAT A CADA UN DELS dels estats x_i
    % C IMAX: NUMERO MAXIM D'ITERACIONS PERMESES EN EL NEWTON
    % C N: DIMENSIO DE L'ESPAI (POS+VEL)
    % C NP:NUMERO DE PUNTS SOBRE L'ORBITA PERIODICA #TODO REDUNDANT CREC
    % C PREC: PRECISIO DEMANADA com a criteri de parada
    % C CAMP: CAMP VECTORIAL QUE INTEGREM
    % C CJAC: INTEGRAL PRIMERA (CONSTANT DE JACOBI) #TODO REDUNDANT?
    % C GRADC: GRADIENT DE LA INTEGRAL PRIMERA (CONSTANT DE JACOBI)
    % C SECCIO: EQUACIO DE  LA SECCIO
    % C GRADS: GRADIENT DE LA SECCIO
    % C mode: ENS INDICA QUE FEM: 0: NO FIXEM ENERGIA NI PERIODE
    % C                           1: FIXEM PERIODE (ELIMINEM EQ. ENERGIA)
    % C                           2: FIXEM ENERGIA (ELIMINEM EQ. PERIODE)
    % C OUTPUT:
    % C XKK=prefop: PUNTS REFINATS SOBRE LA NOVA ORBITA PERIODICA
    % C TEMPS=temps2: VECTOR DE TEMPS ASSOCIATS ALS NOUS PUNTS x_i.
    % C INDN: INDICADOR: 0: TOT BE, NEWTON HA CONVERGIT
    % C                  1: NEWTON NO HA CONVERGIT
    % C****
    '''
    print("refop")
    ite = 0
    xnorm = 1e10 #NORMA2 de correc AUX (no te en compte energia/periode)
    xnorm2 = 1e10 #NORMA 2 de funció F
    ind = 0 #flag diferent a la de newton, no entenc perquè
    flagnewton = 0 #suposo que flag de newton, antic indn
    
    while xnorm2>prec and ite<imax and ind==0:    
        [piterop,tempsit,xnorms,xnorm2s,ind,per2] = iterop(xkk,temps,n,np,per,CAMP,
                                                           CJAC,GRADC,SECCIO,GRADS,
                                                           xnorm,xnorm2,mode)
        xkk = piterop;
        temps = tempsit;
        xnorm = xnorms;
        xnorm2 = xnorm2s;
        per = per2;
        ite =ite+1;
    if ite == imax:
        print("Newton no ha convergit, iteracions",str(ite))
        flagnewton = 1
    if ind == 1:
        print("Error amb Newton a iterop")
        flagnewton = 1
    if xnorm2 < prec:
        print("S'ha convergit amb ", str(ite),"iteracions")
    
    prefop = xkk #TODO gestionar
    return [prefop,temps,flagnewton,per2]


def iterop(xkk,temps,n,np,per,CAMP,CJAC,GRADC,SECCIO,GRADS,xnorm,xnorm2,mode):
    '''
    % C***********************************************************************
    % C UNA ITERACIO EN EL REFINAMENT D'UNA OP. RESOLUCIO DEL SISTEMA DF*aux=F 
    % C ON
    % C DF: MATRIU DIFERENCIAL DEL SISTEMA I TE DIMENSIO NE*(N*NP)
    % C     NE=N*NP+1 SI FIXEM ENERGIA O PERIODE I NE=N*NP SI NO FIXEM RES.
    % C     NP=NUMERO DE PUNTS SOBRE L'ORBITA PERIODICA
    % C     N=DIMENSIO DE L'ESPAI
    % C AUX= VECTOR DE CORRECCIONS
    % C F= VECTOR DE LA IMATGE DELS PUNTS PEL CAMP
    % C INPUT:
    % C XKK: punts q fem servir com a llavor inicial
    % C      es un vector de dimensio (N x NP + 1) si fixem energia o periode
    % C      o (N x NP) si no fixem res i fem norma minima de correccio.
    % C      XKK=[x_1,x_2,..,x_NP, h o tau] on xi es el vector de dim N 
    % C      amb p+v.
    % C TEMPS: TEMPS ASSOCIAT A CADA UN DELS PUNTS x_i
    % C N: DIMENSIO DE L'ESPAI (POS+VEL)
    % C NP:NUMERO DE PUNTS SOBRE L'ORBITA PERIODICA
    % C CAMP: CAMP VECTORIAL QUE INTEGREM
    % C CJAC: CONSTANT DE JACOBI
    % C GRADC: GRADIENT DE LA CONSTANT DE JACOBI
    % C SECCIO: EQUACIO DE  LA SECCIO
    % C GRADS: GRADIENT DE LA SECCIO
    % C IFEM: ENS INDICA QUE FEM: 0: NO FIXEM ENERGIA NI PERIODE
    % C                           1: FIXEM PERIODE (ELIMINEM EQ. ENERGIA)
    % C                           2: FIXEM ENERGIA (ELIMINEM EQ. PERIODE)
    % C OUTPUT:
    % C XKK: PUNTS REFINATS SOBRE LA NOVA ORBITA PERIODICA
    % C TEMPS: VECTOR DE TEMPS ASSOCIATS ALS NOUS PUNTS xi.
    % C XNORM: NORMA DEL VECTOR AUX
    % C XNORM2: NORMA DEL VECTOR F
    % C IND: INDICADOR: 1 -> MATRIU DE SISTEMA SINGULAR, O BE ALTRE ERROR
    % C                      D'ENTRADA (per ex IFEM erroni).
    % C                 0 -> MATRIU DE SISTEMA NO SINGULAR
    % C***********************************************************************
    '''
    piterop = 0
    tempsit = 0
    xnorms = 0
    xnorm2s = 0
    ind = 0
    per2 = 0    
    return [piterop,tempsit,xnorms,xnorm2s,ind,per2]















