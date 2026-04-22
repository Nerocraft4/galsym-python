'''
A set of mathematical helper functions 
'''
import numpy as np

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

def elint(PHI,XK):
# C******************************************************************
# C Incomplete elliptic integrals F & E
# C
# C Method : descending Landen transformation
# C Ref.: Abramowitz & Stegun 17.5 p. 597
# C
# C Input parameters : PHI, XK : amplitude angle and module
# C 0 = XK & 1
# C Output values : F, E : incomplete elliptic integrals of the
# C first and second kind
# C*****************************************************************
    PI=np.pi
    PI2 = 2*PI
    A = 1
    C = XK
    S = C*C
    B = np.sqrt(1-S)

    P = PHI
    P2 = P #TODO?
    T = 0
    XI = 1

    while (C>1e-15):
        XI = XI + XI
        P = P + np.floor(P/PI+.5)*PI + np.arctan(np.tan(P2)*B/A)
        P2 = P%PI2
        C = .5*(A-B)
        A1= .5*(A+B)
        B = np.sqrt(A*B)
        A = A1
        S = S + C*C*XI
        T = T + C*np.sin(P2)
#        IF (C.GT.1D-15) GOTO 10 ###GOTO??? this is ancient
    F = P/(XI*A)
    E = T + F*(1-.5*S)
    return [F,E]

def xlmbd(X,Y,Z,A,B,C):
# C*********************************************************************
# C       Resolution of the equation for L:
# C
# C       X/(A+L) + Y/(B+L) + Z/(C+L) = 1
# C
# C       Restriction :X/A + Y/B + Z/C > 1.
# C***********************************************************************
# C The input parameters are assumed to be positive
# C***********************************************************************
    US3=1/3;
    C2 = (A+B+C-X-Y-Z)*US3;
    C1 = A*(B-Y-Z) + B*(C-X-Z) + C*(A-X-Y);
    C0 = A*B*(C-Z) - (A*Y+X*B)*C;
    P = C1*US3-C2*C2;
    Q = C2**3 + (C0-C1*C2)*.5;
    DE = P**3+Q*Q;

    if DE < 0:
        R = np.sqrt(-P);
        XLA = 2*R*np.cos(np.arccos(-Q/R**3)*US3)-C2;
    else:
        R = -Q+np.sqrt(DE);
        if R < 0:
            U = -(-R)**US3;
        else:
            U = R**US3;
        R = -Q-np.sqrt(DE);
        if R < 0:
            V = -(-R)**US3;
        else:
            V = R**US3;
        XLA = U + V - C2;
    return XLA

from models import pot

def isopotencial(radius, idim, barra, disco, bulge, halo, parsb):
    step = 2*radius/(1.0*idim)
    curva = []
    for i in range(idim):
        for j in range(idim):
            x = -radius+step*i
            y = -radius+step*j
            z = 0
            potef = pot.efectivo(x,y,z,barra, disco, bulge, halo, parsb)
            
            curva.append([x,potef])
    return curva



