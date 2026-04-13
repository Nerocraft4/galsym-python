'''
A set of mathematical helper functions 
'''
import numpy as np
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
