import numpy as np
from numpy import sin, cos, sqrt, arcsin
from maths.helpers import xlmbd, elint

def bar(x,y,z,barra,parsb):
    x = x-barra.xd
    y = y-barra.yd

    US3 = 1/3
    US6 = 1/6

    A2 = barra.a*barra.a
    B2 = barra.b*barra.b
    C2 = barra.c*barra.c

    X2=x*x
    Y2=y*y
    Z2=z*z

    if (X2*parsb.UA2+Y2*parsb.UB2+Z2*parsb.UC2<=1):
        potentra = US6*parsb.CTE*(parsb.V000-6*X2*Y2*Z2*parsb.V111 
                       +X2*(X2*(3*parsb.V200-X2*parsb.V300) 
                            +3*(Y2*(parsb.V110+parsb.V110-Y2*parsb.V120
                                -X2*parsb.V210)-parsb.V100))       
                       +Y2*(Y2*(3*parsb.V020-Y2*parsb.V030) 
                            +3*(Z2*(parsb.V011+parsb.V011-Z2*parsb.V012
                                -Y2*parsb.V021)-parsb.V010)) 
                       +Z2*(Z2*(3*parsb.V002-Z2*parsb.V003) 
                            +3*(X2*(parsb.V101+parsb.V101-X2*parsb.V201
                                -Z2*parsb.V102)-parsb.V001)))
        return potentra

    XL = xlmbd(X2,Y2,Z2,A2,B2,C2);

    UA3 = 1./(A2+XL);
    B3 = B2+XL;
    UB3 = 1./B3;
    UC3 = 1./(C2+XL);
    PHI= arcsin(sqrt(UA3/parsb.UA2C2));
    [F,E]=elint(PHI,parsb.XK);
    D2 = 2*sqrt(UA3*UB3*UC3);

    W100 = 2*(F-E)*parsb.UA2B2*parsb.SUA2C2;
    W001 = (D2*B3 - 2*E*parsb.SUA2C2)*parsb.UB2C2;
    W010 = D2 - W100 - W001;

    W110 = (W010 - W100)*parsb.UA2B2;
    W101 = (W001 - W100)*parsb.UA2C2;
    W011 = (W001 - W010)*parsb.UB2C2;

    W200 = (D2*UA3 - W110 - W101)*US3;
    W020 = (D2*UB3 - W110 - W011)*US3;
    W002 = (D2*UC3 - W011 - W101)*US3;

    W111 = (W011 - W110)*parsb.UA2C2;
    W210 = (W110 - W200)*parsb.UA2B2;
    W201 = (W101 - W200)*parsb.UA2C2;
    W120 = (W020 - W110)*parsb.UA2B2;
    W021 = (W011 - W020)*parsb.UB2C2;
    W102 = (W002 - W101)*parsb.UA2C2;
    W012 = (W002 - W011)*parsb.UB2C2;

    W300 = (D2*UA3*UA3 - W210 - W201)*.2;
    W030 = (D2*UB3*UB3 - W120 - W021)*.2;
    W003 = (D2*UC3*UC3 - W102 - W012)*.2;
    
    potnoentra = US6*parsb.CTE*(2*F*parsb.SUA2C2-6*X2*Y2*Z2*W111
                    +X2*(X2*(3*W200-X2*W300)
                        +3*(Y2*(W110+W110-Y2*W120-X2*W210) -W100)) 
                    +Y2*(Y2*(3*W020-Y2*W030)
                        +3*(Z2*(W011+W011-Z2*W012-Y2*W021)-W010)) 
                    +Z2*(Z2*(3*W002-Z2*W003)
                        +3*(X2*(W101+W101-X2*W201-Z2*W102)-W001)))
    return potnoentra

def disk(x,y,z,disco):
    #Miyamoto
    x = x-disco.xd
    y = y-disco.yd
    
    a = disco.a
    b = disco.b #TODO %Estamos en 2dimensiones. Tomamos BM = B^2, como B=1 aqui no hay problema  pero cuidado si B vale distinto de 1!!!!! #TODO EN TEORIA HO HE FIXEAT ABAIX
    GM = disco.GM 
    R2 = x*x+y*y
    Z2 = z*z
    AUX = R2+(a+sqrt(z*z+b*b))*(a+sqrt(z*z+b*b))
    phi = -GM/sqrt(AUX)
    return phi

def bulge(x,y,z,bulge):
    #Plummer
    x = x-bulge.xd
    y = y-bulge.yd
    
    b = bulge.b
    GM = bulge.GM
    rad = sqrt(x*x+y*y+z*z+b*b)
    phi = -GM/rad
    return phi

def halo(x,y,z,halo):
    #Plummer
    x = x-halo.xd
    y = y-halo.yd
    
    b = halo.b
    GM = halo.GM
    rad = sqrt(x*x+y*y+z*z+b*b)
    phi = -GM/rad
    return phi

def efectivo(x,y,z,mbarra,mdisco,mbulge,mhalo,parsb):
    eps = mbarra.eps
    omega = mbarra.omega
    Q1 = sin(eps)
    Q2 = cos(eps)
    
    pbar = bar(x,y,z,mbarra,parsb)
    pdisk = disk(x,y,z,mdisco)
    pbulge = bulge(x,y,z,mbulge)
    phalo = halo(x,y,z,mhalo)

    potmf = pbar + pdisk + pbulge + phalo
    potef = potmf - 0.5*omega**2*((Q2*x)**2+y**2+(Q1*z)**2)-omega**2*Q1*Q2*x*z
    return potef

