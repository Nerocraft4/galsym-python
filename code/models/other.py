import numpy as np
from numpy import sqrt, sin, cos
from maths.helpers import xlmbd
from models import der1, der2
from models import pot

from numpy.typing import NDArray

'''
Other partial derivatives and useful functions
'''

def derFdelta(delta,xvec,barra,disco,bulge,halo,parsb):
    #TODO CURRENTLY UNUSED
    # Derivada de func2 respecte delta, que és el desplacament de la bulge en x!!! #TODO
    #
    #
    epsilon = barra.eps
    omega = barra.omega
    OMEGA2 = omega*omega

    x = xvec[0]
    y = xvec[1]
    z = xvec[2]
    Q1 = sin(epsilon)
    Q2 = cos(epsilon)

    [pxxb,pyyb,pzzb,pxyb,pxzb,pyzb]=der2.bar(barra,parsb,x,y,z,omega)
    [pxxd,pyyd,pzzd,pxyd,pxzd,pyzd]=der2.disk(disco,x,y,z,omega)
    [pxxbl,pyybl,pzzbl,pxybl,pxzbl,pyzbl]=der2.bulge(bulge,x,y,z,omega)
    [pxxh,pyyh,pzzh,pxyh,pxzh,pyzh]=der2.halo(halo,x,y,z,omega)

    #print("pxy",pxyd,pxyb,pxybl,pxyh)

    pxx=pxxd+pxxb+pxxbl+pxxh
    a = OMEGA2*Q2*Q2 - pxx #TODO estic segur d'això? té sentit, però... idk
    b = OMEGA2*Q1*Q2
    return np.array([0,0,0,a,0,b])

def derl(barra,x2: float,y2: float,z2: float) -> float:
    # C*************************************************************************
    # C calcul de les derivades parcials dL/dx, dL/dy, dL/dz
    # C obtingudes mitjanÃ§ant derivada implicita de la funcio
    # C         x^2          y^2           z^2
    # C p(L)=__________ + __________ + ___________ - 1
    # C        A^2+L        B^2+L         C^2+L
    # C
    # C NOTA! LI ENTRO X2,Y2,Z2!
    # C DBP
    # C*************************************************************************

    A2 = barra.a*barra.a
    B2 = barra.b*barra.b
    C2 = barra.c*barra.c

    xl=xlmbd(x2,y2,z2,A2,B2,C2)
    x=sqrt(x2)
    y=sqrt(y2)
    z=sqrt(z2)
    aux1=x/(A2+xl)
    aux2=y/(B2+xl)
    aux3=z/(C2+xl)
    aux=aux1*aux1+aux2*aux2+aux3*aux3
    dx=2*aux1/aux
    dy=2*aux2/aux
    dz=2*aux3/aux
    return [dx,dy,dz]

def derw(barra,x,y,z,i,j,k) -> float:
    # C***************************************************************************
    # C calcul de
    # C                     1                1                 1
    # C dWijk/dl = - _______________ * _______________ * _______________
    # C              (A^2+l)^(i+1/2)   (B^2+l)^(j+1/2)   (C^2+l)^(k+1/2)
    # C
    # C NOTA! LI ENTRO X2,Y2,Z2!
    # C DBP
    # C***************************************************************************

    A2 = barra.a*barra.a
    B2 = barra.b*barra.b
    C2 = barra.c*barra.c

    xl=xlmbd(x,y,z,A2,B2,C2)
    aux1=sqrt(A2+xl)*((A2+xl)**i)
    aux2=sqrt(B2+xl)*((B2+xl)**j)
    aux3=sqrt(C2+xl)*((C2+xl)**k)
    d=-1/(aux1*aux2*aux3)
    return d

def centro_masas_halo(xydbulge,xydhalo,galparams: list) -> list:
    [barra,disco,bulge,halo,parsb] = galparams

    epsilon = barra.eps

    md = disco.GM
    centrod = np.array([0,0,0])#zero perque està centrat i el tenim com a referència

    mb = barra.GM
    centrob = np.array([0,0,0])

    mesf = bulge.GM
    xd,yd,zd = xydbulge #TODO obligar eventualment a que això siguin vectors 3D
    centroesf = np.array([xd*cos(epsilon),yd,xd*sin(epsilon)+zd])

    mh = halo.GM
    xd,yd,zd = xydhalo
    centroh = np.array([xd*cos(epsilon),yd,xd*sin(epsilon)+zd])

    mt = md+mb+mesf+mh
    cm = (1/mt)*(md*centrod +mb* centrob + mesf* centroesf + mh * centroh)
    return cm

def update_displacement(obj: object, disp: list):
    obj.xd = disp[0]
    obj.yd = disp[1]
    obj.zd = disp[2]

#TODO segurament moure a initializer
from utils.io import extract_galparams, pack_galparams
def setup(galparams: dict, displacements: list) -> dict:
    print("Setting up system")
    [barra, disco, bulge, halo, parsb] = extract_galparams(galparams)
    [despbar,despdis,despbul,desphal] = displacements

    print(displacements)
    print(barra.eps)

    [xcm, ycm, zcm] = centro_masas_halo(despbul,desphal,[barra, disco, bulge, halo, parsb])
    print("Initial step: Centro masas halo",xcm,ycm,zcm)

    despbar = [-xcm,-ycm,-zcm]
    despbul = [despbul[0]-xcm,despbul[1]-ycm,despbul[2]-zcm]

    update_displacement(barra,despbar)
    update_displacement(bulge,despbul)
    update_displacement(disco,despdis)
    update_displacement(halo,desphal)
    return galparams

def update(galparams: dict, displacements: list,
           whichobject: str, whichparam: str, p: float) -> None:
    #print("updating parameters and adjusting to new center of masses")
    [barra, disco, bulge, halo, parsb] = extract_galparams(galparams)
    [despbar,despdis,despbul,desphal] = displacements

    print(whichobject,whichparam,p)

    if whichobject == "halo":
        if whichparam == "xd":
            desphal[0] = p
        elif whichparam == "yd":
            desphal[1] = p
        elif whichparam == "zd":
            desphal[2] = p
        else:
            print("Continuation for object",whichobject,"param",whichparam,"not yet implemented")
            return
    elif whichobject == "bulge":
        if whichparam == "xd":
            despbul[0] = p
        elif whichparam == "yd":
            despbul[1] = p
        elif whichparam == "zd":
            despbul[2] = p
        else:
            print("Continuation for object",whichobject,"param",whichparam,"not yet implemented")
            return
    elif whichobject == "barra":
        if whichparam == "eps":
            barra.eps = p
        else:
            print("Continuation for object",whichobject,"param",whichparam,"not yet implemented")
            return
    else:
        print("Continuation for object",whichobject,"param",whichparam,"not yet implemented")
        return

    [xcm, ycm, zcm] = centro_masas_halo(despbul,desphal,[barra, disco, bulge, halo, parsb])
    print("Centro masas halo",xcm,ycm,zcm)

    despbar = [-xcm,-ycm,-zcm]
    despbul = [despbul[0]-xcm,despbul[1]-ycm,despbul[2]-zcm]

    update_displacement(barra,despbar)
    update_displacement(bulge,despbul)
    update_displacement(disco,despdis)
    update_displacement(halo,desphal)
    return

def CTJAC(xvec: NDArray, pvec: NDArray, params: list) -> float:
    '''
    Rutina que calcula la constant de Jacobi pel cas d'un potencial de
    barra.
    El punt X ha d'estar en coordenades sinodiques (x,y,z,xd,yd,zd)
    La constant de Jacobi ve definida com,
    CTE=1/2*(xd^2+yd^2+zd^2)+phi_e(x,y,z) on phi_e es el pot efectiu
    '''
    [barra,disco,bulge,halo,parsb] = params
    x,y,z = xvec
    xp,yp,zp = pvec

    EPS=barra.eps
    OMEGA = barra.omega
    Q1=sin(EPS)
    Q2=cos(EPS)
    POT = pot.efectivo(x,y,z,params)
    CTJ=POT+0.5*(xp*xp+yp*yp+zp*zp)
    return CTJ

def DCTJAC(xvec: NDArray, pvec: NDArray,params: list) -> float:
    '''
    Rutina que calcula el gradient de la CJAC en un punt
    El punt X ha d'estar en coordenades sinodiques (x,y,z,xd,yd,zd)
    '''
    [barra,disco,bulge,halo,parsb] = params
    x,y,z = xvec
    xp,yp,zp = pvec

    EPS=barra.eps
    omega = barra.omega
    Q1=sin(EPS)
    Q2=cos(EPS)
    POT = der1.efectivo(x,y,z,params)
    dCJAC = [0,0,0,0,0,0]
    dCJAC[0] = POT[0]
    dCJAC[1] = POT[1]
    dCJAC[2] = POT[2]
    dCJAC[3] = xp
    dCJAC[4] = yp
    dCJAC[5] = zp
    return dCJAC
