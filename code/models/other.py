'''
Other partial derivatives and useful functions
'''

import numpy as np
from numpy import sqrt
from maths.helpers import xlmbd 

def derl(barra,x2,y2,z2):
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

def derw(barra,x,y,z,i,j,k):
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

def centro_masas_halo(xdbulge,xydhalo,barra,disco,bulge,halo):
    md = disco.GM
    centrod = np.array([0,0,0])
    md = 0 #zero perque està centrat i el tenim com a referència
    
    mb = barra.GM
    centrob = np.array([0,0,0])

    mesf = bulge.GM
    centroesf = np.array([xdbulge,0,0])

    mh = halo.GM
    centroh = np.array([xydhalo[0], xydhalo[1],0])

    mt = md+mb+mesf+mh
    
    cm = (1/mt)*(md*centrod +mb* centrob + mesf* centroesf + mh * centroh)

    return cm
