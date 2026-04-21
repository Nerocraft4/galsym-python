from models import der1
from models import der2
import numpy as np
from numpy import sin, cos

def DF(xvec,barra,disco,bulge,halo,parsb):
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
    [phxx,phyy,phzz,phxy,phxz,phyz]=der2.halo(halo,x,y,z,omega)

    pxx=pxxd+pxxb+pxxbl+phxx
    pyy=pyyd+pyyb+pyybl+phyy
    pzz=pzzd+pzzb+pzzbl+phzz
    pxy=pxyd+pxyb+pxybl+phxy
    pxz=pxzd+pxzb+pxzbl+phxz
    pyz=pyzd+pyzb+pyzbl+phyz

    a = OMEGA2*Q2*Q2+pxy
    b = OMEGA2*Q1*Q2+pxz
    c = 2*omega*Q2
    d = OMEGA2+pyz
    e = 2*omega*Q1
    f = OMEGA2*Q1*Q1+pzz 

    df =[[0,  0,  0,  1,  0,  0],
         [0,  0,  0,  0,  1,  0],
         [0,  0,  0,  0,  0,  1],
         [a,  pxy,b,  0,  c,  0],
         [pxy,d,  pyz,-c, 0, -e],
         [b,  pyz,f,  0,  e,  0]]

    return(np.array(df))
