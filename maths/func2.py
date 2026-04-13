from ..models import der1
from ..models import der2
import numpy as np

def func2(xvec,barra,disco,bulge,halo,parsb):
    '''
    #TODO encara he de veure què fa exactament aquesta funció
    sembla que va de R^3 a R^3xR^6
    és a dir, envia xvec(3dim) a F(3) i DF(6).
    '''
    epsilon = barra.eps
    omega = barra.omega

    x = xvec[0]
    y = xvec[1]
    z = xvec[2]
    F = [0,0,0]

    q1=np.sin(epsilon);
    q2=np.cos(epsilon);
    OMEGA2=omega*omega;

    [pxb,pyb,pzb]=der1.bar(barra,parsb,x,y,z,omega)
    [pxd,pyd,pzd]=der1.disk(disco,x,y,z,omega)
    [pxbl,pybl,pzbl]=der1.bulge(bulge,x,y,z,omega)
    [pxh,pyh,pzh]=der1.halo(halo,x,y,z,omega)    
    px=pxd+pxb+pxbl+pxh
    py=pyd+pyb+pybl+pyh
    pz=pzd+pzb+pzbl+pzh

    F[0]=OMEGA2*q1*q2*z+OMEGA2*q2*q2*x-px
    F[1]=OMEGA2*y-py
    F[2]=OMEGA2*q1*q2*x+OMEGA2*q1*q1*z-pz
    return F

def func2jac(xvec,barra,disco,bulge,halo,parsb): #TODO potser és inferit i no cal?
    epsilon = barra.eps
    omega = barra.omega

    x = xvec[0]
    y = xvec[1]
    z = xvec[2]
    DF = [[0,0,0],[0,0,0],[0,0,0]]
    q1=np.sin(epsilon);
    q2=np.cos(epsilon);
    OMEGA2=omega*omega;

    [pxxb,pyyb,pzzb,pxyb,pxzb,pyzb]=der2.bar(barra,parsb,x,y,z,omega);
    [pxxd,pyyd,pzzd,pxyd,pxzd,pyzd]=der2.disk(disco,x,y,z,omega);
    [pxxbl,pyybl,pzzbl,pxybl,pxzbl,pyzbl]=der2.bulge(bulge,x,y,z,omega);
    [phxx,phyy,phzz,phxy,phxz,phyz]=der2.halo(halo,x,y,z,omega);

    pxx=pxxd+pxxb+pxxbl+phxx;
    pyy=pyyd+pyyb+pyybl+phyy;
    pzz=pzzd+pzzb+pzzbl+phzz;
    pxy=pxyd+pxyb+pxybl+phxy;
    pxz=pxzd+pxzb+pxzbl+phxz;
    pyz=pyzd+pyzb+pyzbl+phyz;

    DF[0][0]=OMEGA2*q2*q2-pxx;
    DF[0][1]=-pxy;
    DF[0][2]=OMEGA2*q1*q2-pxz;
    DF[1][1]=OMEGA2-pyy;
    DF[1][2]=-pyz;
    DF[2][2]=OMEGA2*q1*q1-pzz;
    #adding this just to test it out
    DF[2][0]=OMEGA2*q1*q2-pxz;
    DF[1][0]=-pxy;
    DF[2][1]=-pyz;
    return DF






