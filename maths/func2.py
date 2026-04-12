from ..models.derivatives import *
import numpy as np

def func2(xvec,barra,parsb):
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

    [pxd,pyd,pzd]=der1disk(barra,parsb,x,y,z,omega)
    [pxb,pyb,pzb]=der1bar(barra,parsb,x,y,z,omega)
    [pxh,pyh,pzh]=der1halo(barra,parsb,x,y,z,omega)
    [pxbl,pybl,pzbl]=der1bulge(barra,parsb,x,y,z,omega)
    px=pxd+pxb+pxbl+pxh
    py=pyd+pyb+pybl+pyh
    pz=pzd+pzb+pzbl+pzh

    F[0]=OMEGA2*q1*q2*z+OMEGA2*q2*q2*x-px
    F[1]=OMEGA2*y-py
    F[2]=OMEGA2*q1*q2*x+OMEGA2*q1*q1*z-pz
    #print("F",F)
    return F

def func2jac(xvec,barra,parsb):
    epsilon = barra.eps
    omega = barra.omega

    x = xvec[0]
    y = xvec[1]
    z = xvec[2]
    DF = [[0,0,0],[0,0,0],[0,0,0]]
    q1=np.sin(epsilon);
    q2=np.cos(epsilon);
    OMEGA2=omega*omega;

    [pxxd,pyyd,pzzd,pxyd,pxzd,pyzd]=der2disk(barra,parsb,x,y,z,omega);
    [pxxb,pyyb,pzzb,pxyb,pxzb,pyzb]=der2bar(barra,parsb,x,y,z,omega);
    [pxxbl,pyybl,pzzbl,pxybl,pxzbl,pyzbl]=der2bulge(barra,parsb,x,y,z,omega);
    [phxx,phyy,phzz,phxy,phxz,phyz]=der2halo(barra,parsb,x,y,z,omega);

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
    
    #print("DF",DF)
    return DF






