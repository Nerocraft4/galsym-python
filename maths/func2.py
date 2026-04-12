from .models.derivatives import *

def func2(xvec,barra):
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
    DF = [[0,0,0],[0,0,0],[0,0,0]]

    q1=sin(eps);
    q2=cos(eps);
    OMEGA2=omega*omega;

    [pxm,pym,pzm]=der1miya(x,y,z,omega)
    [pxb,pyb,pzb]=der1bar(x,y,z,omega)
    [pxh,pyh,pzh]=der1halo(x,y,z,omega)
    [pxbl,pybl,pzbl]=der1bulge(x,y,z,omega)
    px=pxm+pxb+pxbl+pxh
    py=pym+pyb+pybl+pyh
    pz=pzm+pzb+pzbl+pzh

    F[0]=OMEGA2*q1*q2*z+OMEGA2*q2*q2*x-px
    F[1]=OMEGA2*y-py
    F[2]=OMEGA2*q1*q2*x+OMEGA2*q1*q1*z-pz

    [pxxm,pyym,pzzm,pxym,pxzm,pyzm]=der2miya(x,y,z,omega);
    [pxxb,pyyb,pzzb,pxyb,pxzb,pyzb]=der2bar(x,y,z,omega);
    [pxxbl,pyybl,pzzbl,pxybl,pxzbl,pyzbl]=der2bulge(x,y,z,omega);
    [phxx,phyy,phzz,phxy,phxz,phyz]=der2halo(x,y,z,omega);

    pxx=pxxm+pxxb+pxxbl+phxx;
    pyy=pyym+pyyb+pyybl+phyy;
    pzz=pzzm+pzzb+pzzbl+phzz;
    pxy=pxym+pxyb+pxybl+phxy;
    pxz=pxzm+pxzb+pxzbl+phxz;
    pyz=pyzm+pyzb+pyzbl+phyz;

    DF[0][0]=OMEGA2*q2*q2-pxx;
    DF[0][1]=-pxy;
    DF[0][2]=OMEGA2*q1*q2-pxz;
    DF[1][1]=OMEGA2-pyy;
    DF[1][2]=-pyz;
    DF[2][2]=OMEGA2*q1*q1-pzz;

    return np.array([F,DF])
