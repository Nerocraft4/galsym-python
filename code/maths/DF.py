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
    [pxxh,pyyh,pzzh,pxyh,pxzh,pyzh]=der2.halo(halo,x,y,z,omega)

    #print("pxy",pxyd,pxyb,pxybl,pxyh)

    pxx=pxxd+pxxb+pxxbl+pxxh
    pyy=pyyd+pyyb+pyybl+pyyh
    pzz=pzzd+pzzb+pzzbl+pzzh
    pxy=pxyd+pxyb+pxybl+pxyh
    pxz=pxzd+pxzb+pxzbl+pxzh
    pyz=pyzd+pyzb+pyzbl+pyzh

    '''
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

    a = OMEGA2*Q2*Q2-pxy
    b = OMEGA2*Q1*Q2-pxz
    c = 2*omega*Q2
    d = OMEGA2-pyz
    e = 2*omega*Q1
    f = OMEGA2*Q1*Q1-pzz 

    df =[[0,    0,    0,    1,  0,  0],
         [0,    0,    0,    0,  1,  0],
         [0,    0,    0,    0,  0,  1],
         [a,    -pxy, b,    0,  c,  0],
         [-pxy, d,   -pyz,  -c, 0, -e],
         [b,    -pyz, f,    0,  e,  0]]

    print("peq",xvec)
    '''
    A = np.zeros([6,6])
    A[0][3]=1
    A[1][4]=1
    A[2][5]=1
    A[3][0]=OMEGA2*Q2*Q2-(pxxb+pxxd+pxxbl+pxxh)
    A[3][1]=-(pxyb+pxyd+pxybl+pxyh)             
    A[3][2]=OMEGA2*Q1*Q2-(pxzb+pxzd+pxzbl+pxzh)
    A[3][4]=2*omega*Q2
    A[4][0]=A[3][1]                             
    A[4][1]=OMEGA2-(pyyb+pyyd+pyybl+pyyh)
    A[4][2]=-(pyzb+pyzd+pyzbl+pyzh)
    A[4][3]=-A[3][4]
    A[4][5]=-2*omega*Q1
    A[5][0]=A[3][2]
    A[5][1]=A[4][2]
    A[5][2]=OMEGA2*Q1*Q1-(pzzb+pzzd+pzzbl+pzzh)
    A[5][4]=-A[4][5]
    #with np.printoptions(precision=4):
    #    print(A)

    return(np.array(A))
