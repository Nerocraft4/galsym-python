from models import der2
from .DF import DF
from numpy import sin, cos, real, imag, pi, zeros, array
from numpy.linalg import eig
import copy

def aproxlineal(xvec,xkk,params):
    '''
    peq: 1x3 array containing x,y,z coordinates of an equilibrium point
    xkk: default xkk=3e-2
    params: list with all model parameters
    '''
    [barra, disco, bulge, halo, parsb] = params

    hmin = 1e-3
    hmax = 1
    tolrk = 1e-13
    
    t=0 #TODO maybe not necessary
    N=6 #it can be N=42 if using variational. maybe inclue this as a param?
    
    A = DF(xvec,barra,disco,bulge,halo,parsb)
    vaps, veps = eig(A)
    rr = real(vaps)
    ri = imag(vaps)

    rip = []
    for ri_i in ri:
        if ri_i>0:
            rip.append(ri_i)
    xw = min(rip)
    xw2 = xw**2
    xw4 = xw2**2
    
    ti = 0
    tf = 2*pi/xw

    Q3 = -(-A[5][0]**2+A[3][0]*(xw2+A[5][2])+xw4+A[5][2]*xw2)
    Q3 /= xw*(-A[5][4]*A[5][0]+A[3][4]*(A[5][2]+xw2))
    print(Q3)
    S3 = -(-A[5][4]*A[3][0]+A[3][4]*A[5][0]-A[5][4]*xw2)
    S3 /= (-A[5][4]*A[5][0]+A[3][4]*A[5][2]+A[3][4]*xw2)
    
    paprox = [[xvec[0],xvec[1],xvec[2],0,0,0]] #programa original posa "epsilon" com a primer valor de paprox, meh

    m=24 #?

    h = (tf-ti)/m
    t = ti
    x = zeros([6])
    a,b,c = xvec
    ts = []
    for i in range(m):
        x[0] = a + xkk*cos(xw*t)
        x[1] = b + Q3*xkk*sin(xw*t)
        x[2] = c + S3*xkk*cos(xw*t)
        x[3] = -xw*xkk*sin(xw*t)
        x[4] = xw*Q3*xkk*cos(xw*t)
        x[5] = -xw*S3*xkk*sin(xw*t) 
        '''
        if (i == m+1)
            for j=1:6
                xm(j) = -1
            end
            paprox = [paprox -1, xm] #i això simplement posa -1 al final de tot, zzz
        '''
        ts.append(t)
        paprox.append(copy.deepcopy(x))
        t = t+h
    print(paprox)
    paprox = array(paprox)[1:].transpose()
    import matplotlib.pyplot as plt
    plt.scatter(paprox[0],paprox[1])
    plt.show()

