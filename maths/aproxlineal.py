
global barra#, disco, bulge, parsb, ctes_rk78 #TODO revisar això de les globals

def aproxlineal(peq,xkk):
    '''
    peq: 1x3 array containing x,y,z coordinates of an equilibrium point
    xkk: #TODO
    '''
    omega = barra["omega"];
    omega2 = omega*omega;
    eps = barra["eps"]

    hmin = 1e-3
    hmax = 1
    tolrk = 1e-13
    
    t=0 #TODO maybe not necessary
    N=6 #it can be N=42 if using variational. maybe inclue this as a param?

    x,y,z = peq #TODO puc fer això directament?
    print(x,y,z)

    [pxxb,pyyb,pzzb,pxyb,pxzb,pyzb]=der2bar(x,y,z,omega);
    [pxxm,pyym,pzzm,pxym,pxzm,pyzm]=der2miya(x,y,z,omega);
    [pxxbl,pyybl,pzzbl,pxybl,pxzbl,pyzbl]=der2bulge(x,y,z,omega);
    [phxx,phyy,phzz,phxy,phxz,phyz]=der2halo(x,y,z,omega);
    
barra = {"omega":0.3,"eps":0.01} #potser podria fer objectes, no sé si val la pena

aproxlineal(peq=[0,0,1],xkk=1)
