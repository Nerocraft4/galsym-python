from . import der1, der2
import numpy as np
def CAMPSB(xvec: list, params: list):
    '''
    Rutina que calcula el camp vectorial del potencial de barra
    (Miyamoto+Ferrersn2): Pot=POTM+POTB
    Les coordenades son sinodiques Q=(x,y,z,xd,yd,zd) i si N=42
    s'hi afegeixen les variacionals
    '''
    #N=6
    [barra, disco, bulge, halo, parsb] = params
    OMEGA = barra.omega
    OMEGA2 = OMEGA*OMEGA
    EPS = barra.eps

    x = xvec[0]
    y = xvec[1]
    z = xvec[2]
    xp = xvec[3]
    yp = xvec[4]
    zp = xvec[5]

    [PBX,PBY,PBZ]=der1.bar(barra,parsb,x,y,z,OMEGA)
    [PMX,PMY,PMZ]=der1.disk(disco,x,y,z,OMEGA)
    [PBXL,PBYL,PBZL]=der1.bulge(bulge,x,y,z,OMEGA)
    [PHX,PHY,PHZ]=der1.halo(halo,x,y,z,OMEGA)

    PX=PBX+PMX+PBXL+PHX
    PY=PBY+PMY+PBYL+PHY
    PZ=PBZ+PMZ+PBZL+PHZ

    Q1=np.sin(EPS)
    Q2=np.cos(EPS)

    F1=xp
    F2=yp
    F3=zp
    F4=-PX+2*OMEGA*Q2*yp+OMEGA2*Q2*Q2*x+OMEGA2*Q1*Q2*z
    F5=-PY-2*OMEGA*Q2*xp-2*OMEGA*Q1*zp+OMEGA2*y
    F6=-PZ+2*OMEGA*Q1*yp+OMEGA2*Q1*Q2*x+OMEGA2*Q1*Q1*z

    F = [F1,F2,F3,F4,F5,F6] # %Cada columna de F es el campo en un punto
    return np.array(F)

def CAMPSB_var(xvec: list, params: list):
    '''
    Rutina que calcula el camp vectorial del potencial de barra
    (Miyamoto+Ferrersn2): Pot=POTM+POTB
    Les coordenades son sinodiques Q=(x,y,z,xd,yd,zd) i
    s'hi afegeixen les variacionals
    '''
    #N=6
    [barra, disco, bulge, halo,parsb] = params
    OMEGA = barra.omega
    OMEGA2 = OMEGA*OMEGA
    EPS = barra.eps

    x = xvec[0]
    y = xvec[1]
    z = xvec[2]
    xp = xvec[3]
    yp = xvec[4]
    zp = xvec[5]

    [PBX,PBY,PBZ]=der1.bar(barra,parsb,x,y,z,OMEGA)
    [PMX,PMY,PMZ]=der1.disk(disco,x,y,z,OMEGA)
    [PBXL,PBYL,PBZL]=der1.bulge(bulge,x,y,z,OMEGA)
    [PHX,PHY,PHZ]=der1.halo(halo,x,y,z,OMEGA)

    PX=PBX+PMX+PBXL+PHX
    PY=PBY+PMY+PBYL+PHY
    PZ=PBZ+PMZ+PBZL+PHZ

    Q1=np.sin(EPS)
    Q2=np.cos(EPS)

    F1=xp
    F2=yp
    F3=zp
    F4=-PX+2*OMEGA*Q2*yp+OMEGA2*Q2*Q2*x+OMEGA2*Q1*Q2*z
    F5=-PY-2*OMEGA*Q2*xp-2*OMEGA*Q1*zp+OMEGA2*y
    F6=-PZ+2*OMEGA*Q1*yp+OMEGA2*Q1*Q2*x+OMEGA2*Q1*Q1*z

    F = [F1,F2,F3,F4,F5,F6] # %Cada columna de F es el campo en un punto

    [pbxx,pbyy,pbzz,pbxy,pbxz,pbyz]=der2.bar(barra,parsb,x,y,z,OMEGA)
    [pmxx,pmyy,pmzz,pmxy,pmxz,pmyz]=der2.disk(disco,x,y,z,OMEGA)
    [pbxxl,pbyyl,pbzzl,pbxyl,pbxzl,pbyzl]=der2.bulge(bulge,x,y,z,OMEGA)
    [phxx,phyy,phzz,phxy,phxz,phyz]=der2.halo(halo,x,y,z,OMEGA)
    PXX = pbxx + pmxx + pbxxl + phxx
    PYY = pbyy + pmyy + pbyyl + phyy
    PZZ = pbzz + pmzz + pbzzl + phzz
    PXY = pbxy + pmxy + pbxyl + phxy
    PXZ = pbxz + pmxz + pbxzl + phxz
    PYZ = pbyz + pmyz + pbyzl + phyz
    for j in range(6):
        J6 = 6*j
        Fn = [0,0,0,0,0,0]
        Fn[0]=xvec[J6+3] #TODO out of bounds!!
        Fn[1]=xvec[J6+4]
        Fn[2]=xvec[J6+5]
        Fn[3]=(OMEGA2*Q2*Q2-PXX)*xvec[J6+0]-PXY*xvec[J6+1]+(OMEGA2*Q1*Q2-PXZ)*xvec[J6+2]+2*OMEGA*Q2*xvec[J6+4]
        Fn[4]=-PXY*xvec[J6+0]+(OMEGA2-PYY)*xvec[J6+1]-PYZ*xvec[J6+2]-2*OMEGA*Q2*xvec[J6+3]-2*OMEGA*Q1*xvec[J6+5]
        Fn[5]=(OMEGA2*Q1*Q2-PXZ)*xvec[J6+0]-PYZ*xvec[J6+1]+(OMEGA2*Q1*Q1-PZZ)*xvec[J6+2]+2*OMEGA*Q1*xvec[J6+4]
        F.append(Fn)
    return np.array(F)
