from ..maths.helpers import elint, xlmbd
from .other import derl, derw
import numpy as np

def bar(barra,parsb,x,y,z,omega):
    #Ferrers
    x = x-barra.xd
    y = y-barra.yd
    A2 = barra.a*barra.a
    B2 = barra.b*barra.b
    C2 = barra.c*barra.c
    x2 = x*x
    y2 = y*y
    z2 = z*z
    US3=1/3
    
    if x2*parsb.UA2+y2*parsb.UB2+z2*parsb.UC2 <= 1:
        pxxb = (-parsb.CTE * (parsb.V100 
                + x2*(5*x2*parsb.V300 + 6*(y2*parsb.V210 - parsb.V200))
                + y2*(y2*parsb.V120 + 2*(z2*parsb.V111 - parsb.V110)) 
                + z2*(z2*parsb.V102 + 6*x2*parsb.V201 - 2*parsb.V101)))
        pyyb = (-parsb.CTE * (parsb.V010 
                + x2*(x2*parsb.V210 + 6*y2*parsb.V120 - 2*parsb.V110)
                + y2*(5*y2*parsb.V030 + 6*z2*parsb.V021 - 6*parsb.V020)
                + z2*(2*x2*parsb.V111 + z2*parsb.V012 - 2*parsb.V011)))
        pzzb = (-parsb.CTE * (parsb.V001 
                + x2*(x2*parsb.V201 + 2*y2*parsb.V111 - 2*parsb.V101)
                + y2*(y2*parsb.V021 + 6*z2*parsb.V012 - 2*parsb.V011)
                + z2*(5*z2*parsb.V003 + 6*x2*parsb.V102 - 6*parsb.V002)))
        pxyb = (-parsb.CTE*4*x*y*(x2*parsb.V210 
                + y2*parsb.V120 + z2*parsb.V111 - parsb.V110))
        pxzb = (-parsb.CTE*4*x*z*(x2*parsb.V201 
                + y2*parsb.V111 + z2*parsb.V102 - parsb.V101))
        pyzb = (-parsb.CTE*4*y*z*(x2*parsb.V111 
                + y2*parsb.V021 + z2*parsb.V012 - parsb.V011))
    else:
        xl = xlmbd(x2,y2,z2,A2,B2,C2)
        ua3 = 1/(A2+xl)
        b3 = B2+xl
        ub3 = 1/b3
        uc3 = 1/(C2+xl)
        phi= np.arcsin(np.sqrt(ua3/parsb.UA2C2))

        [f,e] = elint(phi,parsb.XK)
        d2 = 2*np.sqrt(ua3*ub3*uc3)

        w100 = 2*(f-e)*parsb.UA2B2*parsb.SUA2C2
        w001 = (d2*b3 - 2*e*parsb.SUA2C2)*parsb.UB2C2
        w010 = d2 - w100 - w001

        w110 = (w010 - w100)*parsb.UA2B2
        w101 = (w001 - w100)*parsb.UA2C2
        w011 = (w001 - w010)*parsb.UB2C2

        w200 = (d2*ua3 - w110 - w101)*US3
        w020 = (d2*ub3 - w110 - w011)*US3
        w002 = (d2*uc3 - w011 - w101)*US3

        w111 = (w011 - w110)*parsb.UA2C2
        w210 = (w110 - w200)*parsb.UA2B2
        w201 = (w101 - w200)*parsb.UA2C2
        w120 = (w020 - w110)*parsb.UA2B2
        w021 = (w011 - w020)*parsb.UB2C2
        w102 = (w002 - w101)*parsb.UA2C2
        w012 = (w002 - w011)*parsb.UB2C2

        w300 = (d2*ua3*ua3 - w210 - w201)*0.2
        w030 = (d2*ub3*ub3 - w120 - w021)*0.2
        w003 = (d2*uc3*uc3 - w102 - w012)*0.2

        [dx,dy,dz]=derl(barra,x2,y2,z2)
        dw100=derw(barra,x2,y2,z2,1,0,0)
        dw010=derw(barra,x2,y2,z2,0,1,0)
        dw001=derw(barra,x2,y2,z2,0,0,1)

        dw110=derw(barra,x2,y2,z2,1,1,0)
        dw101=derw(barra,x2,y2,z2,1,0,1)
        dw011=derw(barra,x2,y2,z2,0,1,1)

        dw200=derw(barra,x2,y2,z2,2,0,0)
        dw020=derw(barra,x2,y2,z2,0,2,0)
        dw002=derw(barra,x2,y2,z2,0,0,2)

        dw111=derw(barra,x2,y2,z2,1,1,1)
        dw210=derw(barra,x2,y2,z2,2,1,0)
        dw201=derw(barra,x2,y2,z2,2,0,1)
        dw120=derw(barra,x2,y2,z2,1,2,0)
        dw021=derw(barra,x2,y2,z2,0,2,1)
        dw102=derw(barra,x2,y2,z2,1,0,2)
        dw012=derw(barra,x2,y2,z2,0,1,2)

        dw300=derw(barra,x2,y2,z2,3,0,0)
        dw030=derw(barra,x2,y2,z2,0,3,0)
        dw003=derw(barra,x2,y2,z2,0,0,3)

        AUX11=-(parsb.CTE * (w100 
                + x2*(5*x2*w300 + 6 * (y2*w210 - w200))
                + y2*(y2*w120 + 2*(z2*w111 - w110)) 
                + z2*(z2*w102 + 6*x2*w201 - 2*w101)))
        AUX21=-(dx*parsb.CTE*x*(dw100
                + x2*(x2*dw300+2*(y2*dw210-dw200))
                + y2*(y2*dw120+2*(z2*dw111-dw110))
                + z2*(z2*dw102+2*(x2*dw201-dw101))))
        pxxb = AUX11+AUX21

        AUX12=-(parsb.CTE * (w010
                + x2 * (x2*w210 + 6*y2*w120 - 2*w110)
                + y2 * (5*y2*w030 + 6*z2*w021 - 6*w020)
                + z2 * (2*x2*w111 + z2*w012 - 2*w011)))
        AUX22=-(dy*parsb.CTE*y*(dw010
                +x2*(x2*dw210+2*(y2*dw120-dw110))
                +y2*(y2*dw030+2*(z2*dw021-dw020))
                +z2*(z2*dw012+2*(x2*dw111-dw011))))
        pyyb = AUX12+AUX22

        AUX13=-(parsb.CTE * (w001
                + x2 * (x2*w201 + 2*y2*w111 - 2*w101)
                + y2 * (y2*w021 + 6*z2*w012 - 2*w011)
                + z2 * (5*z2*w003 + 6*x2*w102 - 6*w002)))
        AUX23=-(dz*parsb.CTE*z*(dw001
                + x2*(x2*dw201+2*(y2*dw111-dw101))
                + y2*(y2*dw021+2*(z2*dw012-dw011))
                + z2*(z2*dw003+2*(x2*dw102-dw002))))
        pzzb = AUX13+AUX23

        AUX14=-parsb.CTE*4*x*y*(x2*w210 + y2*w120 + z2*w111-w110)
        AUX24=-(dy*parsb.CTE*x*(dw100
                + x2*(x2*dw300+2*(y2*dw210-dw200))
                + y2*(y2*dw120+2*(z2*dw111-dw110))
                + z2*(z2*dw102+2*(x2*dw201-dw101))))
        pxyb = AUX14+AUX24

        AUX15=-parsb.CTE*4*x*z*(x2*w201 + y2*w111 + z2*w102-w101)
        AUX25=-(dz*parsb.CTE*x*(dw100
                + x2*(x2*dw300+2*(y2*dw210-dw200))
                + y2*(y2*dw120+2*(z2*dw111-dw110))
                + z2*(z2*dw102+2*(x2*dw201-dw101))))
        pxzb = AUX15+AUX25

        AUX16=-parsb.CTE*4*y*z*(x2*w111 + y2*w021 + z2*w012-w011)
        AUX26=-(dz*parsb.CTE*z*(dw001
                + x2*(x2*dw201+2*(y2*dw111-dw101))
                + y2*(y2*dw021+2*(z2*dw012-dw011))
                + z2*(z2*dw003+2*(x2*dw102-dw002))))
        pyzb = AUX16+AUX26
        
    return [pxxb,pyyb,pzzb,pxyb,pxzb,pyzb]

def disk(disco,x,y,z,omega):
    #Miyamoto
    x = x - disco.xd
    y = y - disco.yd
    AM = disco.a
    BM = disco.b
    GM = disco.GM
    EPS2=0
    QQ1=np.sin(EPS2)
    QQ2=np.cos(EPS2)
    XA=[QQ2*x-QQ1*z,y,QQ1*x+QQ2*z]
    R2 = XA[0]*XA[0]+XA[1]*XA[1]
    Z2 = XA[2]*XA[2]
    Z1 = np.sqrt(Z2+BM)
    AUX1=R2+(AM+Z1)*(AM+Z1)
    AUX2=AUX1*AUX1*AUX1
    AUX=np.sqrt(AUX2)
    AUX3=AUX1*AUX1*AUX1*AUX1*AUX1
    AUXX=np.sqrt(AUX3)
    der1 = 1/AUX1
    der = GM/(AUX)
    der2 = GM/(AUXX)
    AA=np.sqrt((Z2+BM)*(Z2+BM)*(Z2+BM))
    BB=(1+AM/Z1)
    pxxd = der * (1 - 3 * XA[0] * XA[0] * der1)
    pyyd = der * (1 - 3 * XA[1] * XA[1] * der1)
    pzzd =der*(BB*(1-3*Z2*der1*BB)-Z2*AM/AA)
    pxyd = - 3 * XA[0] * XA[1] * der2
    pxzd = - 3 * XA[0] * XA[2] * der2 * (1 + AM / Z1)
    pyzd = - 3 * XA[1] * XA[2] * der2 * (1 + AM / Z1)
    return [pxxd,pyyd,pzzd,pxyd,pxzd,pyzd]

def bulge(bulge,x,y,z,omega):
    x = x - bulge.xd
    y = y - bulge.yd
    AM = 0
    BM = bulge.b
    GM = bulge.GM
    EPS2=0
    QQ1=np.sin(EPS2)
    QQ2=np.cos(EPS2)
    XA=[QQ2*x-QQ1*z,y,QQ1*x+QQ2*z]
    R2 = XA[0]*XA[0]+XA[1]*XA[1]
    Z2 = XA[2]*XA[2]
    Z1 = np.sqrt(Z2+BM)
    AUX1=R2+(AM+Z1)*(AM+Z1)
    AUX2=AUX1*AUX1*AUX1
    AUX=np.sqrt(AUX2)
    AUX3=AUX1*AUX1*AUX1*AUX1*AUX1
    AUXX=np.sqrt(AUX3)
    der1 = 1/AUX1
    der = GM/(AUX)
    der2 = GM/(AUXX)
    AA=np.sqrt((Z2+BM)*(Z2+BM)*(Z2+BM))
    BB=(1+AM/Z1)
    pxxbl = der * (1 - 3 * XA[0] * XA[0] * der1)
    pyybl = der * (1 - 3 * XA[1] * XA[1] * der1)
    pzzbl = der*(BB*(1-3*Z2*der1*BB)-Z2*AM/AA)
    pxybl = - 3 * XA[0] * XA[1] * der2
    pxzbl = - 3 * XA[0] * XA[2] * der2 * (1 + AM / Z1)
    pyzbl = - 3 * XA[1] * XA[2] * der2 * (1 + AM / Z1)
    return [pxxbl,pyybl,pzzbl,pxybl,pxzbl,pyzbl]

def halo(halo,x,y,z,omega):
    #Plummer
    x = x - halo.xd
    y = y - halo.yd
    AM = 0
    BM = halo.b
    GM = halo.GM
    BM2 = BM*BM
    EPS2=0
    QQ1=np.sin(EPS2)
    QQ2=np.cos(EPS2)
    XA=[QQ2*x-QQ1*z,y,QQ1*x+QQ2*z]
    R2 = XA[0]*XA[0]+XA[1]*XA[1]
    Z2 = XA[2]*XA[2]
    Z1 = np.sqrt(Z2+BM2)
    AUX1=R2+(AM+Z1)*(AM+Z1)
    AUX2=AUX1*AUX1*AUX1
    AUX=np.sqrt(AUX2)
    AUX3=AUX1*AUX1*AUX1*AUX1*AUX1
    AUXX=np.sqrt(AUX3)
    der1 = 1/AUX1
    der = GM/(AUX)
    der2 = GM/(AUXX)
    AA=np.sqrt((Z2+BM2)*(Z2+BM2)*(Z2+BM2))
    BB=(1+AM/Z1)
    pxxh = der * (1 - 3 * XA[0] * XA[0] * der1)
    pyyh = der * (1 - 3 * XA[1] * XA[1] * der1)
    pzzh = der*(BB*(1-3*Z2*der1*BB)-Z2*AM/AA)
    pxyh = - 3 * XA[0] * XA[1] * der2
    pxzh = - 3 * XA[0] * XA[2] * der2 * (1 + AM / Z1)
    pyzh = - 3 * XA[1] * XA[2] * der2 * (1 + AM / Z1)
    return [pxxh,pyyh,pzzh,pxyh,pxzh,pyzh]

