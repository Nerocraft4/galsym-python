from maths.helpers import elint, xlmbd
import numpy as np

def bar(barra,parsb,x,y,z,omega):
    #Ferrer's potential
    a,b,c = barra.a,barra.b,barra.c
    A2=a*a
    B2=b*b
    C2=c*c
    x = x-barra.xd
    y = y-barra.yd
    US3 = 1/3
    X2 = x*x
    Y2 = y*y
    Z2 = z*z
    pxb,pyb,pzb = [0,0,0]

    if X2*parsb.UA2+Y2*parsb.UB2+Z2*parsb.UC2 <= 1:
        px = -parsb.CTE*x*(
            (parsb.V100+X2*(X2*parsb.V300+2*(Y2*parsb.V210-parsb.V200))  
            +Y2*(Y2*parsb.V120+2*(Z2*parsb.V111-parsb.V110))             
            +Z2*(Z2*parsb.V102+2*(X2*parsb.V201-parsb.V101))))             
        py = -parsb.CTE*y*(
            (parsb.V010+X2*(X2*parsb.V210+2*(Y2*parsb.V120-parsb.V110))  
            +Y2*(Y2*parsb.V030+2*(Z2*parsb.V021-parsb.V020))             
            +Z2*(Z2*parsb.V012+2*(X2*parsb.V111-parsb.V011))))
        pz = -parsb.CTE*z*(
            (parsb.V001+X2*(X2*parsb.V201+2*(Y2*parsb.V111-parsb.V101))  
            +Y2*(Y2*parsb.V021+2*(Z2*parsb.V012-parsb.V011))          
            +Z2*(Z2*parsb.V003+2*(X2*parsb.V102-parsb.V002))))
    else:
        XL = xlmbd(X2,Y2,Z2,A2,B2,C2)
        UA3 = 1.0/(A2+XL)
        B3 = B2+XL
        UB3 = 1.0/B3
        UC3 = 1.0/(C2+XL)
        PHI= np.arcsin(np.sqrt(UA3/parsb.UA2C2))
        [F,E] = elint(PHI,parsb.XK)
        D2 = 2*np.sqrt(UA3*UB3*UC3)

        W100 = 2*(F-E)*parsb.UA2B2*parsb.SUA2C2
        W001 = (D2*B3 - 2*E*parsb.SUA2C2)*parsb.UB2C2
        W010 = D2 - W100 - W001


        W110 = (W010 - W100)*parsb.UA2B2
        W101 = (W001 - W100)*parsb.UA2C2
        W011 = (W001 - W010)*parsb.UB2C2

        W200 = (D2*UA3 - W110 - W101)*US3
        W020 = (D2*UB3 - W110 - W011)*US3
        W002 = (D2*UC3 - W011 - W101)*US3

        W111 = (W011 - W110)*parsb.UA2C2
        W210 = (W110 - W200)*parsb.UA2B2
        W201 = (W101 - W200)*parsb.UA2C2
        W120 = (W020 - W110)*parsb.UA2B2
        W021 = (W011 - W020)*parsb.UB2C2
        W102 = (W002 - W101)*parsb.UA2C2
        W012 = (W002 - W011)*parsb.UB2C2

        W300 = (D2*UA3*UA3 - W210 - W201)*0.2
        W030 = (D2*UB3*UB3 - W120 - W021)*0.2
        W003 = (D2*UC3*UC3 - W102 - W012)*0.2

        px = -parsb.CTE * x*( 
            (W100+X2*(X2*W300+2*(Y2*W210-W200))
            +Y2*(Y2*W120+2*(Z2*W111-W110))
            +Z2*(Z2*W102+2*(X2*W201-W101))))
        py = -parsb.CTE * y*(
            (W010+X2*(X2*W210+2*(Y2*W120-W110))
            +Y2*(Y2*W030+2*(Z2*W021-W020))
            +Z2*(Z2*W012+2*(X2*W111-W011))))
        pz = -parsb.CTE * z*(
            (W001+X2*(X2*W201+2*(Y2*W111-W101))
            +Y2*(Y2*W021+2*(Z2*W012-W011))
            +Z2*(Z2*W003+2*(X2*W102-W002))))
    return [px,py,pz]

def disk(disco,x,y,z,omega):
    #Miyamoto's potential
    AM = disco.a
    BM = disco.b
    GM = disco.GM
    d2 = x*x+y*y+(AM+np.sqrt(BM*BM+z*z))*(AM+np.sqrt(BM*BM+z*z))
    d3 = d2*np.sqrt(d2)
    c = GM/d3
    pxd = x*c
    pyd = y*c
    pzd =(1+AM/np.sqrt(BM*BM+z*z))*(z*c)
    return [pxd,pyd,pzd]

def bulge(bulge,x,y,z,omega):
    #Plummer potential
    x = x - bulge.xd
    y = y - bulge.yd
    AM = 0
    BM = bulge.b
    GM = bulge.GM
    d2 = x*x+y*y+(AM+np.sqrt(BM*BM+z*z))*(AM+np.sqrt(BM*BM+z*z))
    d3 = d2*np.sqrt(d2)
    c=GM/d3
    pxbl = x*c
    pybl = y*c
    pzbl =(1+AM/np.sqrt(BM*BM+z*z))*(z*c)
    return [pxbl,pybl,pzbl]

def halo(halo,x,y,z,omega):
    #Plummer potential
    x = x - halo.xd
    y = y - halo.yd
    AM = 0
    BM = halo.b
    GM = halo.GM
    d2 = x*x+y*y+(AM+np.sqrt(BM*BM+z*z))*(AM+np.sqrt(BM*BM+z*z))
    d3 = d2*np.sqrt(d2)
    c=GM/d3
    pxh = x*c
    pyh = y*c
    pzh =(1+AM/np.sqrt(BM*BM+z*z))*(z*c)
    return [pxh,pyh,pzh]

