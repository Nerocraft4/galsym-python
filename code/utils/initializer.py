from .classes import Barra, Halo, Bulge, Disk, ParsB
from maths.helpers import elint
import numpy as np
from .io import pack_galparams

def initialize(arxi: str) -> [dict,list]:
    '''
    In charge of initializing galactic parameters.
    Input:
        arxi: file where the data for the Bar, Halo, Bulge and Disk is stored
    Output:
        galparams: dictionary conatining all information from galactic parameters
            {"barra":barra,"disco":disco,"bulge":bulge,"halo":halo,"parsb":parsb}
        displacements: list of 3D displacements of each part of the galaxy
    Potential upgrades: eventually upgrade to a JSON file please xd
    '''
    data = ""
    with open(file=arxi,mode="r") as f:
        data = f.readlines()
    data = [d.strip() for d in data]

    #get displacements for all values
    dbar = [float(x) for x  in data[-4].split(" ")[:3]]
    ddisk = [float(x) for x  in data[-3].split(" ")[:3]]
    dbulge = [float(x) for x  in data[-2].split(" ")[:3]]
    dhalo = [float(x) for x  in data[-1].split(" ")[:3]]

    #start by initializing the Bar
    a,b,c,GM = [float(x) for x in data[1].split(" ")[:4]]
    omega = float(data[2].split(" ")[0])
    eps = float(data[4].split(" ")[1])
    barra = Barra(xd=dbar[0],yd=dbar[1],zd=dbar[2],a=a,b=b,c=c,GM=GM,omega=omega,eps=eps)

    #then the Disk
    a,b,GM = [float(x) for x in data[3].split(" ")[:3]]
    disco = Disk(xd=ddisk[0],yd=ddisk[1],zd=ddisk[2],a=a,b=b,GM=GM)

    #the Bulge
    b, GM = [float(x) for x in data[5].split(" ")[:2]]
    bulge = Bulge(xd=dbulge[0],yd=dbulge[1],zd=dbulge[2],b=b,GM=GM)

    #and the Halo
    b, GM = [float(x) for x in data[6].split(" ")[:2]]
    halo = Halo(xd=dhalo[0],yd=dhalo[1],zd=dhalo[2],b=b,GM=GM)

    #ParsB section
    US3=1/3;
    US6=1/6;
    A2 = barra.a*barra.a
    B2 = barra.b*barra.b
    C2 = barra.c*barra.c
    
    UA2 = 1/A2
    UB2 = 1/B2
    UC2 = 1/C2
    CTE = -barra.GM*105/16
    UA2B2 = 1/(A2 - B2)
    UA2C2 = 1/(A2 - C2)
    UB2C2 = 1/(B2 - C2)
    SUA2C2 = np.sqrt(UA2C2)
    PHI = np.arcsin(np.sqrt(1 - C2*UA2))

    XK = np.sqrt(UA2C2*(A2-B2))

    [F,E] = elint(PHI,XK)

    D2 = 2*np.sqrt(UA2*UB2*UC2)

    #Wijk inside the ellipsoid are constant!

    V000 = 2*F*SUA2C2
    V100 = 2*(F-E)*UA2B2*SUA2C2
    V001 = (D2*B2 - 2*E*SUA2C2)*UB2C2
    V010 = D2 - V100 - V001

    V110 = (V010 - V100)*UA2B2;
    V101 = (V001 - V100)*UA2C2;
    V011 = (V001 - V010)*UB2C2;
    V200 = (D2*UA2 - V110 - V101)*US3;
    V020 = (D2*UB2 - V110 - V011)*US3;
    V002 = (D2*UC2 - V011 - V101)*US3;

    V111 = (V011 - V110)*UA2C2;
    V210 = (V110 - V200)*UA2B2;
    V201 = (V101 - V200)*UA2C2;
    V120 = (V020 - V110)*UA2B2;
    V021 = (V011 - V020)*UB2C2;
    V102 = (V002 - V101)*UA2C2;
    V012 = (V002 - V011)*UB2C2;
    V300 = (D2*UA2*UA2 - V210 - V201)*.2;
    V030 = (D2*UB2*UB2 - V120 - V021)*.2;
    V003 = (D2*UC2*UC2 - V102 - V012)*.2;

    parsb = ParsB(UA2,UB2,UC2,CTE,UA2B2,UA2C2,UB2C2,SUA2C2,XK,
                  V000,V100,V001,V010,V110,V101,V011,V200,V020,V002,
                  V111,V210,V201,V120,V021,V102,V012,V300,V030,V003)

    galparams = pack_galparams([barra,disco,bulge,halo,parsb])
    return galparams, [dbar,ddisk,dbulge,dhalo]
