from numpy import pi, sin, cos, sqrt

def bar(x,y,z,barra):
    '''
    Ferrers
    '''
    x = x-barra.xd
    y = y-barra.yd

    a = barra.a
    b = barra.b
    c = barra.c
    US3 = 1/3
    US6 = 1/6     
    GM = barra.GM
    X2 = x*x
    Y2 = y*y
    Z2 = z*z
    M2 = x*x/(a*a) + y*y/(b*b) + z*z/(c*c)
    rho0 = (105/(32*pi))*(GM/(a*b*c))

    if M2 <= 1:
        dens = rho0*(1-M2)**2
    else:
        dens = 0

    return dens

def disk(x,y,z,disco):
    '''
    Miyamoto
    '''
    x = x - disco.xd
    y = y - disco.yd

    GM=disco.GM
    AM=disco.a
    BM=disco.b #%Estamos en 2dimensiones. Tomamos BM = B^2, como B=1 aqui no hay problema
    #pero cuidado si B vale distinto de 1!!!!!
    # phid = -GMd./sqrt(x.*x+y.*y+(A+sqrt(B*B+z.*z)).*(A+sqrt(B*B+z.*z)))

    #TODO
    '''
    EPS2=0 #? 

    QQ1=sin(EPS2)
    QQ2=cos(EPS2)

    XA(1)=QQ2*x-QQ1*z
    XA(2)=y
    XA(3)=QQ1*x+QQ2*z
    '''

    R2 = x*x+y*y
    Z2 = z*z
    T1 = BM*GM/(4*pi)
    AUX1 = (AM+3*sqrt(Z2+BM))*(AM+sqrt(Z2+BM))**2
    Tsup = AM*R2+AUX1
    AUX2 = R2+(AM+sqrt(Z2+BM))*(AM+sqrt(Z2+BM))
    Tinf = (AUX2**2.5)*((Z2+BM)**1.5)
    dens = T1*Tsup/Tinf
    return dens

def bulge(x,y,z,bulge):
    '''
    Plummer
    '''
    x = x - bulge.xd
    y = y - bulge.yd

    b1 = bulge.b
    GM = bulge.GM
    R2 = x*x+y*y+z*z
    T1 = 3*GM/(4*pi*b1**3)
    T2 = 1 + R2/b1**2
    dens = T1*(T2**(-2.5))
    return dens

def halo(x,y,z,halo):
    '''
    Plummer
    '''
    x = x - halo.xd
    y = y - halo.yd

    b1 = halo.b
    GM = halo.GM
    R2 = x*x+y*y+z*z
    T1 = 3*GM/(4*pi*b1**3)
    T2 = 1 + R2/b1**2
    dens = T1*(T2**(-2.5))
    return dens

def efectiva(x,y,z,mbarra,mdisk,mbulge,mhalo):
    '''
    Returns the "effective" density at a given x,y,z point.
    Originally called "densidad.m" in the legacy code,
    renamed for consistency with pot.efectivo
    '''
    DENSB = bar(x,y,z,mbarra)      #Densidad de barra Ferrers
    DENSD = disk(x,y,z,mdisk)     #Densidad de Miyamoto
    DENSBL = bulge(x,y,z,mbulge)   #Densidad del bulge esferico de Plummer 
    DENSH = halo(x,y,z,mhalo)     #Densidad del halo esferico de Plummer

    dens = DENSH + DENSB + DENSBL + DENSD
    return dens


