from dataclasses import dataclass

@dataclass
class Generic:
    xd: float           #x-displacement from center of masses
    yd: float           #y-displacement ""
    GM: float           #Mass times Grav constant (?)

@dataclass
class Barra(Generic):
    a: float            #semimajor axis
    b: float            #semiminor axis1
    c: float            #semiminor axis2
    omega: float        #omega(?)
    eps: float          #epsilon

@dataclass
class Halo(Generic):
    b: float            #b value (radius?)

@dataclass
class Bulge(Generic):
    b: float            #b value (radius?)

@dataclass
class Disk(Generic):
    a: float            #a (radius?)
    b: float            #b value (radius?)

@dataclass
class ParsB():
    '''
    a: float   not needed i think?
    b: float
    c: float
    '''
    UA2: float
    UB2: float
    UC2: float
    CTE: float
    UA2B2: float
    UA2C2: float
    UB2C2: float
    SUA2C2: float
    XK: float
    V000: float
    
    #unique 1st-degree derivatives?
    V100: float 
    V001: float 
    V010: float 

    #unique 2nd-degree derivatives
    V110: float 
    V101: float 
    V011: float 

    V200: float 
    V020: float 
    V002: float 

    #unique 3rd-degree derivatives
    V111: float 
    V210: float 
    V201: float 
    V120: float 
    V021: float 
    V102: float 
    V012: float 
    V300: float 
    V030: float 
    V003: float


