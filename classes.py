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
