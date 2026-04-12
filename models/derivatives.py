#TODO fix notation (miya->disk)


def der1bar(x,y,z,omega):
    pxb,pyb,pzb = [0,0,0]
    return [pxb,pyb,pzb]

def der1miya(x,y,z,omega):
    pxm,pym,pzm = [0,0,0]
    return [pxm,pym,pzm]

def der1bulge(x,y,z,omega):
    pxbl,pybl,pzbl = [0,0,0]
    return [pxbl,pybl,pzbl]

def der1halo(x,y,z,omega):
    pxh,pyh,pzh = [0,0,0]
    return [pxh,pyh,pzh]

def der2bar(x,y,z,omega):
    pxxb,pyyb,pzzb,pxyb,pxzb,pyzb = [0,0,0,0,0,0]
    return [pxxb,pyyb,pzzb,pxyb,pxzb,pyzb]

def der2miya(x,y,z,omega):
    pxxm,pyym,pzzm,pxym,pxzm,pyzm = [0,0,0,0,0,0]
    return [pxxm,pyym,pzzm,pxym,pxzm,pyzm]

def der2bulge(x,y,z,omega):
    pxxbl,pyybl,pzzbl,pxybl,pxzbl,pyzbl = [0,0,0,0,0,0]
    return [pxxbl,pyybl,pzzbl,pxybl,pxzbl,pyzbl]

def der2halo(x,y,z,omega):
    pxxh,pyyh,pzzh,pxyh,pxzh,pyzh = [0,0,0,0,0,0]
    return [pxxh,pyyh,pzzh,pxyh,pxzh,pyzh]

