from classes import Barra, Halo, Bulge, Disk

def initializer(arxi):
    '''
    arxi: file where the data for the Bar, Halo, Bulge and Disk is stored
    '''
    data = ""
    with open(file=arxi,mode="r") as f:
        data = f.readlines()
    data = [d.strip() for d in data]
    xd,yd = [0,0] #default position

    #start by initializing the Bar
    a,b,c,GM = [float(x) for x in data[0].split(" ")]
    omega = float(data[1].split(" ")[0])
    eps = float(data[3].split(" ")[1])
    barra = Barra(xd=0,yd=0,a=a,b=b,c=c,GM=GM,omega=omega,eps=eps)

    #then the Disk
    a,b,GM = [float(x) for x in data[2].split(" ")]
    disco = Disk(xd=0,yd=0,a=a,b=b,GM=GM)

    #the Bulge
    b, GM = [float(x) for x in data[4].split(" ")]
    bulge = Bulge(xd=0,yd=0,b=b,GM=GM)

    #and the Halo
    b, GM = [float(x) for x in data[5].split(" ")]
    halo = Halo(xd=0,yd=0,b=b,GM=GM)

    #TODO
    parsb = 0 #is this a matrix? or a collection of params?

    return barra, disco, bulge, halo, parsb
