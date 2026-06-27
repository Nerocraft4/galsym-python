from numpy import round

'''
Set of functions to get names of files and folders
'''

def datadir(params: list) -> str:
    [barra, disco, bulge, halo, parsb] = params
    p = r"halomass_"
    p+= str(halo.GM)
    p+= "_bulgemass_"
    p+= str(bulge.GM)
    p+= "_halo_xd_"
    p+= str(halo.xd)
    p+= "_yd_"
    p+= str(halo.yd)
    p+= r"/"
    p+= "bulgexd_"
    p+= str(bulge.xd)
    p+= r"/"
    return p

def ini_guess_pequi_file(params: list) -> str:
    [barra, disco, bulge, halo, parsb] = params
    a = ['peqi_entrada3','_halo_xd_',str(halo.xd),
         '_yd_',str(halo.yd),
         '_bulgexd_',str(bulge.xd)]
    a = "".join(a).replace(".","_")+".dat"
    return a

def arxol(punt: int, xkk: float, fix: str) -> str:#maybe rename
    a = "orlinL"
    a+= str(punt+1)
    a+= "_xkk"
    a+= str(round(xkk,2)).replace(".","_")
    a+= "_"
    a+= fix
    return a
