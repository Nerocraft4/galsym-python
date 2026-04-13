'''
Set of functions to name files and folders
'''

def name_datadir(xydhalo,xdbulge,halo,bulge):
    p = r"halomass_"
    p+= str(halo.GM)
    p+= "_bulgemass_"
    p+= str(bulge.GM)
    p+= "_halo_xd_"
    p+= str(xydhalo[0])
    p+= "_yd_"
    p+= str(xydhalo[1])
    p+= r"/"
    p+= "bulgexd_"
    p+= str(xdbulge)
    p+= r"/"
    return p

def name_ini_guess_pequi_file(halo,xdbulge):
    a = ['peqi_entrada3','_halo_xd_',str(halo.xd),
         '_yd_',str(halo.yd),
         '_bulgexd_',str(xdbulge)]
    a = "".join(a).replace(".","_")+".dat"
    return a
