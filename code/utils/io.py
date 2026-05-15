from .namer import ini_guess_pequi_file, datadir
import os
import numpy as np
from numpy.typing import NDArray

def read_ini_peqs(inputdir: str, params: dict) -> NDArray[np.float_]:
    my_ini_guess_pequi_file = inputdir+"/"+ini_guess_pequi_file(extract_galparams(params))
    print("Attempting to read initial guess for eq points from",my_ini_guess_pequi_file)
    ini_peqs_data = ""
    try:
        print("Found initial seed for eq points")
        print("Reading given points")
        with open(file=my_ini_guess_pequi_file, mode="r") as f:
            ini_peqs_data = f.readlines() #5 L_i punts de 6D cadascun (p,q)\in R^6     
    except:
        print("Initial guess for eq points file not found")
        print("Reading default points. Might not converge")
        with open(file="DEFAULT", mode="r") as f:
            ini_peqs_data = f.readlines()
    ini_peqs_data = [x.strip().split(" ") for x in ini_peqs_data]
    ini_peqs = []
    for i in ini_peqs_data:
        peq_i = []
        for j in i:
            if j=="":
                continue
            peq_i.append(float(j))
        ini_peqs.append(peq_i[:3])
    return np.array(ini_peqs)

def make_data_dir(data: str, params: dict) -> None:
    mydatadir = data + "/" + datadir(extract_galparams(params))    
    print("Creating new directory in",mydatadir)
    try:
        os.makedirs(mydatadir)
        print("New directory successfully created")
    except:
        print("The directory already existed")

def extract_galparams(galparams: dict) -> list:
    disco = galparams["disco"]
    barra = galparams["barra"]
    bulge = galparams["bulge"]
    halo  = galparams["halo"]
    parsb = galparams["parsb"]
    return [barra,disco,bulge,halo,parsb]

def pack_galparams(galparams: list) -> dict:
    [barra,disco,bulge,halo,parsb] = galparams
    galparams_dict = {"barra":barra,"disco":disco,"bulge":bulge,"halo":halo,"parsb":parsb}
    return galparams_dict
