def read_ini_peqs_file(myfile):
    ini_peqs_data = ""
    with open(file=myfile, mode="r") as f:
        ini_peqs_data = f.readlines() #5 L_i punts de 6D cadascun (p,q)\in R^6
    ini_peqs_data = [x.strip().split(" ") for x in ini_peqs_data]
    ini_peqs = []
    for i in ini_peqs_data:
        peq_i = []
        for j in i:
            if j=="":
                continue
            peq_i.append(float(j))
        ini_peqs.append(peq_i)
    return ini_peqs
