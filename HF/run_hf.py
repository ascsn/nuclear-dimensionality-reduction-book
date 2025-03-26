import numpy as np
import matplotlib.pyplot as plt
import math
import os
import subprocess

def NM_to_t(nm_params):
    c13=1.0/3.0
    c23=2.0/3.0
    #VMASS =   1.249838574232270 UNEDF1
    VMASS = 1.25000101440918243
    kfconst =(1.50*math.pi**2)**(1.0/3.0)    # (3Pi^2/2)^(1/3)
    CK      = 3.0/5.0*kfconst**2
    hbzero   =   20.735530000000000
    mpi=   138.03/197.3
    #E,K,SMASS,RHO,ASS,LASS, Crdr, CrdJ
    E = nm_params["E"]
    K = nm_params["K"]
    SMASS = nm_params["SMASS"]
    RHO = nm_params["RHO"]
    ASS = nm_params["ASS"]
    LASS = nm_params["LASS"]
    Crdr = nm_params["Crdr"]
    CrdJ = nm_params["CrdJ"]
    tauc=CK*RHO**c23
    u=(kfconst/mpi)*RHO**c13
    rho2=RHO**2
    Crho = [0.0,0.0]
    Cdrho = [0.0,0.0]
    Ctau = [0.0,0.0]
    # set/calculate sigma
    sigma=((1.0/3.0)*(-K+tauc*hbzero*(-3.0+4.0*SMASS)-9.0*E)) \
             /(tauc*hbzero*(-3.0+2.0*SMASS)+3.0*E)
    
    Crho[0]=(c13*(tauc*hbzero*(-3.0+(2.0-3.0*sigma)*SMASS) \
        +3.0*(1.0+sigma)*E))/(sigma*RHO)
    Cdrho[0]=(c13*RHO**(-1.0-sigma)*(tauc*hbzero*(3.0-2.0*SMASS)\
        -3.0*E))/sigma
    Ctau[0]=(hbzero*(SMASS-1.0))/RHO
    
    Crho[1]=(27.0*ASS*(1.0+sigma)-9.0*LASS \
        +5.0*tauc*hbzero*(5.0-6.0*VMASS+3.0*sigma*(-4.0+3.0*VMASS)) \
        +RHO*(40.0*tauc*Ctau[0]-60.0*tauc*sigma*Ctau[0]))/(27.0*sigma*RHO)
    Cdrho[1]=-(RHO**(-1.0-sigma)*(27.0*ASS-9.0*LASS \
        +5.0*tauc*hbzero*(5.0-6.0*VMASS) \
        +40.0*tauc*RHO*Ctau[0]))/(27.0*sigma)
    Ctau[1]=(hbzero-hbzero*VMASS+RHO*Ctau[0])/RHO
    
    t0  =  (8.0/3)*Crho[0]
    t1  =  4.0/3.0*(Ctau[0]-4.0*Crdr[0])
    t2  =  4.0/3.0*(3.0*Ctau[0]-6.0*Ctau[1]+4.0*Crdr[0]-8.0*Crdr[1])
    t3  =  16.0*Cdrho[0]
    x0  = -0.50*(3.0*Crho[1]/Crho[0]+1.0)
    x1  =  2.0*(-Ctau[0]-3.0*Ctau[1]+4.0*Crdr[0]+12.0*Crdr[1])/t1/3.0
    x2  = -2.0*(3.0*Ctau[0]-15.0*Ctau[1]+4.0*Crdr[0]-20.0*Crdr[1])/t2/3.0
    x3  = -0.50*(3.0*Cdrho[1]/Cdrho[0]+1.0)
    # t4  =  (CrdJ[1]-CrdJ[0]) / 2.0
    # t4p = -CrdJ[1]*4.0
    t_params = {"t0":t0,"x0":x0,"t1":t1,"x1":x1,"t2":t2,"x2":x2,"t3":t3,"x3":x3,"cso0":CrdJ[0],"cso1":CrdJ[1],"sig":1/sigma}
    #print(t_params)
    return t_params


def make_input(params,restart,directory):
    in_template = "in.template"
    f1 = open(in_template, 'r')
    
    int_string = "&interaction t0={t0},x0={x0},t1={t1},x1={x1},t2={t2},x2={x2},t3={t3},x3={x3},cso0={cso0},cso1={cso1},sig={sig},icoul=1,icm=0,j2terms=.false. / \n".format(\
    t0=params["t0"],x0=params["x0"],t1=params["t1"],x1=params["x1"],t2=params["t2"],x2=params["x2"],t3=params["t3"],x3=params["x3"],cso0=params["cso0"],cso1=params["cso1"],sig=params["sig"])
    if(restart):
        io_string = "&io printwf=F,printdens=F,restartwf=T,printhpsi=T /"
    else:
        io_string = "&io printwf=F,printdens=T,restartwf=F,printhpsi=F /"
    in_out = directory+"/in"
    f2 = open(in_out, 'w')
    f1.seek(0)
    f2.write(f1.read())
    f2.write(int_string)
    f2.write(io_string)
    f2.write("\n")

    f1.close()
    f2.close()

    return

def run_HF(tag,restart):
    noNaN = 1
    neut = 0
    if(restart == 0):
        os.system("cp in results/in."+tag)
        #print(os.system("../HF/./main"))
        output = subprocess.run(["./main"], capture_output=True)
        #print(output)
        #print(subprocess.run("../HF/./main",shell=True))
        os.system("mv out results/out."+tag)
        # os.system("mv wf_numpy.bin results/wf_numpy.bin."+tag)
        # noNaN = os.system("grep -q NaN results/out."+tag)
        # neut = os.system('grep -q "Total Neutrons =   28.000" results/out.'+tag)
        # neut = neut+os.system('grep -q "Total Neutrons =   27.999" results/out.'+tag)
    else:
        output = subprocess.run(["./main"], capture_output=True)
        #os.system("mv out results/out_rst."+tag)
        #os.system("mv hpsi.bin results/hpsi.bin."+tag)
        #print(output)
        #noNaN = os.system("grep -q NaN results/out_rst."+tag)
    # if (noNaN == 0):
    #     print("NaNs in run {tag}!".format(tag=tag))
    #     return 1
    # if(neut > 256):
    #     print("Wrong number of neutrons in run {tag}!".format(tag=tag))
    #     return 1

    return 0


nmdefault_params = {"E":-15.8,"K":220.0,"SMASS":0.992423332283364,"RHO":0.158706769332587,"ASS":28.986789057772100,"LASS":40.004790480413600,"Crdr":[-45.1351310222373,-145.382167908057],"CrdJ":[-74.0263331764599,-35.6582611147917]}

restart = False
directory = '.'
tag = 'default'

make_input(NM_to_t(nmdefault_params),restart,directory)
run_HF(tag,restart)