#!/hetghome/hetgsoft/anaconda3/bin/python3.9
#sed -i 's/\r//g' cal_insta_map.py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math as m
import math
from scipy.optimize import fsolve
from scipy import linalg
import csv
import pandas as pd
import sys
import numpy as np
import os
from scipy.integrate import RK45
import time

#######################################
#Lmu = 200
LLmu = [300]
AR = [0.7]
theta = m.pi/2
alpha=5/3
ar = 0.7
omg = 1
omga = -omg
Omg = 0
dz = 0.01
hz =5
v = 1
#########################################################
deg = theta/2/m.pi*360
vxp=m.sin(theta)*v
vxn=m.sin(theta)*(-v)
vz=m.cos(theta)*v
zwpi = False
######################################################
def hvv(vxp,vxn,mu):
    if zwpi == True:
        return mu*(1-(vz*vz+vxp*vxn))/math.pi/2
    else:
        return mu*(1-(vz*vz+vxp*vxn))#2*mu*(E**3*g*(1-(vz*vz+vxp*vxn)))#/(2*m.pi)**3
def hvva(vxp,vxn,mu):
    if zwpi == True:
        return mu*(1-(vz*vz+vxp*vxn))*(-1)*ar/math.pi/2
    else:
        return mu*(1-(vz*vz+vxp*vxn))*(-1)*ar
def fmu(z):
    mu= 10**4*math.exp(-0.3*z)
    return mu
def Lda(vxp,vxn,vz,mu):
    #2*mu*(E**3*g*(1-(vz*vz+vxp*vxn))*(1-ar))#/(2*m.pi)**3       #2*mu*(1-(vz*vz+vxp*vxn))*(E**2)*(1-ar)/((2*m.pi)**1)
    return mu*(1-ar)*(1-(vz*vz+vxp*vxn))
def lda0(mu):
    return Lmu*mu
#########################################################
def fim_eig(z,modk):
    LHSM = []
    #sys.exit()
    print(omg,'omg')
    print(omga,'omga')
    for i in range(4):
        A = []
        for j in range(4):
            A.append(0+0j)
        LHSM.append(A)
    for i in range(len(LHSM)):
        for j in range(len(LHSM[i])):
            if i == j:
                if i%4 ==0:
                    LHSM[i][j] = LHSM[i][j]+(-Omg+omg+vxp*modk)+Lda(vxp,vxn,vz,fmu(z))+lda0(fmu(z))
                if i%4 == 1:
                    LHSM[i][j] = LHSM[i][j]+(-Omg+omg+vxn*modk)+Lda(vxp,vxn,vz,fmu(z))+lda0(fmu(z))
                if i%4 == 2:
                    LHSM[i][j] = LHSM[i][j]+(-Omg+omga+vxp*modk)+Lda(vxp,vxn,vz,fmu(z))+lda0(fmu(z))
                if i%4 == 3:
                    LHSM[i][j] = LHSM[i][j]+(-Omg+omga+vxn*modk)+Lda(vxp,vxn,vz,fmu(z))+lda0(fmu(z))
        j =  i   #see test part
        if i%4 ==0:
            if j+1 < 4:
                LHSM[i][j+1] = LHSM[i][j+1]-hvv(vxp,vxn,fmu(z))
            if j+3 < 4:
                LHSM[i][j+3] = LHSM[i][j+3]-hvva(vxp,vxn,fmu(z))
        if i%4 == 2:
            if j-1 > -1:
                LHSM[i][j-1] = LHSM[i][j-1]-hvv(vxp,vxn,fmu(z))
            if j+1 < 4:
                LHSM[i][j+1] = LHSM[i][j+1]-hvva(vxp,vxn,fmu(z))
        if i%4 == 1:
            if j-1 > -1:
                LHSM[i][j-1] = LHSM[i][j-1]-hvv(vxp,vxn,fmu(z))
            if j+1 < 4:
                LHSM[i][j+1] = LHSM[i][j+1]-hvva(vxp,vxn,fmu(z))
        if i%4 == 3:
            if j-1 > -1:
                LHSM[i][j-1] = LHSM[i][j-1]-hvva(vxp,vxn,fmu(z))
            if j-3 > -1:    
                LHSM[i][j-3] = LHSM[i][j-3]-hvv(vxp,vxn,fmu(z))
    NCFM = LHSM
    NCM = np.array(NCFM)# / vz
    #print("zwpi",zwpi)
    l, vec= linalg.eig(NCM)
    Al = []
    Cii = []#"check if in"
    #iniQ = []
    for j in range(len(l)):
        Al.append(abs(((l[j]).imag)))
    Ml = np.max(Al)
    print(Al)
    """for j in range(len(Al)):
        if Al[j] == np.max(Al):
            if len(Cii) == 0:
                iniQ = vec[:,j]"""
    #sys.exit()
    #print("max eig val",Ml)
    #ML.append(Ml)
    del LHSM
    return Ml
Lmu=300
omg=1
omga=-1
print(fim_eig((np.log(100)-np.log(10000))/(-0.3), 0))
omg=-1
omga=1
print(fim_eig((np.log(100)-np.log(10000))/(-0.3),0))
sys.exit()
for g in AR:
    ar = g
    for h in LLmu:
        Lmu = h
        Lmodk = -10000

        Hmodk = 10000
        hmu=1000
        lmu=10 #low mu
        Lz = (np.log(hmu)-np.log(10000))/(-0.3)
        Hz = (np.log(lmu)-np.log(10000))/(-0.3)
        print(Lz,Hz)
        #sys.exit()
        DATA = []
        modk = Lmodk
        #fn = "/hetghome/jordan/Cose_nu/research_experience_program/COSEnu/f_insta_map/lmu_"+str(int(Lmodk/(10**(5))))+"_hmu_"+str(int(Hmodk/(10**(5))))+"_Lz_"+str(Lz)+"_Hz_"+str(Hz)+"_Lmu_"+str(Lmu)+"_ar_"+str(ar)+"_deg_"+str(deg)
        fn = "/hetghome/jordan/Cose_nu/research_experience_program/COSEnu/f_insta_map/Lk_"+str(int(Lmodk/(10**(4))))+"_Hk_"+str(int(Hmodk/(10**(4))))+"_lmu_"+str(lmu)+"_hmu_"+str(hmu)+"_Lmu_"+str(Lmu)+"_ar_"+str(ar)+"_deg_"+str(deg)        
        fnc = fn+'_f'
        if os.path.isdir(fnc):
            print(fn)
            print('exist')
            #continue
        else:
            os.mkdir(fnc)
            TRASH = 1
        ni=10000
        nj=10000
        for i in range(ni):
            z = Lz
            st = time.time()
            for j in range(nj):
                ival = fim_eig(z,modk)
                mu=lmu + (hmu-lmu)* j /nj
                z= (np.log(mu)-np.log(10000))/(-0.3)
                DATA.append([z,modk,ival])
                #z = Lz+(Hz-Lz) * j /nj
                #print(z,modk,ival)
                if z > Hz+10**-3:
                    break
            modk = Lmodk + (Hmodk-Lmodk) * i/ni
            print("time --------------------------------------------------",time.time()-st)
            if modk > Hmodk+1:
                break
        np.save(fn,DATA)
        #np.save("/hetghome/jordan/NeuOsc/PAP_kmin/f_insta_map/Lk_"+str(int(Lmodk/(10**(-5))))+"_Hk_"+str(int(Hmodk/(10**(-5))))+"_Lz_"+str(Lz)+"_Hz_"+str(Hz),DATA)
    