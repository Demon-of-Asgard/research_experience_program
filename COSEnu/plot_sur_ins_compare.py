#!/hetghome/hetgsoft/anaconda3/bin/python3.9
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math as m
import math
from scipy.fft import fft, fftfreq, fftshift
"""from scipy.optimize import fsolve
from scipy import linalg"""
import csv
import pandas as pd
import sys
import numpy as np
import os
from scipy.integrate import RK45
import time
import seaborn as sns
from os.path import exists
from matplotlib import rc
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
rc('font',**{'family':'serif','serif':['Roman'],'size':20})
plt.rc('legend', fontsize=10) 
plt.rc('xtick', labelsize=15)    
plt.rc('ytick', labelsize=15)
rc('text', usetex=True)
#################################################################
from configparser import ConfigParser
########################################
#insta growth rate search
infn = '/hetghome/jordan/Cose_nu/research_experience_program/COSEnu/f_insta_map/Lk_-1_Hk_1_lmu_10_hmu_1000_Lmu_300_ar_0.7_deg_90.0.npy'
######################################
path = os.getcwd()
#print(path)
nax=0
pst=True #physical time
for h in os.listdir(path):
    if 'flat' not in h:
        print('run123 conti')
        continue
    """if 'hetg' not in h:
        print('hetg pui conti')
        continue"""
    if 'N1' not in h:
        print('no N1 conti')
        continue
    if 'pb' not in h:
        print('pb conti')
        continue
    cc=0
    if '_mu_' in h:
        #print(h)
        if 'pb' in h:
            pb=h[h.find('pb')+2:]
        if 'flat' in h:
            mu=h[h.find('mu_')+3:h.find('_Lmu')]
        else:
            mu=h[h.find('mu_')+3:h.find('_ar')]
        xli = 1200
        if os.path.isdir(path+'/'+h)==False:
            continue
        for i in os.listdir(path+'/'+h):
            #print(i)
            if os.path.isdir(path+'/'+h+'/'+i):
                ###########################################
                confn = path+'/'+h+'/'+i+'/job.config'
                print(confn)
                f=open(confn,'r')
                for line in f.readlines():
                    if "dz" in line:
                        dz=float(line[line.find(":")+1:])
                    if "nz" in line:
                        nz=float(line[line.find(":")+1:])
                        #print(line)
                        #print(dz)
                kmax = max(fftshift(fftfreq(int(nz),dz)))*float(mu)
                print('kmax',kmax)
                print(fftshift(fftfreq(int(nz),dz)))
                #sys.exit()
                ############################################
                print(path+'/'+h+'/'+i)
                #pign = path+'/'+h+'/'+i+'/'+i+'_ins_comp.jpg'
                if pst==True:
                    pign='/hetghome/jordan/Cose_nu/research_experience_program/COSEnu/pig_flat_rap_yvymys_N1_mu_100_Lmu_300_ar_0.7_pb.jpg'
                else:
                    pign='/hetghome/jordan/Cose_nu/research_experience_program/COSEnu/pig_flat_rap_yvymys_N1_mu_100_Lmu_300_ar_0.7_pb_pst.jpg'
                if os.path.exists(pign):
                    print('pig exit conti')
                    #continue
                ##################
                if 'pb' not in h:
                    GR = []
                    INDA = np.load(infn)
                    for j in INDA:
                        xm = 10**4*m.exp(-0.3*j[0])
                        if abs(xm-float(mu))<0.1:
                            GR.append(j[2])
                    mgr = max(GR)
                    print('mgr',mgr)
                    ####
                    GR = []
                    #INDA = np.load(infn)
                    for j in INDA:
                        xm = 10**4*m.exp(-0.3*j[0])
                        if abs(xm-float(mu))<0.1 and abs(j[1])<kmax:
                            GR.append(j[2])
                    mgrk = max(GR)
                    print('mgrk',mgrk)
                ###################################################
                Para = i.split("_")
                print(Para)
                Nz = float(Para[0])
                Nvz = float(Para[1])
                CFL = float(Para[2])
                
                #print(path+'/'+h+'/'+i)
                #print(i,'i')
                #DATA=np.fromfile(path+'/'+h+'/'+i+'/'+i+'_survival_probability.dat')
                #print(DATA)
                fn = path+'/'+h+'/'+i+'/'+i+'_survival_probability.dat'
                print(fn)
                if os.path.isfile(fn)==False:
                    print(fn)
                    print('no csv conti')
                    continue
                DATA = pd.read_csv(fn, sep='\t',names=["time","Pee","Pbee"])
                #DATA = pd.read_csv(path+'/'+i+'/'+i+'_conserved_quantities.dat')
                DATA = DATA.drop(index=0)
                toflo = {'time':float}
                DATA = DATA.astype(toflo)
                #pd.to_numeric(DATA["time"])
                """for j in DATA["time"]:
                    print(j)
                    print(type(j))"""
                DATA["physical time"]=DATA["time"].multiply(CFL*dz)
                DATA["1-Pee"]=1-DATA["Pee"]
                DATA["1-Pbee"]=1-DATA["Pbee"]
                if len(DATA.index[DATA['1-Pee']>0].tolist()) == 0:
                    plt.plot([1,2,3],[2,3,4])
                    plt.title('this system is stable so no fig')
                    plt.savefig(pign)
                    continue
                mint = DATA.at[min(DATA.index[DATA['1-Pee']>0].tolist()),"physical time"]
                DATA["aligned physical time"]=DATA["physical time"]-float(mint)
                ming = DATA.at[min(DATA.index[DATA['1-Pee']>0].tolist()),"1-Pee"]
                print('ming',ming)
                N=max(DATA.index[DATA["physical time"]<=xli].tolist())#,"physical time"#                                 [DATA[DATA["physical time"]<xli]]
                maxt= DATA.at[(DATA['1-Pee'].iloc[:N]).idxmax(),"physical time"]#DATA.at[DATA.idxmax(axis=0),"physical time"]
                print('maxt',maxt)
                print('mint',mint)
                #plt.title("dz,dt")
                #print(DATA.head(),'DATA head')
                #print(DATA.columns,'DATA[[0]]')
                print(DATA)
                if pst==True:
                    if nax!=0:
                        DATA.plot(x="physical time",y=["1-Pee"],label=[str(CFL)+r',$\mathrm{-log}_{10}\epsilon=$'+pb],ax=ax,c='C'+str(nax))
                    else:
                        ax=DATA.plot(x="physical time",y=["1-Pee"],label=[str(CFL)+r',$\mathrm{-log}_{10}\epsilon=$'+pb],c='C'+str(nax))
                else:
                    if nax!=0:
                        DATA.plot(x="aligned physical time",y=["1-Pee"],label=[str(CFL)+r',$\mathrm{-log}_{10}\epsilon=$'+pb],ax=ax,c='C'+str(nax))
                    else:
                        ax=DATA.plot(x="aligned physical time",y=["1-Pee"],label=[str(CFL)+r',$\mathrm{-log}_{10}\epsilon=$'+pb],c='C'+str(nax))
                nax+=1
                
                #plt.yscale('symlog')
                
                #plt.ylim(10**-8,0.6)
                
                #plt.savefig(pign)
                #cc+=1
                #plt.show()
                #sys.exit()
                #DATA.plot(x="time",y="Pbee")
                """if 'tur' in h:
                    plt.title("Averaged survival probability (tur) "+i)
                if 'non' in h:
                    plt.title("Averaged survival probability (w/o tur) "+i)"""
                #DATA.plot()
    plt.yscale('log')
    plt.xlim(0,xli)
    #plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.legend()
    plt.savefig(pign)

                #plt.close('all')
                #plt.show()
                #plt.plot()