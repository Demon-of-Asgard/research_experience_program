#!/hetghome/hetgsoft/anaconda3/bin/python3.9
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.fft import fft, fftfreq
import math as m
import math
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
#############################################################
from configparser import ConfigParser
########################################
path = os.getcwd()
print(path)
for h in os.listdir(path):
    if 'zlm' not in h:
        print('run123 conti')
        continue
    """if 'hetg' not in h:
        print('hetg pui conti')
        continue"""
    cc=0
    if 'jpg' in h:
        continue
    if '_non_mu_' in h:
        #print(h)
        mu=h[h.find('mu_')+3:h.find('_ar')]
        PT = []
        PEE = []
        Nz=0
        
        for i in os.listdir(path+'/'+h):
            if os.path.isfile(path+'/'+h+'/'+i):
                continue
            print(i)
            print('------------------')
            #dira = 
            if 'jpg' in i:
                continue
            """if '3000' not in i:
                print('3000 conti' )
                continue"""
            if os.path.isdir(path+'/'+h+'/'+i):
                if '48' in i:
                    print('48 conti')
                    continue
                print(path+'/'+h+'/'+i)
                pign = path+'/'+h+'/'+i+'/fft_off_'+i+'.jpg'
                if os.path.exists(pign)==True:
                    print('pig exi conti')
                    continue
                """if '/hetghome/jordan/Cose_nu/research_experience_program/COSEnu/fd_hetg_10_non_mu_100_ar_0.7_zlm_1_ini/9000_2_0.2' !=path+'/'+h+'/'+i:
                    print('test conti')
                    print(h)
                    print(i)
                    print('>>>>>>>>>>>>>>>>>>>>>>')
                    continue"""
                ###########################################
                confn = path+'/'+h+'/'+i+'/job.config'
                f=open(confn,'r')
                for line in f.readlines():
                    if "dz" in line:
                        dz=float(line[line.find(":")+1:])
                        #print(line)
                        #print(dz)
                #sys.exit()
                ############################################
                Para = i.split("_")
                print('Para i',Para)
                Nz = float(Para[0])
                Nvz = float(Para[1])
                CFL = float(Para[2])
                csvn = path+'/'+h+'/'+i+'/fft_off_'+i+'.csv'
                if os.path.exists(csvn)==False:
                    print('no csv conti')
                    continue
                DATA = pd.read_csv(csvn)
                #print(DATA)
                #print(DATA.columns[0])
                
                Time = list(map(float,DATA.columns[1:]))
                #print(Time)
                #sys.exit()
                QKS = []
                LBK = []
                ptnum = 6
                pti = int(DATA.shape[0]/ptnum)
                for j in range(ptnum):
                    print(pti*j)
                    QKS.append((DATA.loc[j*pti,:].values.flatten().tolist())[1:])
                    LBK.append(round(DATA.loc[j*pti,'k'],1))
                cc=0
                for j in range(ptnum):
                    plt.scatter(Time,QKS[j],s=0.9,c='C'+str(cc),label='k='+str(LBK[j]))
                    cc+=1
                
                plt.yscale('log')
                plt.legend(loc='lower right')
                #plt.show()
                #sys.exit()
                plt.savefig(pign,dpi=500)
                plt.close('all')
                #plt.show()
                #sys.exit()
                print('-------------------------------------------------------------------------------------------------')
                #Nz = 1
                #Nz = float(Para[0])
