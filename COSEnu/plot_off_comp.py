#!/hetghome/hetgsoft/anaconda3/bin/python3.9
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
    if 'run3' not in h:
        print('run123 conti')
        continue
    if 'hetg' not in h:
        print('hetg pui conti')
        continue
    if 'N1' not in h:
        print('N1  conti')
        continue
    cc=0
    if 'jpg' in h:
        continue
    if '_mu' in h:
        #print(h)
        mu=h[h.find('mu_')+3:h.find('_Lmu')]
        PT = []
        PEE = []
        Nz=0
        pign = '1234'
        for i in os.listdir(path+'/'+h):
            if os.path.isfile(path+'/'+h+'/'+i):
                continue
            print(i)
            #dira = 
            if 'jpg' in i:
                continue
            if os.path.isdir(path+'/'+h+'/'+i):
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
                for k in os.listdir(path+'/'+h):
                    if os.path.isfile(path+'/'+h+'/'+k):
                        continue
                    Para = k.split("_")
                    print('Para k',Para)
                    if 'jpg' in k:
                        continue
                    if Nz==0:
                        print('nz 0')
                    else:
                        """if Nz==1:
                            Nz = float(Para[0])"""
                        if Nz != float(Para[0]):
                            print(Nz,'pre Nz')
                            print(float(Para[0]))
                            print(' not same')
                            #Nz = float(Para[0])
                            continue   
                        else:
                            print(Nz,'pre Nz')
                            print(float(Para[0]))
                            print('same')
                    #print(i)
                    #dira = 
                    if os.path.isdir(path+'/'+h+'/'+k):
                        Nz = float(Para[0])
                        Nvz = float(Para[1])
                        CFL = float(Para[2])
                        #print(path+'/'+h+'/'+k)
                        if os.path.exists(pign):
                            print('pig exit conti')
                            #continue
                        #print(path+'/'+h+'/'+i)
                        #print(i,'i')
                        #DATA=np.fromfile(path+'/'+h+'/'+i+'/'+i+'_survival_probability.dat')
                        #print(DATA)
                        Time = []
                        AOp= [[],[],[],[]]
                        AOpb = []
                        AOn= []
                        AOnb = []
                        for l in os.listdir(path+'/'+h+'/'+k):
                            if 'dat' in l:
                                if '_rho' in l:
                                    print(l)
                                    fn = path+'/'+h+'/'+k+'/'+l
                                    itr=float(l[l.find("_rho_")+5:l.find(".dat")])
                                    DATA = pd.read_csv(fn, sep='\t',names=["Z", "vz", "rho_ee", "rho_xx", "Re[rho_ex]", "Im[rho_ex]", "brho_ee", "brho_xx", "Re[brho_ex]", "Im[brho_ex]"])
                                    DATA["sqr"]=DATA["Re[rho_ex]"] * DATA["Re[rho_ex]"] + DATA["Im[rho_ex]"] * DATA["Im[rho_ex]"] 
                                    DATA["sqrb"]=DATA["Re[brho_ex]"] * DATA["Re[brho_ex]"] + DATA["Im[brho_ex]"] * DATA["Im[brho_ex]"]
                                    #DATA = pd.read_csv(path+'/'+i+'/'+i+'_conserved_quantities.dat')
                                    #ave_off=(DATA["sqr"].mean())**0.5
                                    #print(DATA)
                                    for m in range(int(Nvz)):
                                        if m ==1:
                                            AOp[2].append((DATA.loc[m,"sqr"])**0.5)#+np.sin(itr * CFL*dz*100)*((DATA.loc[m,"sqr"])**0.5*0.01))
                                            AOp[3].append((DATA.loc[m,"sqrb"])**0.5)
                                        if m==0:
                                            AOp[0].append((DATA.loc[m,"sqr"])**0.5)
                                            AOp[1].append((DATA.loc[m,"sqrb"])**0.5)
                                    #print(ave_off)
                                    Time.append(itr * CFL*dz)
                                    #AOp.append(ave_off)
                                    #sys.exit()
                                    """DATA = DATA.drop(index=0)
                                    toflo = {'time':float}
                                    DATA = DATA.astype(toflo)"""
                        sTime = sorted(Time)
                        sTimeid=sorted(range(len(Time)), key=lambda l: Time[l])
                        print(len(AOp[2]),'AOp[2]')
                        print(len(sTimeid),'sTimeid')
                        print(h,k)
                        lab = [r'$|\rho^{ex}_{v_{-}}|$',r'$|\bar{\rho}^{ex}_{v_{-}}|$',r'$|\rho^{ex}_{v_{+}}|$',r'$|\bar{\rho}^{ex}_{v_{+}}|$']
                        sAOp= [[],[],[],[]]
                        
                        for l in range(len(AOp)):
                            for m in sTimeid:
                                #print(l,m)
                                sAOp[l].append((AOp[l])[m])
                            if l<2:
                                plt.plot(sTime,sAOp[l],linestyle='dashed',linewidth=2.5,label=lab[l],c='C'+str(l))
                            else:
                                plt.plot(sTime,sAOp[l],label=lab[l],c='C'+str(l))
                        """plt.scatter(Time,AOpb,s=0.7,label=,c='C'+str(1))
                        plt.scatter(Time,AOn,s=0.7,label=,c='C'+str(2))
                        plt.scatter(Time,AOnb,s=0.7,label=,c='C'+str(3))"""
                        
                        plt.yscale('log')
                        #plt.xlim(150,151)
                        #plt.ylim(0.005,0.01)
                        plt.legend()
                        #plt.show()
                        pign = path+'/'+h+'/fig_'+k+'_each_off_compare.jpg'
                        plt.savefig(pign)
                        plt.close('all')
                        ##sys.exit()
                        #cc+=1
                        
                        #DATA.plot(x="time",y="Pbee")
                        """if 'tur' in h:
                            plt.title("Averaged survival probability (tur) "+i)
                        if 'non' in h:
                            plt.title("Averaged survival probability (w/o tur) "+i)"""
                        #DATA.plot()
                        #
                        #plt.show()
                        #plt.plot()
                        print('cc',cc)
                    """if cc>2:
                        print('break')
                        break"""
                """pign = path+'/'+h+'/'+h+'_dz_'+str(dz)[:8]+'_off_compaare.jpg'
                print('pign', pign)
                plt.title('$\omega=$'+(str(1/float(mu)))[:5])
                plt.legend()
                #plt.ylabel('')
                #plt.legend()
                plt.yscale('log')
                print('show')
                #plt.show()
                #sys.exit()
                plt.savefig(pign)
                plt.close('all')"""
                print('-------------------------------------------------------------------------------------------------')
                #Nz = 1
                #Nz = float(Para[0])
