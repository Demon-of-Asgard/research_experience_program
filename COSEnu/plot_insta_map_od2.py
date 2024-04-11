#!/hetghome/hetgsoft/anaconda3/bin/python3.9
#sed -i 's/\r//g' plot_insta_map_km_prd.py
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
mude = 0.3
sml = 4
omg=1
mu0=10000
"""for h in range(2):
    h=1
    if h ==0:
        small = True
    else:
        small = False"""
pth = '/hetghome/jordan/Cose_nu/research_experience_program/COSEnu/f_insta_map/'
for hi in os.listdir(pth):
    if '_f' in hi:
        continue
    if '_fpg' in hi:
        continue
    print('hi',hi)
    """if hi != 'Lk_-6_Hk_6_Lz_0_Hz_15_Lmu_300_ar_0.7_deg_120.0.npy':
        print('hi conti')
        continue"""
    if 'png' not in hi:
        if 'jpg' not in hi:
            dtn = hi
            dtn = "Lk_0_Hk_0_lmu_10_hmu_150_Lmu_30_ar_0.7_deg_90.0_gz_od2.npy"
            #dtn= 'Lk_-6_Hk_6_Lz_0_Hz_15_Lmu_300_ar_0.7_deg_60.0.npy'
            if 'ar' in hi:
                Lmu = (dtn[dtn.find('_Lmu_')+5:dtn.find('_ar_')])
            else:
                Lmu = (dtn[dtn.find('_Lmu_')+5:dtn.find('.npy')])
            ar = (dtn[dtn.find('_ar_')+4:dtn.find('_deg_')])
            deg = (dtn[dtn.find('_deg_')+5:dtn.find('.npy')])
            MUT = [10000,8500,7000,6000,5000,4000]
            MUT = [200,175,150,125,100,75]
            MUT = [130,120,110,100,90,80,70]
            km = -150#/float(mu)
            sfn = "./f_insta_map/insta_map_bi_Lmu_"+Lmu+"_ar_"+ar+"_deg_"+deg+"_rgmuk_"+str(max(MUT))+"_"+str(min(MUT))+"_"+str(abs(km))+"_tick_od2.jpg"
            """print((sfn))   
            print(exists(sfn))
            print((nfn))
            print(exists(nfn))"""
            if exists(sfn)== True:
                print(sfn)
                print('pig exist')
                continue
            #if (exists(sfn) == False) and (exists(nfn) == False):
            apa = float(ar)
            if apa==1:
                print('apa conti')
                continue
            #dtn = "Lk_-6_Hk_6_Lz_0_Hz_15_Lmu_300_ar_0.7_deg_45.0.npy"
            DATA = np.load(pth+dtn)
            #DATA = np.load(pth+'/Lk_-6_Hk_6_Lz_0_Hz_15_Lmu_300_ar_0.7_deg_60.0.npy')
            #A = DATA[:50]+DATA[-50:]
            #print(A[:,[0]])
            #sys.exit()
            """X = DATA[:,[0]]
            Y = DATA[:,[1]]
            x,y = np.meshgrid(X,Y)"""
            cm = plt.cm.get_cmap('RdYlBu')
            #X, Y, Z = grid(i[0], i[1], i[2])
            #plt.contourf(X, Y, Z)
            XM = []
            X = []
            Y = []
            Z = []
            cou = 0
            jup=0
            XT = []
            
            
            for i in DATA:
                if jup%1 != 0:
                    jup+=1
                    continue
                xm = 10**4*m.exp(-0.3*i[0])
                if max(MUT)> xm and xm> min(MUT):
                    if i[1] > km and i[1] < -km:
                        XM.append(xm)
                        X.append(i[0])
                        Y.append(i[1])
                        Z.append(i[2])
                #print(i[0])
                if cou%5000000 == 0:
                    print(cou)
                cou +=1
                jup+=1
                """if cou%5000000 == 5000000-1:
                    break"""
            #fig = plt.figure()
            #ax1 = fig.add_subplot(111)
            fig, axs = plt.subplots(1, 1,layout='constrained')
            ax1=axs#[0,0]
            ####
            """fig, ax1 = plt.subplots()
            ax2 = ax1.twiny()
            ax2.set_xlabel()"""
            #sc = ax1.scatter(XM,Y,c = Z,cmap=cm, vmin=0, vmax=18,s=0.7)
            sc = ax1.scatter(XM,Y,c = Z,cmap=cm,s=0.7)
            #sc = ax1.scatter(X,Y,c = Z,cmap=cm,s=0.7)
            #divider = make_axes_locatable(plt.gca())
            #cax = divider.append_axes("right", "5%", pad="3%")
            #plt.colorbar(sc,cax=cax)
            #plt.savefig('/hetghome/jordan/NeuOsc/PAP_kmin/f_insta_map/test.jpg',dpi=400)
            #sys.exit()
            #cbar.set_label('growth rate')
            #ax1.set_xlabel("$z$[km]")
            ax1.set_ylabel("$k$[km$^{-1}$]")
            ax1.tick_params(left=True,right=True)
            LX = int(len(X)/5)
            """XT = []
            MUT = []
            for i in range(5):
                XT.append(X[i*LX])
                MUT.append(int(10000*m.exp(-mude*X[i*LX])))"""
            
            for i in range(len(MUT)):
                for j in range(len(X)):
                    if abs(10000*m.exp(-mude*X[j])-MUT[i])<1:
                        XT.append(XM[j])
                        print(MUT[i],'append')
                        break
                #print(MUT[i])
                #print('sth is wrong')
            
            #XT.append(int(XM[-1]))
            #MUT.append(int(XM[-1]))
            print('XT',XT)
            print('MUT',MUT)
            ####################################################################plot with the analysis solution approximation 
            if True == False:
                XA = []
                YA = []
                YA1 = []
                
                for i in range(10000):
                    z = i*sml/10000
                    XA.append(z)
                    #XA.append(z)
                    YA.append(2.25*2/3**0.5*(4*(1+apa)/(1-apa)**3)**(-1)/omg*(mu0*m.exp(-mude*z))**2)
                    YA1.append(-2.25*2/3**0.5*(4*(1+apa)/(1-apa)**3)**(-1)/omg*(mu0*m.exp(-mude*z))**2)
                #print(YA)
                #print(XA)
                plt.plot(XA,YA,c='black',linewidth=2)
                plt.plot(XA,YA1,c='black',linewidth=2)
            ######################################################################
            ax1.set_xlim(xmin=MUT[0], xmax=MUT[-1])
            #fig.colorbar(sc,cax=cax)
            #ax2 = ax1.twiny()
            #ax2.set_xlim(ax1.get_xlim())
            ax1.set_title(r'$ \theta $='+str(deg)+r', $\beta$'+'='+str(Lmu)+r', $\alpha$='+str(ar))
            #ax1.set_title(r'$ \theta $='+str(deg)+', $\u03B2$'+'='+str(Lmu)+', $\u03B1$'+str(ar))
                
            ax1.set_xticks(XT)
            ax1.set_xticklabels(MUT)
            ax1.set_xlabel("$\omega$")
            #ax2.spines['right'].set_visible(False)
            fig.colorbar(sc,shrink=1.0)
            #plt.xticks(XT,MUT)
            #plt.close('all')
            #plt.plot(XA,YA,linewidth=5)
            #plt.ylim(-600000,600000)
            plt.ylim(km,-km)
            #plt.show()
            plt.savefig(sfn,dpi=400)
            sys.exit()


            """for i in DATA:
                #print(i[2])
                if i[2] == 0:
                    B.append(i)
                elif i[2] <1:
                    G.append(i)
                elif i[2] <2:
                    Y.append(i)
                elif i[2] <3:
                    O.append(i)
                else:
                    R.append(i)
            cc = 0
            col = ['blue','green','yellow','darkorange','red']
            for i in [B,G,Y,O,R]:
                X = []
                #print(R)
                for j in i:
                    X.append(j[0])
                Y = []
                for j in i:
                    Y.append(j[1])
                if cc == 0:
                    plt.scatter(X,Y,s = 0.1,c = col[cc],label = 'eig-val='+str(cc))
                elif cc == 4:
                    plt.scatter(X,Y,s = 0.1,c = col[cc],label = 'eig-val>'+str(cc))
                else:
                    plt.scatter(X,Y,s = 0.1,c = col[cc],label = 'eig-val~'+str(cc))
                cc +=1"""
            """plt.scatter(R[:,[0]],[i[1]],c = 'blue', s = 0.1)
            plt.scatter([i[0]],[i[1]],c = 'green', s = 0.1)
            plt.scatter([i[0]],[i[1]],c = 'yellow', s = 0.1)
            plt.scatter([i[0]],[i[1]],c = 'darkorange', s = 0.1)
            plt.scatter([i[0]],[i[1]],c = 'red', s = 0.1)"""
            
            #plt.title("instability map of different mode k")
            #plt.xlabel()
            
            #plt.tight_layout(rect=(0,0,1,1))
            #plt.savefig('/hetghome/jordan/NeuOsc/PAP_kmin/f_insta_map/test.jpg')
            
            #plt.show()
            #sys.exit()
            #plt.legend(loc = 'upper right')
            #sfn = "./f_insta_map/insta_map_small_Lmu_"+Lmu+".jpg"
            if sfn == nfn:
                if small == True:
                    plt.savefig(sfn,dpi=400)
                else:
                    plt.savefig(sfn,dpi=400)
            else:
                if small == True:
                    plt.savefig(nfn,dpi=400)
                else:
                    plt.savefig(nfn,dpi=400)
            
            #plt.show()
            plt.close('all')