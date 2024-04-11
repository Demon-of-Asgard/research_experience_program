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

#########################################
path = os.getcwd()
#print(path)
for h in os.listdir(path):
    print(h)
    if 'N1' not in h:
        print('N1 conti')
        continue
    if 'RK45' not in h:
        print('RK45 conti')
        continue
    if 'non'  in h:
        print('non conti')
        continue
    print(h)
    if 'f_' not in h:
        print('f_ conti')
        continue
    ome=h[h.find('omg_')+4:h.find('_ar_')]
    fot = 'fat8_'
    #fot = 'asy98_'
    #for i in os.listdir(path+'/'+h):
    if os.path.exists(path+'/'+h+'/'+fot+'RK45_T.npy') == False:
        print('no file')
        continue
    T = np.load(path+'/'+h+'/'+fot+'RK45_T.npy')
    Qn = np.load(path+'/'+h+'/'+fot+'RK45_Qn.npy')
    plt.plot(T,Qn)
    plt.title(r'$\omega$='+ome)
    #plt.xlim(-1,5)
    plt.yscale('log')
    plt.savefig(path+'/'+h+'/'+fot+'RK45.jpg')
    plt.close('all')
    #plt.show()
    #sys.exit()