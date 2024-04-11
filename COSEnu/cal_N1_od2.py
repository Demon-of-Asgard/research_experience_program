#!/hetghome/hetgsoft/anaconda3/bin/python3.9
#sed -i 's/\r//g' RK45_wo_tur_kmin.py
#chmod 777 RK45_wo_tur_kmin.py
from cProfile import label
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
import cmath 
acu = 1e-13
ts = 0.01
####################
#theta = m.pi/2
ar = 0.7
omg = 0.01
omga = -omg
v=1
mu=1
t1=600
t0=0
##############################################################
"""vxp=m.sin(theta)*v
vxn=m.sin(theta)*(-v)
vz=m.cos(theta)*v"""
#################################################
def absq(x):
    return cmath.sqrt(1-abs(x)**2)
def tbif(t,y):
    dqdt=np.array([0+0j,0+0j,0+0j,0+0j])
    #dqdt=y
    #print('0-0-0-0-0-0-0-0-0-0',dqdt[0],dqdt[1],dqdt[2],dqdt[3])
    #print(type(dqdt[0]))
    dqdt[0]=-omg*y[0]+0.75*0.5*mu*(absq(y[1])*y[0]-y[1]*absq(y[0]))-0.75*0.5*ar*mu*(absq(y[3])*y[0]-y[3]*absq(y[0]))
    dqdt[1]=-omg*y[1]+0.75*0.5*mu*(absq(y[0])*y[1]-y[0]*absq(y[1]))-0.75*0.5*ar*mu*(absq(y[2])*y[1]-y[2]*absq(y[1]))
    dqdt[2]=-omga*y[2]+0.75*0.5*mu*(absq(y[1])*y[2]-y[1]*absq(y[2]))-0.75*0.5*ar*mu*(absq(y[3])*y[2]-y[3]*absq(y[2]))
    dqdt[3]=-omga*y[3]+0.75*0.5*mu*(absq(y[0])*y[3]-y[0]*absq(y[3]))-0.75*0.5*ar*mu*(absq(y[2])*y[3]-y[2]*absq(y[3]))
    """dqdt[0]=-omg*y[0]+0.75*0.5*mu*(y[0]-y[1])-0.75*0.5*ar*mu*(y[0]-y[3])
    dqdt[1]=-omg*y[1]+0.75*0.5*mu*(y[1]-y[0])-0.75*0.5*ar*mu*(y[1]-y[2])
    dqdt[2]=-omga*y[2]+0.75*0.5*mu*(y[2]-y[1])-0.75*0.5*ar*mu*(y[2]-y[3])
    dqdt[3]=-omga*y[3]+0.75*0.5*mu*(y[3]-y[0])-0.75*0.5*ar*mu*(y[3]-y[2])"""
    #dqdt[0]+=-omg*y[0]
    #dqdt[1]+=-omg*y[1]
    #dqdt[2]+=-omga*y[2]
    #dqdt[3]+=-omga*y[3]
    #print('---------',dqdt[0],dqdt[1],dqdt[2],dqdt[3])
    dqdt=(-1j)*dqdt
    #print('>>>>>>>>>',dqdt[0],dqdt[1],dqdt[2],dqdt[3])
    return dqdt

fn='./f_RK45_N1_omg_'+str(omg)+'_ar_'+str(ar)+'/'
if os.path.isdir(fn) == False:    
    os.mkdir(fn)
#Q0 = np.array([1.0+0j,1.0+0j,1.0+0j,1.0+0j])/2
Q0 = np.array([10**-9+0j,10**-8+0j,10**-9+0j,10**-8+0j])/2
print(tbif(0,Q0))
#sys.exit()
"""def tbif(t,y):
    #dqdt=np.array([0+0j,0+0j,0+0j,0+0j])
    dqdt=y
    dqdt[0]=y[0]
    dqdt[1]=y[1]
    return dqdt
Q0 = np.array([10**-2+0j,10**-2+0j])"""
Sol = RK45(tbif,t0,Q0,t1,max_step=ts, rtol=acu, atol=acu * 10**-3, vectorized=False, first_step=None)
T =  []
Qt = []
Qn = []
T.append(t0)
Qt.append(Q0)
Qn.append(np.linalg.norm(Q0))
count_z  = 0
fot = 'asy_'
fot = 'asy98_'
for k in range(10**12):
    Sol.step()
    #print(Sol.t,Sol.y)
    if count_z%1000 ==0:
        T.append(Sol.t)
        Qn.append((np.linalg.norm(Sol.y)))
        Qt.append(((Sol.y)))
        print(Sol.t,'append')
        print(Sol.y,'soly')
        np.save(fn+fot+'RK45_T',T)
        np.save(fn+fot+'RK45_Qt',Qt)
        np.save(fn+fot+'RK45_Qn',Qn)
    count_z +=1
    if Sol.status == 'finished':
        break

"""np.save(fn+fot+'RK45_T',T)
np.save(fn+fot+'RK45_Qt',Qt)
np.save(fn+fot+'RK45_Qn',Qn)
"""