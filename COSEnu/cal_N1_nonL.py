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
acu = 1e-16
ts = 0.01
####################
theta = m.pi/2
ar = 0.7
omg = 0.01
omga = -omg
v=1
mu=1
t1=600
t0=0
##############################################################
vxp=m.sin(theta)*v
vxn=m.sin(theta)*(-v)
vz=m.cos(theta)*v
#################################################
def absq(x):
    return (1-abs(x)**2)**0.5
def comu(a,b):
    return np.matmul(a,b)-np.matmul(b,a)
Pz=np.array([[1,0],[0,-1]])
eps=10**-8
epsp=10**-7
rpp0=np.array([[1/2*(absq(eps)),1/2*eps],[1/2*eps,-1/2*(absq(eps))]])
rpn0=np.array([[1/2*(absq(eps)),1/2*eps],[1/2*eps,-1/2*(absq(eps))]])
rnp0=ar * np.array([[1/2*(absq(epsp)),1/2*epsp],[1/2*epsp,-1/2*(absq(epsp))]])
rnn0=ar * np.array([[1/2*(absq(epsp)),1/2*epsp],[1/2*epsp,-1/2*(absq(epsp))]])
y0=np.squeeze(np.array([rpp0,rpn0,rnp0,rnp0],dtype=complex).reshape((1,16)))
def tbif(t,y):
    #y[:4]=y.reshape()
    y=np.array([y[0:4].reshape((2,2)),y[4:8].reshape((2,2)),y[8:12].reshape((2,2)),y[12:16].reshape((2,2))])
    #print(y)
    #print(y[0])
    s = np.zeros((2,2),dtype=complex)
    dqdt=np.array([s,s,s,s])
    dqdt[0]=(omg*comu(-Pz,y[0])+0.75*comu(y[1]-np.conjugate(y[3]),y[0]))
    dqdt[1]=(omg*comu(-Pz,y[1])+0.75*comu(y[0]-np.conjugate(y[2]),y[1]))
    dqdt[2]=(omga*comu(-Pz,y[2])+0.75*comu(y[1]-np.conjugate(y[3]),y[2]))
    dqdt[3]=(omga*comu(-Pz,y[3])+0.75*comu(y[0]-np.conjugate(y[2]),y[3]))
    #print('------------------------------')
    #print(dqdt)
    rtn=np.squeeze(dqdt.reshape((1,16))/1j)
    #print(rtn)
    return rtn
"""def tbif(t,y):
    #dqdt=np.array([0+0j,0+0j,0+0j,0+0j])
    dqdt=y
    dqdt[0]=omg*dqdt[0]+0.75*0.5*mu*(absq(dqdt[1])*dqdt[0]-dqdt[1]*absq(dqdt[0]))+0.75*0.5*ar*mu*(absq(dqdt[3])*dqdt[0]-dqdt[3]*absq(dqdt[0]))
    dqdt[1]=omg*dqdt[1]+0.75*0.5*mu*(absq(dqdt[0])*dqdt[1]-dqdt[0]*absq(dqdt[1]))+0.75*0.5*ar*mu*(absq(dqdt[2])*dqdt[1]-dqdt[2]*absq(dqdt[1]))
    dqdt[2]=omga*dqdt[2]+0.75*0.5*ar*mu*(absq(dqdt[1])*dqdt[2]-dqdt[1]*absq(dqdt[2]))+0.75*0.5*ar*ar*mu*(absq(dqdt[3])*dqdt[2]-dqdt[3]*absq(dqdt[2]))
    dqdt[3]=omga*dqdt[3]+0.75*0.5*ar*mu*(absq(dqdt[0])*dqdt[3]-dqdt[0]*absq(dqdt[3]))+0.75*0.5*ar*ar*mu*(absq(dqdt[2])*dqdt[3]-dqdt[2]*absq(dqdt[3]))
    dqdt=1/(1j)*dqdt
    return dqdt"""
#print(tbif(0,y0))
#sys.exit()
fn='./f_nonL_RK45_N1_omg_'+str(omg)+'_ar_'+str(ar)+'/'
if os.path.isdir(fn) == False:    
    os.mkdir(fn)
#Q0 = np.array([1.0+0j,1.0+0j,1.0+0j,1.0+0j])/2
#Q0 = np.array([10**-8+0j,10**-8+0j,10**-8+0j,10**-8+0j])/2
print(y0)
Sol = RK45(tbif,t0,y0,t1,max_step=ts, rtol=acu, atol=acu * 10**-3, vectorized=False, first_step=None)
T =  []
Qt = []
Qn = []
Qo = []
T.append(t0)
Qt.append(y0)
Qn.append(np.linalg.norm((rpp0)[0][0]))
Qo.append(np.linalg.norm((rpp0)[0][1]))
count_z  = 0
for k in range(10**12):
    Sol.step()
    if count_z%1000 ==0:
        T.append(Sol.t)
        Qn.append((np.linalg.norm(((Sol.y)[0]))))
        Qo.append((np.linalg.norm(((Sol.y)[1]))))
        Qt.append(((Sol.y)))
        print(Sol.t,'append')
    count_z +=1
    if Sol.status == 'finished':
        break
np.save(fn+'nonL_RK45_T',T)
np.save(fn+'nonL_RK45_Qt',Qt)
np.save(fn+'nonL_RK45_Qn',Qn)
np.save(fn+'nonL_RK45_Qo',Qo)
