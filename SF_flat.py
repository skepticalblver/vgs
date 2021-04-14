# 2017.10.24 18:31:50 CDT
#Embedded file name: /home/hcx2925/chem_pl/N-body_data/step_func.py
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
import glob
from math import log10
from matplotlib.pyplot import cm
import math
import random
from scipy import stats
import matplotlib as mpl
#fig, ax1 = plt.subplots(1, 1)
bins = np.linspace(4, 22, 20)
au = []
nindex = []
ratio = []
multiAU = []
with open('aorig.dat', 'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        nindex.append(line.split(None, 2)[0])
        au.append(line.split(None, 2)[1])
        multiAU.append([line.split(None, 2)[0], line.split(None, 2)[1]])

def water_mf(r):   #Raymond et al. 2009, Abe et al. 2000
    if r < 2:
        res =  "<<CC_water>>"
    elif 2 < r < 2.5:
        res =  "<<CC_water>>"
    else:
        res = "<<CC_water>>" #0.05 
    return res


def argon_mf(r):
    if r < 2.:
        res = 1e-14  
    else:
        res = 1e-10
    return res


def carbon_mf(r):
    if r < 2.5:
        res =  "<<CC_carbon>>"   #Marty et 2012
    elif 2 < r < 2.5:
        res =  "<<CC_carbon>>"  
    else:
        res = "<<CC_carbon>>" #0.02    #Thomas etl. 1993,Kerridge et 1985
    return res


def nitrogen_mf(r):
    if r < 2.0:
        res = "<<CC_nitrogen>>"
    elif 2 < r < 2.5:
        res =  "<<CC_nitrogen>>"
    else:
        res = "<<CC_nitrogen>>"     #0.003, Fegley and Shcaefer, Kerridge 1985
    return res



def neon_mf(r):
    if r < 2:
        res = 1e-11
    elif 2 < r < 2.5:
        res = 3e-11
    else:
        res = 4e-10
    return res


def xenon_mf(r):
    if r < 2:
        res = 1.302e-15
    elif 2 < r < 2.5:
        res = 6.55e-12
    else:
        res = 1.31e-11
    return res


def helium_mf(r):
    if r < 2.5:
        res = 5e-14
    else:
        res = 1e-11
    return res

au = sorted(au)
wmfList = []
for i in range(len(au)):
    wmfList.append(helium_mf(float(au[i])))

#ax1.plot(au, wmfList, lw='2', c='g')
#ax1.annotate('data from Marty et al. (2008)\nand Marty (2012)', xy=(4, 1e-12), fontsize=14)
#ax1.set_xlabel('orbital semi-major axis (AU)', fontsize=14)
#ax1.set_ylabel('He mass fraction', fontsize=14)
#fig.subplots_adjust(hspace=0.01, wspace=0.04, bottom=0.1, top=0.9, left=0.2, right=0.98)
#ax1.set_yscale('log')

# decompiled 1 files: 1 okay, 0 failed, 0 verify failed
# 2017.10.24 18:31:51 CDT
