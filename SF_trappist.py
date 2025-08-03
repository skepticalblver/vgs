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
#from scipy import stats
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
    if r < 0.03:
        res = 1e-05
    elif 0.03 < r < 0.05:
        res = 3e-3
    else:
        res = 1e-2 #5e-4
    return res


def carbon_mf(r):   #Hirschman et al. 2016
    if r < 0.035:
        res = 1e-6   #Marty et 2012
    elif 0.03 < r < 0.05:
        res = 2e-4  
    else:
        res = 1e-2    #Thomas etl. 1993,Kerridge et 1985
    return res


def nitrogen_mf(r):
    if r < 0.035:
        res = 1e-7
    elif 0.035 < r < 0.05:
        res = 1e-5
    else:
        res = 5e-4   #Fegley * Schaefer
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
