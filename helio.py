#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure()

ax1 = plt.subplot(1, 1, 1)


atm_pres = 5.15e21 #grams
psurf_pres = 1.013 #bars

nitrogen_pres = atm_pres*0.767 #mass fraction*atm_pres
water_pres =  atm_pres*0.005
carbon_pres = atm_pres*0.004
argon_pres = atm_pres*0.01

ocean_pres = 1.4e24
ci_water = 0.15
ci_nitrogen =  0.1
ci_carbon = 0.05

au = [0,2,2.5,5]
water = [0.05,0.001,1e-5,1e-5]
water = water[::-1]
water = [x/ci_water for x in water]

nitrogen = [0.03,2.5e-3,5e-4,5e-4]
nitrogen = nitrogen[::-1]
nitrogen = [x/ci_nitrogen for x in nitrogen]

carbon = [0.02,0.01,0.001,0.001]
carbon = carbon[::-1]
carbon = [x/ci_carbon for x in carbon]

ax1.step(au,nitrogen,'r',lw='2',label='N$_2$')
ax1.step(au,water,'b',lw='2',label='H$_2$O')
ax1.step(au,carbon,'g',lw='2',label='CO$_2$')


ax1.legend(loc='best', frameon=False,fontsize=14)
#ax2.legend(loc='best', frameon=False,fontsize=14)
#ax1.set_ylim(10,1e6)
#ax1.set_yscale('log')
ax1.set_xlim(0,5)
ax1.set_ylim(0,0.45)
ax1.set_xlabel('heliocentric distance (AU)', fontsize=14)
#ax1.set_ylabel('surface pressure (bar) ', fontsize=14)
ax1.set_ylabel('Volatile Mass Fraction ($MF_{\\rm CI Chrondrite}$)', fontsize=14)
#fig.subplots_adjust(hspace=.01, wspace=.04, bottom=0.1, top=0.9, left=0.2, right=0.98)
#fig.set_size_inches(14, 7.5, forward=True)

plt.savefig('helio.pdf', dpi=200)                          

