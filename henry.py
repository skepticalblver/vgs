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


def kH_t(T,gas):   #in mol/L*atm 
	T_ref = 298   #standard state temp
	if gas == 'nitrogen':
		kH_ref = 6.1e-4  
		C = 1300
	elif gas == 'carbon':
		kH_ref = 3.4e-2
		C = 2400

	elif gas == 'water':   #hydrogen
		kH_ref = 7.8e-4
		C = 500
	elif gas == 'helium':
		kH_ref = 3.7e-4
		C = 230
	elif gas == 'neon':
		kH_ref = 4.5e-4
		C = 490		

	elif gas ==  'argon':
		kH_ref = 1.4e-3
		C = 1300

	elif gas  == 'xenon':
		kH_ref = 4.3e-3 
		C = 2200

	kH_t = kH_ref*math.exp(C*(1/T - 1/T_ref))
	kH_t_cgs = kH_t/3e3  #mol/L*atm to mol/g*bar for rocky melt
	return kH_t_cgs


def kH_tf(T,gas):   #in mol/L*atm 
	T_ref = 298   #standard state temp
	if gas == 'nitrogen':
		kH_ref = 6.1e-4  
		C = 1300
	elif gas == 'carbon':
		kH_ref = 3.4e-2
		C = 2400

	elif gas == 'water':   #hydrogen
		kH_ref = 7.8e-4
		C = 500
	elif gas == 'helium':
		kH_ref = 3.7e-4
		C = 230
	elif gas == 'neon':
		kH_ref = 4.5e-4
		C = 490		

	elif gas ==  'argon':
		kH_ref = 1.4e-3
		C = 1300

	elif gas  == 'xenon':
		kH_ref = 4.3e-3 
		C = 2200

	kH_t = kH_ref
	kH_t_cgs = kH_t*3. /1e6 #mol/L*atm to mol*cm3/g*bar for "rocky" melt
	return kH_t_cgs



