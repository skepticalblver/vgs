#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import numpy as np 
import SF_flat
import SF_flat_rp
import henry
import random
import csv
import get_radius
import glob
import os

Mearth = 5.97e27  #in grams
Rearth = 6.3e8 # in cm
SAearth = 5.1e18 #in cm2
Vearth = 5.1e18*15000000  #in cm3
R_gas = 83.144   #(cm3*bar)/(K*mol)
cgrav  =  6.67e-8   # cm3/g/s2  gravitaitonal constant
g_accel_earth = 988 #cm/s^2

#fig, ax1 = plt.subplots(1, 1)
#read original nopl
au=[]
nindex=[]
ratio=[]
multiAU=[]
#read in initial object orbital distance file
with open( "aorig.dat", 'r') as f:  
	lines =  f.read().splitlines()

	for line in lines:
		
		nindex.append(line.split(None,2)[0])   
		au.append(line.split(None,2)[1])	
		multiAU.append([line.split(None,2)[0],line.split(None,2)[1]])	

final_indexL=[]
final_massL=[]
#read in final object mass and distance file
with open( "finalplanets.txt", 'r') as f:
        lines =  f.read().splitlines()

        for line in lines:

                final_indexL.append(line.split(None,6)[1])
#                final_massL.append(line.split(None,6)[2]*Mearth)

final_mass = 0.7068*Mearth

#read in asteroid file
ab_i = []
ab_radii = []
with open( "ABsizes.dat", 'r') as f:
        lines =  f.read().splitlines()

        for line in lines:

                ab_i.append(line.split(None,8)[0])
                ab_radii.append(line.split(None,8)[6])


ab_radii = [float(i) for i in ab_radii]


#read in asteroid file

#with open( "MODEL_OUT_emb6.csv", 'r') as f:
    #    lines =  f.read().splitlines()

   #     for line in lines:

  #              ab_i.append(line.split(None,8)[0])
 #               ab_radii.append(line.split(None,8)[6])


#with open( "MODEL_OUT_emb8.csv", 'r') as f:
    #    lines =  f.read().splitlines()

   #     for line in lines:

  #              ab_i.append(line.split(None,8)[0])
 #               ab_radii.append(line.split(None,8)[6])




r = 6000 #radius of Earth mantle in km
T_liquidus = -1.16e-7*r**3 + 0.0014*r**2 - 6.382*r + 1.444e4  #liquidus temperature

newdata=[]

i=0
Mcurr=0
Vesc= 11e5   #cm/s earth escape speed

Ts =1000 # initital surface temp in K
henry_water_save= 0.
henry_argon_save = 0.
henry_nitrogen_save= 0.
henry_carbon_save= 0.
mantle_grams = 0.
h_frac = 1.
Cp = 1e7 #specific heat of surface layer in ergs
ingass_grams = 0.
impactEsc= 0.
impactEsc_save = 0.
kappa_0  = 0.1  # absorption coefficient a the surface, in cm2/g, from Matsui & Abe 1986
Teq = 265  #Teq of Earth in K
sigma = 5.6704e-5  #Stefan-Boltzman's constant,  erg/cm2/s/K4
CompArray =  np.empty(5)

neon_grams = 0



CC_water = [1e-3]  #0.05, 1e-5
CC_nitrogen = [1e-6] #1e-4, 1e-7
CC_carbon =  [1e-5] #0.01, 1e-6
for z in range(len(CC_water)):
	#initialize list to store various atmophiles      
	neon_List=[]
	nitrogen_List=[]
	water_List=[]
	carbon_List=[]
	argon_List=[]
	atm_List=[]
	atm2_List=[]
	impactEsc_List=[]
	waterfrac_List=[]
	carbonfrac_List=[]
	nitrogenfrac_List=[]
	mantle_List=[]
	mantle2_List=[]
	esc_List=[]
	esc2_List=[]
	rpl_List=[]
	rpl2_List=[]
	mpl_List=[]
	atmpl_List=[]
	atmpl2_List=[]
	esctot_List=[]
	projectau_List=[]
	accret_List=[]
	accret2_List=[]
	mantle2_List=[]
	waterm_List=[]
	carbonm_List=[]
	ingassed_List=[]
	ingassed2_List=[]
	nitrogenm_List=[]
	argonm_List=[]
	ts_List=[]
	time_List=[]
	earthm_List=[]
	Psurf_List=[]
	mantle_List=[]
	fracEarth_List=[]
	fracEarth2_List=[]
	carbonc_List=[]



        j = open('SF_flat.py', 'r')
        g = j.read()
        j.close()
        g = g.replace("<<CC_water>>",str(CC_water[z]))
        g = g.replace("<<CC_nitrogen>>",str(CC_nitrogen[z]))
        g = g.replace("<<CC_carbon>>",str(CC_carbon[z]))
        h = open('SF_flat_rp.py', 'w')
        h.write(g)
        h.close()


        with open("planetgrowth.out",'r') as f:
		lines_pl =  f.read().splitlines()[1:]
		for line in lines_pl:  #read in entire file and separate by column
			time = float(line.split(None,8)[0])
			target_i = float(line.split(None,8)[1])
			target_m = float(line.split(None,8)[2])*Mearth
			project_i =  float(line.split(None,8)[3])
			project_m  = float(line.split(None,8)[4])*Mearth   #in grams
			impact_v  = float(line.split(None,8)[5])*1e5     #in cm/s
			impact_theta = (float(line.split(None,8)[7])/180.)*np.pi   #convert from degres to radians

			
			for i  in range(len(au)):  #assign each impactor and projector to an index from aorig
				if target_i == float(nindex[i]):
					target_au = float(au[i])
				if project_i == float(nindex[i]):
					project_au = float(au[i])

			if time  == 2.52380E+02:  #initialize values of each gas speices
				if target_i == 9.:
					print "hello0"
					target_r = get_radius.r_cm(target_m)
					g_accel = cgrav*target_m/(target_r)**2 
					neon_a = 1e18 #target_m*SF.neon_mf(target_au)*0.1
				
					water_a=1e18 # target_m*SF.water_mf(target_au)*0.005
					carbon_a=1e18 #target_m*SF.carbon_mf(target_au)*0.005
					argon_a=  1e18# target_m*SF.argon_mf(target_au)*0.005
					nitrogen_a= 1e18#target_m*SF.nitrogen_mf(target_au)*0.005
									
					atm_grams= water_a +  carbon_a + argon_a + nitrogen_a

					Psurf  = atm_grams*g_accel_earth/SAearth
					Psurf_bar =Psurf*1e-6				#mantle_grams = target_m*SF.water_mf(target_au) +target_m*SF.carbon_mf(target_au)
					
					water_m= target_m*float(SF_flat_rp.water_mf(target_au))
					carbon_m=  target_m*float(SF_flat_rp.carbon_mf(target_au))
					argon_m= target_m*float(SF_flat_rp.argon_mf(target_au))
					nitrogen_m= target_m*float(SF_flat_rp.nitrogen_mf(target_au))
					accret_tot = water_m+ carbon_m +argon_m +nitrogen_m

					mantle_grams =   water_m + carbon_m+ nitrogen_m + argon_m

					impactEsc = 0
					
					Psurf_List.append(round(Psurf_bar,3))
					waterfrac_List.append('{:0.3e}'.format(water_a))
					fracEarth_List.append(round(target_m/final_mass,3))
					nitrogenfrac_List.append('{:0.3e}'.format(nitrogen_a))
					carbonfrac_List.append('{:0.3e}'.format(carbon_a))
					argon_List.append('{:0.3e}'.format(argon_a))
					atm_List.append('{:0.3e}'.format(atm_grams))
					mantle_List.append('{:0.3e}'.format(mantle_grams))
					time_List.append(round(float(time),3))
					esctot_List.append(round(float(impactEsc),3))
					waterm_List.append('{:0.3e}'.format(water_m))
					carbonm_List.append('{:0.3e}'.format(carbon_m))
					ingassed_List.append('{:0.3e}'.format(ingass_grams))
					nitrogenm_List.append('{:0.3e}'.format(nitrogen_m))
					argonm_List.append('{:0.3e}'.format(argon_m))
					ts_List.append('{:0.3e}'.format(Ts))
					accret_List.append('{:0.3e}'.format(accret_tot))
					print'escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2)
			#starts calculation here
			if target_i ==  9.:  #change target index to desired embryo

				target_r = get_radius.r_cm(target_m)  #use function from SA Jacobson
				g_accel = cgrav*target_m/(target_r)**2   #cm/s^2

				if int(project_i) <= 40: # for embryos
					print "embryo"

					##########################deposit/accret gas$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
								
					water_a= water_a + project_m*float(SF_flat_rp.water_mf(project_au))*0.1
					carbon_a= carbon_a + project_m*float(SF_flat_rp.carbon_mf(project_au))*0.1
					argon_a= argon_a + project_m*float(SF_flat_rp.argon_mf(project_au))*0.1
					nitrogen_a= nitrogen_a + project_m*float(SF_flat_rp.nitrogen_mf(project_au))*0.1
					accret_tot = Mpl*float(SF_flat_rp.water_mf(project_au)) + Mpl*float(SF_flat_rp.carbon_mf(project_au)) + Mpl*float(SF_flat_rp.argon_mf(project_au))+ Mpl*float(SF_flat_rp.nitrogen_mf(project_au))

					atm_grams= water_a +nitrogen_a + carbon_a +argon_a

					water_mol_a = water_a/18.
					carbon_mol_a = carbon_a/12.   # g/(g/mol)
					#argon_mol = argon_grams/40
					nitrogen_mol_a = nitrogen_a/14
					argon_mol_a = argon_a/40

					total_mol =  water_mol_a + carbon_mol_a + nitrogen_mol_a + argon_mol_a

					water_molfrac = water_mol_a/total_mol
					carbon_molfrac =  carbon_mol_a/total_mol
					nitrogen_molfrac = nitrogen_mol_a/total_mol
					argon_molfrac = argon_mol_a/total_mol

					#########################impact loss from Schitling et al. (2015)###################
					target_r = get_radius.r_cm(target_m)
					Vesc = math.sqrt((2*cgrav*target_m)/target_r)
					kai = (impact_v*project_m)/(Vesc*target_m)
					impact_frac = 0.4*kai + 1.8*kai**2  - 1.2*(kai)**3

					impactEsc = impact_frac*atm_grams

					print "impact_frac=", impact_frac		
					impactEsc_save =  impactEsc_save + impactEsc  #save atm loss res
					water_a  = water_a -impactEsc*water_molfrac 
					carbon_a = carbon_a - impactEsc*carbon_molfrac 
					nitrogen_a = nitrogen_a - impactEsc*nitrogen_molfrac 
					argon_a = argon_a- impactEsc*argon_molfrac 

					atm_grams= water_a +nitrogen_a +carbon_a  +argon_a

                                        #########################hydrodynamics escape, energy-limited##############
                                        Feuv=29.7*(time/1e9)**(-1.23)*(1.23)**(-2)
                                        m_elim = (np.pi*Feuv*target_r**3)/(cgrav*target_m)

                                        water_a  = water_a - m_elim*water_molfrac
                                        carbon_a = carbon_a - m_elim*carbon_molfrac
                                        nitrogen_a = nitrogen_a - m_elim*nitrogen_molfrac
                                        argon_a = argon_a- m_elim*argon_molfrac

                                        atm_grams= water_a +nitrogen_a +carbon_a  +argon_a

					##########################Henry's Law#########
					Psurf  = atm_grams*g_accel/SA_planet
					Psurf_bar = Psurf*1e-6
						#Tsurf = (Psurf*Vearth)/(water_mol*R_gas)
					henry_water_eq= (henry.kH_t(Ta,'water')*Psurf_bar*water_molfrac)*water_mm*target_m*0.6
					henry_carbon_eq  = (henry.kH_t(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm*target_m*0.6
					henry_nitrogen_eq = (henry.kH_t(Ta,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm*target_m*0.6
					henry_argon_eq = (henry.kH_t(Ta,'argon')*Psurf_bar*argon_molfrac)*argon_mm*target_m*0.6

					henry_water_grams= float(henry_water_eq) - float(henry_water_save)
					henry_carbon_grams = float(henry_carbon_eq) - float(henry_carbon_save)
					henry_nitrogen_grams = float(henry_nitrogen_eq) - float(henry_nitrogen_save)
					henry_argon_grams = float(henry_argon_eq) - float(henry_argon_save)
					ingass_grams = henry_carbon_grams + henry_nitrogen_grams + henry_water_grams + henry_argon_grams	
					
						#print "water_mol=", (henry.kH_t(1000,'nitrogen')*Psurf_bar*nitrogen_molfrac)


					water_m = water_m +henry_water_grams
					carbon_m = carbon_m + henry_carbon_grams
					nitrogen_m = nitrogen_m + henry_nitrogen_grams
					argon_m = argon_m  + henry_argon_grams

					water_m = max(0,water_m)
					carbon_m = max(0,carbon_m)
					nitrogen_m = max(0,nitrogen_m)
					argon_m = max(0,argon_m)
					mantle_grams =   water_m + carbon_m+ nitrogen_m + argon_m


					henry_water_save = water_m
					henry_carbon_save = carbon_m
					henry_nitrogen_save = nitrogen_m
					henry_argon_save =  argon_m

						#print 'target_m2=',henry_water_eq
					if str(henry_water_eq) == 'nan':
						sys.exit()

					
					water_a  = water_a -  henry_water_grams -swind_loss
					carbon_a = carbon_a -  henry_carbon_grams
					nitrogen_a = nitrogen_a   -henry_nitrogen_grams
					argon_a = argon_a  -henry_argon_grams

						#print 'target_m_2=', water_grams,impactEsc*water_molfrac,henry_water_grams,mantle_grams
					water_a = max(0,water_a)
					carbon_a = max(0,carbon_a)
					nitrogen_a = max(0,nitrogen_a)
					argon_a = max(0,argon_a)
						

					atm_grams= water_a +nitrogen_a + carbon_a + argon_a
						
						#print 'target_m=',atm_grams,water_a,nitrogen_a,carbon_a , argon_a,Mcap
					Psurf  = atm_grams*g_accel/SA_planet
					Psurf_bar =Psurf*1e-6


					Psurf_List.append(round(Psurf_bar,3))  
					waterfrac_List.append('{:0.3e}'.format(water_a))
					fracEarth_List.append(round(target_m/final_mass,3))
					nitrogenfrac_List.append('{:0.3e}'.format(nitrogen_a))
					carbonfrac_List.append('{:0.3e}'.format(carbon_a))
					argon_List.append('{:0.3e}'.format(argon_a))
					atm_List.append('{:0.3e}'.format(atm_grams))
					mantle_List.append('{:0.3e}'.format(mantle_grams))
					time_List.append(round(float(time),3))
					esctot_List.append(round(float(impactEsc),3))
					waterm_List.append('{:0.3e}'.format(water_m))
					carbonm_List.append('{:0.3e}'.format(carbon_m))
					ingassed_List.append('{:0.3e}'.format(ingass_grams))
					nitrogenm_List.append('{:0.3e}'.format(nitrogen_m))
					argonm_List.append('{:0.3e}'.format(argon_m))
					ts_List.append('{:0.3e}'.format(Ts))
					accret_List.append('{:0.3e}'.format(accret_tot))

					print'escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2)

			
				elif int(project_i) > 40:   #for super-planetesiamls assembled from plantesiamls.

					print "planetesimal"

				
					target_r = get_radius.r_cm(target_m)  #use function from SA Jacobson
					g_accel = cgrav*target_m/(target_r)**2   #cm/s^2
					SA_planet = 4*np.pi*target_r**2
					while Mcurr <= project_m:
								
						Rpl =  random.choice(ab_radii)*1e5  #from km to cm
						Rhopl = 3.
						Mpl = Rhopl*(4./3*np.pi*Rpl**3)
				
					
						Mcurr  = Mcurr + Mpl

						##########################deposit/accret gas$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
						#neon_a= neon_a + Mpl*SF.neon_mf(project_au)
						#water_a= water_a + Mpl*SF.water_mf(project_au)
						#carbon_a= carbon_a + Mpl*SF.carbon_mf(project_au)
						#argon_a= argon_a + Mpl*SF.argon_mf(project_au)
						#nitrogen_a= nitrogen_a + Mpl*SF.nitrogen_mf(project_au)
						#atm_grams= water_a +nitrogen_a +carbon_a + argon_a
				
						
						water_m= water_m + Mpl*float(SF_flat_rp.water_mf(project_au))
						carbon_m= carbon_m + Mpl*float(SF_flat_rp.carbon_mf(project_au))
						argon_m= argon_m + Mpl*float(SF_flat_rp.argon_mf(project_au))
						nitrogen_m= nitrogen_m + Mpl*float(SF_flat_rp.nitrogen_mf(project_au))
						
						accret_tot = Mpl*float(SF_flat_rp.water_mf(project_au)) + Mpl*float(SF_flat_rp.carbon_mf(project_au)) + Mpl*float(SF_flat_rp.argon_mf(project_au))+ Mpl*float(SF_flat_rp.nitrogen_mf(project_au))
						#print 'target_m=',atm_grams,water_grams ,nitrogen_grams,carbon_grams , argon_grams
						mantle_grams =   water_m  + carbon_m  + argon_m + nitrogen_m

						water_mm = 18.  # in g /mol
						carbon_mm = 12.
						nitrogen_mm = 28.
						argon_mm = 40.				
					
						#print "water_mol=", henry_water_grams
						water_mol_a = water_a/18.
						
						carbon_mol_a = carbon_a/12.   # g/(g/mol)
			
						nitrogen_mol_a = nitrogen_a/28.
						argon_mol_a =  argon_a/40.
						
						
						total_mol_a =  water_mol_a + carbon_mol_a + nitrogen_mol_a +argon_mol_a
						if total_mol_a > 0:

							water_molfrac = water_mol_a/total_mol_a
							carbon_molfrac = carbon_mol_a/total_mol_a
							nitrogen_molfrac = nitrogen_mol_a/total_mol_a
							argon_molfrac =  argon_mol_a/total_mol_a
						Psurf  = atm_grams*g_accel/SA_planet
						

						#########################impact loss from Schitling et al. (2015)###################
						R_gas =  8.134e7  #units of cm3 barye/K/mol
						R_gas_specific = 287.06   # J/(kg*K)
						Ta= 1200        #base temp of atm, 1500
						rho_surf =(100000*Psurf_bar/(R_gas_specific*Ta))*0.001 #atmosphere base density. kg/m3 to g/cm3
						#rho_surf = 1.2*0.001
						target_r = get_radius.r_cm(target_m)  #use function from SA Jacobson
						g_accel = cgrav*target_m/(target_r)**2   #cm/s^2
						kerg = 1.38e-23  #J/K boltzman's constant
						
						
						mu =  (carbon_a/atm_grams)*carbon_mm+  (water_a/atm_grams)*water_mm+ (nitrogen_a/atm_grams)*nitrogen_mm+ (argon_a/atm_grams)*argon_mm       #calculate mean molecular weight 
						mu = max(1,mu)
						m_H = 1.67e-26 # hydrogen atomic mass in kg	
						H = 1e2*(kerg*Ta)/(mu*m_H*g_accel*0.01)    #convert from m to cm 
						
						H_earth = H
						Mcap = 2*np.pi*rho_surf*H_earth**2*target_r  
						

						r_min =  (3*rho_surf/Rhopl)**(1/3.)*H_earth  #1e5 for earth
						r_max = (3*(2*np.pi)**(1/2)*rho_surf/4.*Rhopl)**(1/3.)*(H_earth*target_r)**(1/2.) #25e5 for earth
						
						m_min = 4*np.pi*rho_surf*H
						m_max = math.sqrt(2)*rho_surf*(np.pi*H*target_r)**(3/2.)			
								
						if r_min < Rpl < r_max:
						#if m_min < Mpl < m_max:
							
							#impactEsc = Mcap*(math.fabs(math.log10(Rpl/r_min)))**3.5 #empiracal fit
							impactEsc = Mpl*((r_min/(2*Rpl))*(1- (r_min/Rpl)**2))
						elif Rpl < r_min:
						#elif Mpl < m_min:
							impactEsc  = 0.

						else:
							impactEsc = Mcap
						impactEsc = min(Mcap,impactEsc)

						impactEsc_save =  impactEsc_save + impactEsc  #save atm loss res
					#	impactEsc = Mpl*((r_min/2*Rpl)*(1-(r_min/Rpl)**2))
						rpl_List.append(float(Rpl))
						esc_List.append(float(impactEsc))
						mpl_List.append(float(Mpl))
						atmpl_List.append(float(atm_grams))
						###########################Calculate Henry's Law and Psurf#############
						water_a  = water_a -impactEsc*water_molfrac 
						carbon_a = carbon_a - impactEsc*carbon_molfrac 
						nitrogen_a = nitrogen_a - impactEsc*nitrogen_molfrac 
						argon_a = argon_a- impactEsc*argon_molfrac
					

						atm_grams= water_a +nitrogen_a +carbon_a  +argon_a
						atm_grams  = max(1,atm_grams)  #avoid zero atm values

                                        	#########################hydrodynamics escape, energy-limited##############
                                        	Feuv=29.7*(time/1e9)**(-1.23)*(1.23)**(-2)
                                        	m_elim = (np.pi*Feuv*target_r**3)/(cgrav*target_m)

                                        	water_a  = water_a - m_elim*water_molfrac
                                	        carbon_a = carbon_a - m_elim*carbon_molfrac
                        	                nitrogen_a = nitrogen_a - m_elim*nitrogen_molfrac
                	                        argon_a = argon_a- m_elim*argon_molfrac

        	                                atm_grams= water_a +nitrogen_a +carbon_a  +argon_a
	

						Psurf  = atm_grams*g_accel/SA_planet
						Psurf_bar = Psurf*1e-6
						#Tsurf = (Psurf*Vearth)/(water_mol*R_gas)
						
						henry_water_eq= (henry.kH_t(Ta,'water')*Psurf_bar*water_molfrac)*water_mm*target_m*0.6
						henry_carbon_eq  = (henry.kH_t(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm*target_m*0.6
						henry_nitrogen_eq = (henry.kH_t(Ta,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm*target_m*0.6
						henry_argon_eq = (henry.kH_t(Ta,'argon')*Psurf_bar*argon_molfrac)*argon_mm*target_m*0.6
						
						henry_water_grams= float(henry_water_eq) - float(henry_water_save)
						henry_carbon_grams = float(henry_carbon_eq) - float(henry_carbon_save)
						henry_nitrogen_grams = float(henry_nitrogen_eq) - float(henry_nitrogen_save)
						henry_argon_grams = float(henry_argon_eq) - float(henry_argon_save)
						
						ingass_grams = henry_carbon_grams + henry_nitrogen_grams + henry_water_grams + henry_argon_grams
						
						water_m = water_m +henry_water_grams
						carbon_m = carbon_m + henry_carbon_grams
						nitrogen_m = nitrogen_m + henry_nitrogen_grams
						argon_m = argon_m  + henry_argon_grams
						
						water_m = max(0,water_m)
						carbon_m = max(0,carbon_m)
						nitrogen_m = max(0,nitrogen_m)
						argon_m = max(0,argon_m)

						mantle_grams =   water_m + carbon_m + nitrogen_m+ argon_m #renormalize datra


						henry_water_save = water_m
						henry_carbon_save = carbon_m
						henry_nitrogen_save = nitrogen_m
						henry_argon_save =  argon_m
						
					
						if water_a/atm_grams >= 0.1:
				
							swind_loss = 0.8*water_a
						else:
							swind_loss  =0
						water_a  = water_a -  henry_water_grams #-swind_loss
						carbon_a = carbon_a -  henry_carbon_grams
						nitrogen_a = nitrogen_a   -henry_nitrogen_grams
						argon_a = argon_a  -henry_argon_grams
						water_a = max(0,water_a)
						carbon_a = max(0,carbon_a)
						nitrogen_a = max(0,nitrogen_a)
						argon_a = max(0,argon_a)

						atm_grams= water_a +nitrogen_a + carbon_a + argon_a
						atm_grams  = max(1,atm_grams)
					#	print 'target_m=',atm_grams,water_a,nitrogen_a,Mcap, impactEsc,henry_water_grams
						Psurf  = atm_grams*g_accel/SA_planet
						Psurf_bar =Psurf*1e-6
						
						
						#print F_atm, sigma,Teq,tau,Ts,atm_grams
						#if math.isnan(F_atm):
						#	sys.exit()
						if (Rpl > 1e7 or random.uniform(0,1) < 5e-4):
							fracEarth2_List.append(round(target_m/final_mass,3))
							rpl2_List.append(float(Rpl))
							esc2_List.append('{:0.3e}'.format(impactEsc))
							ingassed2_List.append('{:0.3e}'.format(ingass_grams))
							atm2_List.append('{:0.3e}'.format(atm_grams))
							mantle2_List.append('{:0.3e}'.format(mantle_grams))
							accret2_List.append('{:0.3e}'.format(accret_tot))	
					###################Append to various output lists###############################
					

					
					Psurf_List.append(round(Psurf_bar,3))
					waterfrac_List.append('{:0.3e}'.format(water_a))
					fracEarth_List.append(round(target_m/final_mass,3))
					nitrogenfrac_List.append('{:0.3e}'.format(nitrogen_a))
					carbonfrac_List.append('{:0.3e}'.format(carbon_a))
					argon_List.append('{:0.3e}'.format(argon_a))
					atm_List.append('{:0.3e}'.format(atm_grams))
					mantle_List.append('{:0.3e}'.format(mantle_grams))
					time_List.append(round(float(time),3))
					esctot_List.append('{:0.3e}'.format(impactEsc))
					waterm_List.append('{:0.3e}'.format(water_m))
					carbonm_List.append('{:0.3e}'.format(carbon_m))
					ingassed_List.append('{:0.3e}'.format(ingass_grams))
					nitrogenm_List.append('{:0.3e}'.format(nitrogen_m))
					argonm_List.append('{:0.3e}'.format(argon_m))
					ts_List.append('{:0.3e}'.format(Ts))
					accret_List.append('{:0.3e}'.format(accret_tot))
						#amantle_List.append(mantle_grams)
				      # print 'water_frac=', water_grams/atm_grams,'nitro_frac=',nitrogen_grams/atm_grams,'Psurf=', Psurf_bar,'mantle=', mantle_grams
					print'escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2), 'Ts = ', Ts

				Mcurr = 0
				mearth = target_m*Mearth
				earthm_List.append(mearth)
			#	time_List.append(float(time))

			#print('water_mol=',water_mol,'Psurf = ',Psurf,'atm_grams=', atm_grams,'mantle_grams=',mantle_grams)
			#withAU1 = line[:31] + str(target_au)  + line[31:]
			#withAU2 = withAU1[:-24] + str(project_au)  + withAU1[-24:]
			#newdata.append(withAU2)

	

	######################save data################
	
	with open('plcomp/OUT_ATM_ec_' + str(CC_water[z]), 'w') as f2:
		writer = csv.writer(f2, delimiter='\t')
		#writer.write(zip("Earth", "Psurf", "atm", "water", "carbon", "nitrogen", "argon", "mantle"))
		f2.write("%s %s %s %s %s %s %s %s %s" % ("Earth   ", "Time  ",  "Psurf   ", "atm   ", "water   ", "carbon   ", "nitrogen   ", "accretot", "esctot \n"))  
		writer.writerows(zip(fracEarth_List,time_List,Psurf_List,atm_List,waterfrac_List,carbonfrac_List,nitrogenfrac_List,accret_List,esctot_List))

	print'escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2), 'Ts = ', Ts

	with open('plcomp/OUT_MANT_ec_' + str(CC_water[z]), 'w') as f4:
		writer = csv.writer(f4, delimiter='\t')
 		f4.write("%s %s %s %s %s %s %s" % ("Earth", "Time", "Bulk mantle ", "Water  ", "carbon ", "Nitrogen ", "Argon \n"))
		writer.writerows(zip(fracEarth_List, time_List, mantle_List,waterm_List,carbonm_List, nitrogenm_List,argonm_List))

        with open('plcomp/atmdot_ec_'+ str(CC_water[z]), 'w') as f5:
                writer = csv.writer(f5, delimiter='\t')
                f5.write("%s %s %s %s %s %s %s" % ("Earth", "Radius", "Gas Escaped", "Ingassed", "Atm grams", "Mantle grams", "Accret tot\n"))
                writer.writerows(zip(fracEarth2_List, rpl2_List, esc2_List, ingassed2_List, atm2_List, mantle2_List,accret2_List))

                        
