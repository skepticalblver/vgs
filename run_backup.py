#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import numpy as np 
import SF
import henry
import random
import csv
import get_radius

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

#read in asteroid file
ab_i = []
ab_radii = []
with open( "ABsizes.dat", 'r') as f:
        lines =  f.read().splitlines()

        for line in lines:

                ab_i.append(line.split(None,8)[0])
                ab_radii.append(line.split(None,8)[6])


ab_radii = [float(i)/2 for i in ab_radii]


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

henry_water_save= 0.
henry_argon_save = 0.
henry_nitrogen_save= 0.
henry_carbon_save= 0.
mantle_grams = 0.

impactEsc_save = 0.
#initialize list to store various atmophiles
neon_List=[]
nitrogen_List=[]
water_List=[]
carbon_List=[]
argon_List=[]
atm_List=[]
impactEsc_List=[]
waterfrac_List=[]
carbonfrac_List=[]
nitrogenfrac_List=[]
mantle_List=[]
esc_List=[]
rpl_List=[]
mpl_List=[]
atmpl_List=[]
esctot_List=[]
projectau_List=[]
#reservous output lists
time_List=[]
earthm_List=[]
Psurf_List=[]
mantle_List=[]
fracEarth_List=[]

CompArray =  np.empty(5)

neon_grams = 0
with open("planetgrowth.out", 'r') as f:
	lines_pl =  f.read().splitlines()[1:]
	for line in lines_pl:  #read in entire file and separate by column
		time = float(line.split(None,8)[0])
		target_i = float(line.split(None,8)[1])
		target_m = float(line.split(None,8)[2])*Mearth
		project_i =  float(line.split(None,8)[3])
		project_m  = float(line.split(None,8)[4])*Mearth   #in grams
		impact_v  = float(line.split(None,8)[5])*1e5     #in cm/s
		impact_theta = (float(line.split(None,8)[7])/180.)*np.pi   #concert from degres to radians

		
		for i  in range(len(au)):  #assign each impactor and projector to an index from aorig
			if target_i == float(nindex[i]):
				target_au = float(au[i])
			if project_i == float(nindex[i]):
				project_au = float(au[i])

		if time  == 4.18068E+03:  #initialize values of each gas speices
			if target_i == 7.:
				print "hello0"
				neon_a = 1e10 #target_m*SF.neon_mf(target_au)*0.1
                        
                        	water_a=1e10 # target_m*SF.water_mf(target_au)*0.005
                       		carbon_a=1e10 #target_m*SF.carbon_mf(target_au)*0.005
                        	argon_a=  1e10# target_m*SF.argon_mf(target_au)*0.005
                	        nitrogen_a= 1e10#target_m*SF.nitrogen_mf(target_au)*0.005
								
        	                atm_grams= water_a +  carbon_a + argon_a + nitrogen_a


                                Psurf  = atm_grams*g_accel_earth/SAearth
                                Psurf_bar =Psurf*1e-6


				#mantle_grams = target_m*SF.water_mf(target_au) +target_m*SF.carbon_mf(target_au)
				neon_m= target_m*SF.neon_mf(target_au)
				water_m= target_m*SF.water_mf(target_au)
				carbon_m= target_m*SF.carbon_mf(target_au)
				argon_m=  target_m*SF.argon_mf(target_au)
				nitrogen_m= target_m*SF.nitrogen_mf(target_au)

				mantle_grams =   water_m + carbon_m+ nitrogen_m + argon_m

				impactEsc = 0
				
                                Psurf_List.append(round(Psurf_bar,3))   #convert from Ba to bar
                                waterfrac_List.append('{:0.3e}'.format(water_a))
                                fracEarth_List.append(round(target_m/Mearth,3))
                                nitrogenfrac_List.append('{:0.3e}'.format(nitrogen_a))
                                carbonfrac_List.append('{:0.3e}'.format(carbon_a))
                                argon_List.append('{:0.3e}'.format(argon_a))
                                atm_List.append('{:0.3e}'.format(atm_grams))
                                mantle_List.append('{:0.3e}'.format(mantle_grams))
                                time_List.append(round(float(time),3))
                                esctot_List.append(round(float(impactEsc),3))
		
				print'escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2)
		#starts calculation here
		if target_i ==  7.:  #change target index to desired embryo

		
			if int(project_i) <= 40: # for embryos
				print "embryo"

				##########################deposit/accret gas$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
				neon_a= neon_a + project_m*SF.neon_mf(project_au)
				water_a= water_a + project_m*SF.water_mf(project_au)
				carbon_a= carbon_a + project_m*SF.carbon_mf(project_au)
				argon_a= argon_a + project_m*SF.argon_mf(project_au)
				nitrogen_a= nitrogen_a + project_m*SF.nitrogen_mf(project_au)
				
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

				###########################Calculate Henry's Law and Psurf#############

				water_a  = water_a -impactEsc*water_molfrac 
				carbon_a = carbon_a - impactEsc*carbon_molfrac 
				nitrogen_a = nitrogen_a - impactEsc*nitrogen_molfrac 
				argon_a = argon_a- impactEsc*argon_molfrac 

				atm_grams= water_a +nitrogen_a +carbon_a  +argon_a


				Psurf  = atm_grams*g_accel_earth/SAearth
				Psurf_bar = Psurf*1e-6
					#Tsurf = (Psurf*Vearth)/(water_mol*R_gas)
				henry_water_eq= (henry.kH_t(1000,'water')*Psurf_bar*water_molfrac)*water_mm*target_m*0.5*1e-3
				henry_carbon_eq  = (henry.kH_t(1000,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm*target_m*0.6*1e-3
				henry_nitrogen_eq = (henry.kH_t(1000,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm*target_m*0.6*1e-3
				henry_argon_eq = (henry.kH_t(1000,'argon')*Psurf_bar*argon_molfrac)*argon_mm*target_m*0.6*1e-3

				henry_water_grams= float(henry_water_eq) - float(henry_water_save)
				henry_carbon_grams = float(henry_carbon_eq) - float(henry_carbon_save)
				henry_nitrogen_grams = float(henry_nitrogen_eq) - float(henry_nitrogen_save)
				henry_argon_grams = float(henry_argon_eq) - float(henry_argon_save)
					

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
				henry_nitrogen_save = carbon_m
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
				Psurf  = atm_grams*g_accel_earth/SAearth
				Psurf_bar =Psurf*1e-6

                                Psurf_List.append(round(Psurf_bar,3))  
                                waterfrac_List.append('{:0.3e}'.format(water_a))
                                fracEarth_List.append(round(target_m/Mearth,3))
                                nitrogenfrac_List.append('{:0.3e}'.format(nitrogen_a))
                                carbonfrac_List.append('{:0.3e}'.format(carbon_a))
                                argon_List.append('{:0.3e}'.format(argon_a))
                                atm_List.append('{:0.3e}'.format(atm_grams))
                                mantle_List.append('{:0.3e}'.format(mantle_grams))
                                time_List.append(round(float(time),3))
                                esctot_List.append(round(float(impactEsc),3))

				rpl_List.append(float(Rpl))
				esc_List.append(float(impactEsc))
				mpl_List.append(float(Mpl))
				atmpl_List.append(float(atm_grams))

				print'escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2)

		
			elif int(project_i) > 40:   #for super-planetesiamls assembled from plantesiamls.

				print "planetesimal"
		
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
			
                                        neon_m= neon_m + Mpl*SF.neon_mf(project_au)
                                        water_m= water_m + Mpl*SF.water_mf(project_au)*0.1
                                        carbon_m= carbon_m + Mpl*SF.carbon_mf(project_au)*0.05
                                        argon_m= argon_m + Mpl*SF.argon_mf(project_au)
                                        nitrogen_m= nitrogen_m + Mpl*SF.nitrogen_mf(project_au)

					#print 'target_m=',atm_grams,water_grams ,nitrogen_grams,carbon_grams , argon_grams
					mantle_grams =  neon_m  + water_m  + carbon_m  + argon_m + nitrogen_m

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
					Psurf  = atm_grams*g_accel_earth/SAearth


					#########################impact loss from Schitling et al. (2015)###################
					R_gas =  8.134e7  #units of cm3 barye/K/mol
					Ta= 1000         #base temp of atm
					rho_surf = Psurf/(R_gas*Ta) #atmosphere base density
					target_r = get_radius.r_cm(target_m)  #use function from SA Jacobson
					g_accel = cgrav*target_m/(target_r)**2   #cm/s^2
					kerg = 1.38e-23  #J/K boltzman's constant
					mu =  (carbon_a/atm_grams)*carbon_mm+  (water_a/atm_grams)*water_mm+ (nitrogen_a/atm_grams)*nitrogen_mm+ (argon_a/atm_grams)*argon_mm       #calculate mean molecular weight 
					m_H = 1.67e-26 # hydrogen atomic mass in kg	
					H = 1e2*(kerg*Ta)/(mu*m_H*g_accel*0.01)    #convert from m to cm 
					
					H_earth = H
					Mcap = 2*np.pi*rho_surf*H_earth**2*target_r  


					r_min =  (3*rho_surf/Rhopl)**(1/3.)*H_earth  #1e5 for earth
					r_max = (3*(2*np.pi)**(1/2)*rho_surf/4.*Rhopl)**(1/3.)*(H_earth*target_r)**(1/2.) #25e5 for earth
									
					if r_min < Rpl < r_max:
						Rpl_km = Rpl/1e5
						impactEsc = Mcap*(math.fabs(math.log10(Rpl_km)))**1.8  #empiracal fit
						#impactEsc = Mpl*((r_min/(2*Rpl))*(1-(r_min/Rpl)**2))
					elif Rpl < r_min:
						impactEsc  = 0.

					else:
						impactEsc = Mcap#Mpl*random.uniform(1e-4,0.1)
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

					Psurf  = atm_grams*g_accel_earth/SAearth
					Psurf_bar = Psurf*1e-6
					#Tsurf = (Psurf*Vearth)/(water_mol*R_gas)
					henry_water_eq= (henry.kH_t(1000,'water')*Psurf_bar*water_molfrac)*water_mm*target_m*0.5*1e-3
					henry_carbon_eq  = (henry.kH_t(1000,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm*target_m*0.6*1e-3
					henry_nitrogen_eq = (henry.kH_t(1000,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm*target_m*0.6*1e-3
					henry_argon_eq = (henry.kH_t(1000,'argon')*Psurf_bar*argon_molfrac)*argon_mm*target_m*0.6*1e-3

					henry_water_grams= float(henry_water_eq) - float(henry_water_save)
					henry_carbon_grams = float(henry_carbon_eq) - float(henry_carbon_save)
					henry_nitrogen_grams = float(henry_nitrogen_eq) - float(henry_nitrogen_save)
					henry_argon_grams = float(henry_argon_eq) - float(henry_argon_save)
					

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
					henry_nitrogen_save = carbon_m
					henry_argon_save =  argon_m

					#print 'target_m2=',henry_water_eq
					if str(henry_water_eq) == 'nan':
						sys.exit()

				
					if water_a/atm_grams >= 0.1:
			
						swind_loss = 0.7*water_a
					else:
						swind_loss  =0
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
					Psurf  = atm_grams*g_accel_earth/SAearth
					Psurf_bar =Psurf*1e-6
					#print'escape=','{:0.3e}'.format(water_m),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2)
					#sys.exit()


				###################Append to various output lists###############################

				Psurf_List.append(round(Psurf_bar,3))   #convert from Ba to bar
				waterfrac_List.append('{:0.3e}'.format(water_a))
				fracEarth_List.append(round(target_m/Mearth,3))
				nitrogenfrac_List.append('{:0.3e}'.format(nitrogen_a))
				carbonfrac_List.append('{:0.3e}'.format(carbon_a))
				argon_List.append('{:0.3e}'.format(argon_a))
				atm_List.append('{:0.3e}'.format(atm_grams))
				mantle_List.append('{:0.3e}'.format(mantle_grams))
				time_List.append(round(float(time),3))
				esctot_List.append(round(float(impactEsc),3))

					#amantle_List.append(mantle_grams)
				#nitrogen_List.append(nitrogen_grams)
				#water_List.append(water_grams)
				#carbon_List.append(carbon_grams)
				#argon_List.append(argon_grams)
						
				#atm_List.append(atm_grams)

				
				#fracEarth_List.append(target_m/Mearth)

				#print('escape=',impactEsc,'Psurf = ', Psurf,'atm_grams=', atm_grams)

                               # print 'water_frac=', water_grams/atm_grams,'nitro_frac=',nitrogen_grams/atm_grams,'Psurf=', Psurf_bar,'mantle=', mantle_grams
				print'escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2), 'Scale height = ', round(H,2)

			Mcurr = 0
			mearth = target_m*Mearth
			earthm_List.append(mearth)
		#	time_List.append(float(time))

		#print('water_mol=',water_mol,'Psurf = ',Psurf,'atm_grams=', atm_grams,'mantle_grams=',mantle_grams)
		#withAU1 = line[:31] + str(target_au)  + line[31:]
		#withAU2 = withAU1[:-24] + str(project_au)  + withAU1[-24:]
		#newdata.append(withAU2)

f.close()

######################save data################

with open('MODEL_OUT_VENUS.csv', 'w') as f:
	writer = csv.writer(f, delimiter='\t')
	#writer.write(zip("Earth", "Psurf", "atm", "water", "carbon", "nitrogen", "argon", "mantle"))
	f.write("%s %s %s %s %s %s %s %s" % ("Earth   ", "Psurf   ", "atm   ", "water   ", "carbon   ", "nitrogen   ", "argon   ", "mantle   \n"))  
	writer.writerows(zip(fracEarth_List,Psurf_List,atm_List,waterfrac_List,carbonfrac_List,nitrogenfrac_List,argon_List,mantle_List,esctot_List))


#with open('pl_dist.csv', 'w') as f2:
#	writer = csv.writer(f2, delimiter='\t')
#	writer.writerows(zip(rpl_List,mpl_List,esc_List,atmpl_List))

