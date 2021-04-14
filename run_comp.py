#!/usr/bin/env python
import numpy as np
import math
import sys
import SF_variable
import SF_variable_rp    #step function parameterization for volatile frac
import henry  #solubility constants
import random
import csv
import get_radius  #func to calc embryo radius 
import glob
import run_embryo_f

Mearth = 5.97e27  #in grams
Rearth = 6.3e8 # in cm
SAearth = 5.1e18 #in cm2
Vearth = 5.1e18*15000000  #in cm3
R_gas = 83.144   #(cm3*bar)/(K*mol)
cgrav  =  6.67e-8   # cm3/g/s2  gravitaitonal constant
g_accel_earth = 988 #cm/s^2
e_m = 9e16  # specific energy, J kg

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
henry_nitrogen_save= 0.
henry_carbon_save= 0.
mantle_grams = 0.
h_frac = 1.
Cp = 1e7 #specific heat of surface layer in ergs
ingass_grams = 0.
impactEsc= 0.
impactEsc_save = 0.
D_c = 3000  #partition coefficient for carbon
kappa_0  = 0.1  # absorption coefficient a the surface, in cm2/g, from Matsui & Abe 1986
Teq = 265  #Teq of Earth in K
sigma = 5.6704e-5  #Stefan-Boltzman's constant,  erg/cm2/s/K4
project_rho = 3.

earth_list=[9,9,6,20,16,18,6,12,19,28,9,15,8,7,23,5]

emb_list=[40,40,31,69,69,90,90,90,58,110,110,110,110,71,71,71]

CompArray =  np.empty(5)


counter = 0
filepath = 'nbody_sims/1to1-0.5-8'
CC_water = [0.2,0.15,0.1,0.08,0.05,0]  #default = 0.05
CC_nitrogen = [1e-2,7e-3,5e-3,1e-3,5e-4,1e-4,0] #1e-3
CC_carbon =  [0.15,0.1,5e-2,1e-2,5e-3,0.01,0] #0.01
for z in range(0,len(CC_water)):
	print(CC_water[z])

	initial = True #if true, run initialization8to1-0.5-10
	#initialize list to store various atmophiles      
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
	carbonc_List=[]
	ingassed2_List=[]
	nitrogenm_List=[]
	ts_List=[]
	time_List=[]
	earthm_List=[]
	Psurf_List=[]
	mantle_List=[]
	fracEarth_List=[]
	fracEarth2_List=[]

	au=[] 
	nindex=[] 
	ratio=[] 
	multiAU=[] 
    
    #replace empty strings in SF_var with designited values from the list
	j = open('SF_variable.py', 'r')
	g = j.read()
	j.close()
	g = g.replace("<<CC_water>>",str(CC_water[z]))
	g = g.replace("<<CC_nitrogen>>",str(CC_nitrogen[5]))
	g = g.replace("<<CC_carbon>>",str(CC_carbon[5]))
	h = open('SF_variable_rp.py', 'w')
	h.write(g)
	h.close()

	#manually lok up desired planet and embryo/planetesimal divide
	target_i_file = earth_list[counter]
	emb_divide = emb_list[counter]
	GI_mass = 0.
	cc_mass = 0.
	#read in final object mass and distance file
	with open( filepath + "/finalplanets.txt", 'r') as f1:
            lines =  f1.read().splitlines()
            
            for line in lines:
                index = float(line.split(None,6)[1])
                mass = float(line.split(None,6)[2])
                     

                if int(index) == target_i_file:

                        final_mass  =  mass

	#read in initial object orbital distance file
	with open( filepath + "/aorig.dat", 'r') as aorig:
            lines =  aorig.read().splitlines()
            for line in lines:
                nindex.append(line.split(None,2)[0])
                au.append(line.split(None,2)[1])
                multiAU.append([line.split(None,2)[0],line.split(None,2)[1]])
		
		
	
	with open(filepath + "/planetgrowth.out",'r') as f:

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

			if initial:  #initialize values of each gas speices
				if target_i == target_i_file:
					print ("get target index")
					target_r = get_radius.r_cm(target_m)
					g_accel = cgrav*target_m/(target_r)**2 
				
					water_a=1e18 # target_m*float(float(SF_variable_rp.water_mf(target_au))
					carbon_a=1e18 #target_m*float(SF_variable_rp.carbon_mf(target_au)*0.005
					nitrogen_a= 1e18#target_m*float(SF_variable_rp.nitrogen_mf(target_au)*0.005
					carbon_c = 1e16				
					atm_grams= water_a +  carbon_a + nitrogen_a
					

					Psurf  = atm_grams*g_accel_earth/SAearth
					Psurf_bar =Psurf*1e-6				#mantle_grams = target_m*float(SF_variable_rp.water_mf(target_au) +target_m*float(SF_variable_rp.carbon_mf(target_au)
					water_m= target_m*SF_variable_rp.water_mf(target_au)
					carbon_m=  target_m*SF_variable_rp.carbon_mf(target_au)
					nitrogen_m= target_m*SF_variable_rp.nitrogen_mf(target_au)
					accret_tot = water_m+ carbon_m  +nitrogen_m

					mantle_grams =   water_m + carbon_m+ nitrogen_m 

					impactEsc = 0
					
                    #save initialized values
					Psurf_List.append(round(Psurf_bar,3))
					waterfrac_List.append('{:0.3e}'.format(water_a))
					fracEarth_List.append(round(target_m/Mearth,3))
					nitrogenfrac_List.append('{:0.3e}'.format(nitrogen_a))
					carbonfrac_List.append('{:0.3e}'.format(carbon_a))
					atm_List.append('{:0.3e}'.format(atm_grams))
					mantle_List.append('{:0.3e}'.format(mantle_grams))
					time_List.append(round(float(time),3))
					esctot_List.append(round(float(impactEsc),3))
					waterm_List.append('{:0.3e}'.format(water_m))
					carbonm_List.append('{:0.3e}'.format(carbon_m))
					ingassed_List.append('{:0.3e}'.format(ingass_grams))
					nitrogenm_List.append('{:0.3e}'.format(nitrogen_m))
					ts_List.append('{:0.3e}'.format(Ts))
					accret_List.append('{:0.3e}'.format(accret_tot))
					#print('escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2))
                    
			#starts main calculation here
			if target_i ==  target_i_file:  #change target index to desired embryo
		
				Ta =1200
				target_r = get_radius.r_cm(target_m)  #function to calc radius from SA Jacobson
				g_accel = cgrav*target_m/(target_r)**2   #cm/s^2
				SA_planet = 4*np.pi*target_r**2


				initial = False
				if int(project_i) <= emb_divide: # for embryos
					print ("embryo")
					print(carbon_a,nitrogen_a,carbon_m,nitrogen_m,carbon_c)
					##########################deposit/accret gas$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
					emb_water, emb_carbon, emb_nitrogen = run_embryo_f.calc(project_i,emb_divide)
					#water_a= water_a + emb_water
					#carbon_a= carbon_a + emb_carbon
					#nitrogen_a= nitrogen_a + emb_nitrogen
					print(emb_water, emb_carbon, emb_nitrogen)
					#accret_tot = emb_water + emb_carbon + emb_nitrogen

					#atm_grams= water_a +nitrogen_a + carbon_a 
                    
					water_a = water_a + emb_water
					carbon_a = carbon_a + emb_carbon
					nitrogen_a = nitrogen_a + emb_nitrogen
					accret_tot = emb_water + emb_carbon + emb_nitrogen

					atm_grams =   water_a + carbon_a+ nitrogen_a 


					water_mm = 18.  # in g /mol
					carbon_mm = 12.
					nitrogen_mm = 28.

					water_mol_a = water_a/18.
					carbon_mol_a = carbon_a/12.   # g/(g/mol)
					nitrogen_mol_a = nitrogen_a/14

					total_mol =  water_mol_a + carbon_mol_a + nitrogen_mol_a

					water_molfrac = water_mol_a/total_mol
					carbon_molfrac =  carbon_mol_a/total_mol
					nitrogen_molfrac = nitrogen_mol_a/total_mol
                    
                    #calculate radius and other vars needed for mass loss scheme
					target_r = get_radius.r_cm(target_m)
					Vesc = math.sqrt((2*cgrav*target_m)/target_r)
					kai = (impact_v*project_m)/(Vesc*target_m)

					#########################impact loss from Schitling et al. (2015)###################
                    
					impact_frac = 0.4*kai + 1.8*kai**2  - 1.2*(kai)**3

					impactEsc = impact_frac*atm_grams
		
					impactEsc_save =  impactEsc_save + impactEsc  #save atm loss res

					water_a  = water_a -impactEsc*water_molfrac 
					carbon_a = carbon_a - impactEsc*carbon_molfrac 
					nitrogen_a = nitrogen_a - impactEsc*nitrogen_molfrac  

					atm_grams= water_a +nitrogen_a +carbon_a 
                #########################hydrodynamics escape, energy-limited##############
					Feuv=29.7*(time/1e9)**(-1.23)*(1.23)**(-2)
					m_elim = (np.pi*Feuv*target_r**3)/(cgrav*target_m)

					water_a  = water_a - m_elim*water_molfrac
					carbon_a = carbon_a - m_elim*carbon_molfrac
					nitrogen_a = nitrogen_a - m_elim*nitrogen_molfrac

					atm_grams= water_a +nitrogen_a +carbon_a
					atm_grams  = max(1,atm_grams)

					#############################Henry's Law###################
					Psurf  = atm_grams*g_accel/SA_planet
					Psurf_bar = Psurf*1e-6
					#Tsurf = (Psurf*Vearth)/(water_mol*R_gas)

					water_m_eq= (henry.kH_t(Ta,'water')*Psurf_bar*water_molfrac)*water_mm*target_m*0.6
					carbon_m_eq  = (henry.kH_t(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm*target_m*0.6
					nitrogen_m_eq = (henry.kH_t(Ta,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm*target_m*0.6

					if water_m_eq > water_m and (henry.kH_t(Ta,'water')*Psurf_bar*water_molfrac)*water_mm > water_m_eq/(target_m*0.6):						
							water_a = water_a -(water_m_eq -water_m)
							water_m = water_m_eq
					else:
							water_a = water_a +(water_m -  water_m_eq)
							water_m = water_m_eq
                            
					if carbon_m_eq > carbon_m and (henry.kH_t(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm > carbon_m_eq/(target_m*0.6):						
							carbon_a = carbon_a -(carbon_m_eq -carbon_m)
							carbon_m = carbon_m_eq
					else:
							carbon_a = carbon_a +(carbon_m - carbon_m_eq)
							carbon_m = carbon_m_eq
                            
					if nitrogen_m_eq > nitrogen_m and (henry.kH_t(Ta,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm > nitrogen_m_eq/(target_m*0.6):						
							nitrogen_a = nitrogen_a -(nitrogen_m_eq - nitrogen_m)
							nitrogen_m = nitrogen_m_eq
					else:
							nitrogen_a = nitrogen_a +(nitrogen_m - nitrogen_m_eq)
							nitrogen_m = nitrogen_m_eq
												
					water_a = max(10.,water_a)
					carbon_a = max(10.,carbon_a)
					nitrogen_a = max(10.,nitrogen_a)
                        
					atm_grams= water_a +nitrogen_a +carbon_a 

					water_m = max(10.,water_m)
					carbon_m = max(10.,carbon_m)
					nitrogen_m = max(10.,nitrogen_m)

					mantle_grams =   water_m + carbon_m + nitrogen_m
                    
					#######core partition#####
					Rhopl = 3.
					m_react = (0.15*(4/3*np.pi*target_r**3)*((8*np.pi*cgrav*Rhopl*target_r**2)/(3*e_m)))*project_rho  #from Deguen+2011
					carbon_c_eq = (D_c*(carbon_m/m_react))*((1/3*target_m)*(m_react/(target_m*2/3)))
					if carbon_c_eq > carbon_c and (D_c*(carbon_m/m_react)) > carbon_c_eq/((1/3*target_m)*(m_react/(target_m*2/3))):						
						carbon_m = carbon_m -(carbon_c_eq -carbon_c)
						carbon_c = carbon_c_eq

                        
					atm_grams= water_a +nitrogen_a + carbon_a 
						#print 'target_m=',atm_grams,water_a,nitrogen_a,carbon_a , argon_a,Mcap
					Psurf  = atm_grams*g_accel/SA_planet
					Psurf_bar =Psurf*1e-6
				
					last_GI_time = time
					GI_mass = GI_mass + project_m
					print(carbon_a,nitrogen_a,carbon_m,nitrogen_m)

					Psurf_List.append(round(Psurf_bar,3))  
					waterfrac_List.append('{:0.3e}'.format(water_a))
					fracEarth_List.append(round(target_m/Mearth,3)) #target_m/final_mass
					nitrogenfrac_List.append('{:0.3e}'.format(nitrogen_a))
					carbonfrac_List.append('{:0.3e}'.format(carbon_a))
					atm_List.append('{:0.3e}'.format(atm_grams))
					mantle_List.append('{:0.3e}'.format(mantle_grams))
					time_List.append(round(float(time),3))
					esctot_List.append(round(float(impactEsc),3))
					waterm_List.append('{:0.3e}'.format(water_m))
					carbonm_List.append('{:0.3e}'.format(carbon_m))
					carbonc_List.append('{:0.3e}'.format(carbon_c))
					ingassed_List.append('{:0.3e}'.format(ingass_grams))
					nitrogenm_List.append('{:0.3e}'.format(nitrogen_m))
					ts_List.append('{:0.3e}'.format(Ts))
					accret_List.append('{:0.3e}'.format(accret_tot))
                    		

				else:   #for super-planetesiamls assembled from plantesiamls.

					print ("planetesimal")
				
					target_r = get_radius.r_cm(target_m)  #use function from SA Jacobson
					g_accel = cgrav*target_m/(target_r)**2   #cm/s^2
					SA_planet = 4*np.pi*target_r**2

					#gather mass of CCs
					if project_au >= 2.5:
						cc_mass = cc_mass +project_m
					while Mcurr <= project_m:
								
						Rpl =  random.choice(ab_radii)*1e5  #from km to cm
						Rhopl = 3.
						Mpl = Rhopl*(4./3*np.pi*Rpl**3)		
									
						Mcurr  = Mcurr + Mpl


						##########################deposit/accret gas$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
						#water_a= water_a + Mpl*float(SF_variable_rp.water_mf(project_au)
						#carbon_a= carbon_a + Mpl*float(SF_variable_rp.carbon_mf(project_au)
						#nitrogen_a= nitrogen_a + Mpl*float(SF_variable_rp.nitrogen_mf(project_au)
						#atm_grams= water_a +nitrogen_a +carbon_a + argon_a
				
						water_a= water_a + Mpl*float(SF_variable_rp.water_mf(project_au))
						carbon_a= carbon_a + Mpl*float(SF_variable_rp.carbon_mf(project_au))
						nitrogen_a= nitrogen_a + Mpl*float(SF_variable_rp.nitrogen_mf(project_au))
						
						accret_tot = Mpl*float(SF_variable_rp.water_mf(project_au)) + Mpl*float(SF_variable_rp.carbon_mf(project_au)) +  Mpl*float(SF_variable_rp.nitrogen_mf(project_au))
						#print 'target_m=',atm_grams,water_grams ,nitrogen_grams,carbon_grams
						atm_grams =   water_a  + carbon_a  +  nitrogen_a

						water_mm = 18.  # in g /mol
						carbon_mm = 12.
						nitrogen_mm = 28.				
					
						#print "water_mol=", henry_water_grams
						water_mol_a = water_a/18.
						
						carbon_mol_a = carbon_a/12.   # g/(g/mol)
			
						nitrogen_mol_a = nitrogen_a/28.
						
						
						total_mol_a =  water_mol_a + carbon_mol_a + nitrogen_mol_a
						if total_mol_a > 0:

							water_molfrac = water_mol_a/total_mol_a
							carbon_molfrac = carbon_mol_a/total_mol_a
							nitrogen_molfrac = nitrogen_mol_a/total_mol_a
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
						
						atm_grams = max(100, atm_grams)
						mu =  (carbon_a/atm_grams)*carbon_mm+  (water_a/atm_grams)*water_mm+ (nitrogen_a/atm_grams)*nitrogen_mm      #calculate mean molecular weight 
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
						
						water_a  = water_a -impactEsc*water_molfrac 
						carbon_a = carbon_a - impactEsc*carbon_molfrac 
						nitrogen_a = nitrogen_a - impactEsc*nitrogen_molfrac 
					

						atm_grams= water_a +nitrogen_a +carbon_a 
						  #avoid zero atm values

                       #########################hydrodynamics escape, energy-limited##############
						Feuv=29.7*(time/1e9)**(-1.23)*(1.23)**(-2)
						m_elim = (np.pi*Feuv*target_r**3)/(cgrav*target_m)

						water_a  = water_a - m_elim*water_molfrac
						carbon_a = carbon_a - m_elim*carbon_molfrac
						nitrogen_a = nitrogen_a - m_elim*nitrogen_molfrac
        	                                
						atm_grams= water_a +nitrogen_a +carbon_a 

						atm_grams  = max(1,atm_grams)
						##############Henry's Law####################
						Psurf  = atm_grams*g_accel/SA_planet
						Psurf_bar = Psurf*1e-6
						#Tsurf = (Psurf*Vearth)/(water_mol*R_gas)
						
						water_m_eq= (henry.kH_tf(Ta,'water')*Psurf_bar*water_molfrac)*water_mm*target_m*0.6  #molarity * molar mass = grams/litter, (molarity * mass)* solution mass = grams
						carbon_m_eq  = (henry.kH_tf(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm*target_m*0.6
						nitrogen_m_eq = (henry.kH_tf(Ta,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm*target_m*0.6


						if water_m_eq > water_m: #and (henry.kH_t(Ta,'water')*Psurf_bar*water_molfrac)*water_mm > water_m_eq/(target_m*0.6):						
							water_a = water_a -(water_m_eq -water_m)
							water_m = water_m_eq
						else:
							water_a = water_a +(water_m - water_m_eq)
							water_m = water_m_eq
                            
						if carbon_m_eq > carbon_m: #and (henry.kH_t(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm > carbon_m_eq/(target_m*0.6):						
							carbon_a = carbon_a -(carbon_m_eq -carbon_m)
							carbon_m = carbon_m_eq
						else:
							carbon_a = carbon_a +(carbon_m - carbon_m_eq)
							carbon_m = carbon_m_eq
                            
						if nitrogen_m_eq > nitrogen_m: #and (henry.kH_t(Ta,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm > nitrogen_m_eq/(target_m*0.6):						
							nitrogen_a = nitrogen_a -(nitrogen_m_eq - nitrogen_m)
							nitrogen_m = nitrogen_m_eq
						else:
							nitrogen_a = nitrogen_a +(nitrogen_m - nitrogen_m_eq)
							nitrogen_m = nitrogen_m_eq
                            
												
						water_a = max(10.,water_a)
						carbon_a = max(10.,carbon_a)
						nitrogen_a = max(10.,nitrogen_a)
                        
						atm_grams= water_a +nitrogen_a +carbon_a 
                        
						water_m = max(10.,water_m)
						carbon_m = max(10.,carbon_m)
						nitrogen_m = max(10.,nitrogen_m)

						mantle_grams = water_m + carbon_m + nitrogen_m

						#######core partition#####
						m_react = (0.15*(4/3*np.pi*target_r**3)*((8*np.pi*cgrav*Rhopl*target_r**2)/(3*e_m)))*project_rho  #from Deguen+2011
						carbon_c_eq = (D_c*(carbon_m/m_react))*(1/3*project_m)*(m_react/(target_m*2/3))
						if carbon_c_eq > carbon_c: #and (D_c*(carbon_m/m_react)) > carbon_c_eq/((1/3*target_m)*(m_react/(target_m*2/3))):						
							#carbon_m = carbon_m -(carbon_c_eq -carbon_c)
							carbon_c = carbon_c_eq
						else:
							#carbon_m = carbon_m +(carbon_c - carbon_c_eq)
							carbon_c = carbon_c_eq
                    
						water_m = max(10.,water_m)
						carbon_m = max(10.,carbon_m)
						nitrogen_m = max(10.,nitrogen_m)
						
    
						mantle_grams =   water_m + carbon_m + nitrogen_m
						#print 'target_m=',atm_grams,water_a,nitrogen_a,carbon_a
                        
						Psurf  = atm_grams*g_accel/SA_planet
						Psurf_bar =Psurf*1e-6
						
						
						#print F_atm, sigma,Teq,tau,Ts,atm_grams
						#if math.isnan(F_atm):
						#	sys.exit()
						if (Rpl > 1e7 or random.uniform(0,1) < 1e-4):
							#print(carbon_a,nitrogen_a,carbon_m)

							fracEarth2_List.append(round(target_m/Mearth,3))  #Mearth
							rpl2_List.append(float(Rpl))
							esc2_List.append('{:0.3e}'.format(impactEsc))
							ingassed2_List.append('{:0.3e}'.format(ingass_grams))
							atm2_List.append('{:0.3e}'.format(atm_grams))
							mantle2_List.append('{:0.3e}'.format(mantle_grams))
							accret2_List.append('{:0.3e}'.format(accret_tot))	
					###################Append to various output lists###############################
					

					
				#	kappa = ((kappa_0*g_accel)/(3*Psurf))**(1/2.)
				#	tau = (3*kappa*atm_grams)/(8*np.pi*target_r**2)
				#	F_atm = 2*sigma*(Ts**4-Teq**4)/(tau +2)
			#		Ts = (F_atm/sigma + Teq**4)**(1/4.)
			#		if len(time_List) == 1.:
				#		deltaM = target_m
				#		deltaT = (time)*3.154e7
				#	else:		
				#		deltaM = project_m
					#	deltaT = (time - time_List[-1])*3.154e7
					#Ts= (0.5*deltaM*impact_v**2/2. - 4*np.pi*target_r**2*F_atm*deltaT)/(Cp*deltaM)   +  300
				#	Ts = ((tau + 2)*(h_frac*deltaM*impact_v**2)/(16*sigma*np.pi*target_r**2*deltaT) + Teq**4)**(1/4#.)
				#	if Ts <=0.:
					#	Ts = 500
				
					Psurf_List.append(round(Psurf_bar,3))
					waterfrac_List.append('{:0.3e}'.format(water_a))
					fracEarth_List.append(round(target_m/Mearth,3))  #target_m/Mearth
					nitrogenfrac_List.append('{:0.3e}'.format(nitrogen_a))
					carbonfrac_List.append('{:0.3e}'.format(carbon_a))
					atm_List.append('{:0.3e}'.format(atm_grams))
					mantle_List.append('{:0.3e}'.format(mantle_grams))
					time_List.append(round(float(time),3))
					esctot_List.append(round(float(impactEsc),3))
					waterm_List.append('{:0.3e}'.format(water_m))
					carbonm_List.append('{:0.3e}'.format(carbon_m))
					carbonc_List.append('{:0.3e}'.format(carbon_c))
					ingassed_List.append('{:0.3e}'.format(ingass_grams))
					nitrogenm_List.append('{:0.3e}'.format(nitrogen_m))
					ts_List.append('{:0.3e}'.format(Ts))
					accret_List.append('{:0.3e}'.format(accret_tot))
						#amantle_List.append(mantle_grams)
					#nitrogen_List.append(nitrogen_grams)
					#water_List.append(water_grams)
					#carbon_List.append(carbon_grams)
							
					#atm_List.append(atm_grams)

					
					#fracEarth_List.append(target_m/Mearth)

					#print('escape=',impactEsc,'Psurf = ', Psurf,'atm_grams=', atm_grams)

				       # print 'water_frac=', water_grams/atm_grams,'nitro_frac=',nitrogen_grams/atm_grams,'Psurf=', Psurf_bar,'mantle=', mantle_grams
					#print('escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2))

				Mcurr = 0
				mearth = target_m*Mearth
				earthm_List.append(mearth)
			#	time_List.append(float(time))

			#print('water_mol=',water_mol,'Psurf = ',Psurf,'atm_grams=', atm_grams,'mantle_grams=',mantle_grams)
			#withAU1 = line[:31] + str(target_au)  + line[31:]
			#withAU2 = withAU1[:-24] + str(project_au)  + withAU1[-24:]
			#newdata.append(withAU2)

	######################save data################


	with open('plcomp/OUT_ATM_cc_' + str(CC_water[z]), 'w') as f2:
		writer = csv.writer(f2, delimiter='\t')
		#writer.write(zip("Earth", "Psurf", "atm", "water", "carbon", "nitrogen",  "mantle"))
		f2.write("%s %s %s %s %s %s %s" % ("Earth   " +str(final_mass) +" "+ str(emb_divide),  "Psurf   ", "atm   ", "water   ", "carbon   ", "nitrogen", "time  final GIt = " + str(last_GI_time) + "  final GIm =" + str(GI_mass) +" tot CC =" + str(cc_mass) + "\n") ) 
		writer.writerows(zip(fracEarth_List,Psurf_List,atm_List,waterfrac_List,carbonfrac_List,nitrogenfrac_List,time_List))

	with open('plcomp/OUT_MANT_cc_' + str(CC_water[z]), 'w') as f4:
		writer = csv.writer(f4, delimiter='\t')
		f4.write("%s %s %s %s %s %s %s" % ("Earth", "CoreC", "BulkMantle", "Water  ", "carbon ", "Nitrogen", "time \n"))
		writer.writerows(zip(fracEarth_List, carbonc_List,  mantle_List, waterm_List,carbonm_List, nitrogenm_List,time_List))


	#with open('nbody_sims/n_out/OUT_FLUX_' + outname, 'w') as f6:
	#	writer = csv.writer(f6, delimiter='\t')
		#f6.write("%s %s %s %s %s %s %s" % ("Earth", "CoreC", "BulkMantle", "Water  ", "carbon ", "Nitrogen", "time \n"))
		#
      #  writer.writerows(zip(fracEarth_List, carbonc_List,  mantle_List, waterm_List,carbonm_List, nitrogenm_List,time_List))

