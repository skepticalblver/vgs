#!/usr/bin/env python
import numpy as np
import math
import SF    #step function parameterization for volatile frac
import henry  #solubility constants
import random
import get_radius  #func to calc embryo radius 
import glob


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

Vesc= 11e5   #cm/s earth escape speed

Ts =1000 # initital surface temp in K
h_frac = 1.
Cp = 1e7 #specific heat of surface layer in ergs
kappa_0  = 0.1  # absorption coefficient a the surface, in cm2/g, from Matsui & Abe 1986
Teq = 265  #Teq of Earth in K
sigma = 5.6704e-5  #Stefan-Boltzman's constant,  erg/cm2/s/K4

earth_list=[9,9,6,20,16,18,6,12,19,28,9,15,8,7,23,5]

emb_list=[45,45,45,69,69,90,90,90,58,110,110,110,110,71,71,71]

CompArray =  np.empty(5)

def calc(input_index,input_div):
    print("in function")
    Mcurr=0
    counter = 0
    henry_water_save= 0.
    henry_nitrogen_save= 0.
    henry_carbon_save= 0.
    mantle_grams = 0.
    ingass_grams = 0.
    impactEsc= 0.
    impactEsc_save = 0.
    water_tot=0.
    carbon_tot=0.
    nitrogen_tot=0.
    target_rho= 3.
    Rhopl=3.
    D_c = 3000  #partition coefficient for carbon
    for filename in sorted(glob.glob('trap_sims/sim*')):
    	initial = True #if true, run initialization8to1-0.5-10
    	#initialize list to store various atmophiles      

    
    	au=[] 
    	nindex=[] 
    	multiAU=[] 
    
    	#manually lok up desired planet and embryo/planetesimal divide
    	target_i_file = input_index
    	emb_divide = input_div
    	GI_mass = 0.
    	cc_mass = 0.
    	#read in final object mass and distance file
  
    	#read in initial object orbital distance file
    	with open( filename + "/a_orig.dat", 'r') as aorig:
                lines =  aorig.read().splitlines()
                for line in lines:
                    nindex.append(line.split(None,2)[0])
                    au.append(line.split(None,2)[1])
                    multiAU.append([line.split(None,2)[0],line.split(None,2)[1]])
    		
    		
    	
    	with open(filename + "/planet_growth.out",'r') as f:
    
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
    					
    					target_r = get_radius.r_cm(target_m)
    					g_accel = cgrav*target_m/(target_r)**2 
    				
    					water_a=1e18 # target_m*SF.water_mf(target_au)*0.005
    					carbon_a=1e18 #target_m*SF.carbon_mf(target_au)*0.005
    					nitrogen_a= 1e18#target_m*SF.nitrogen_mf(target_au)*0.005
    					carbon_c = 0.   									
    					atm_grams= water_a +  carbon_a + nitrogen_a
    
    					Psurf  = atm_grams*g_accel_earth/SAearth
    					Psurf_bar =Psurf*1e-6				#mantle_grams = target_m*SF.water_mf(target_au) +target_m*SF.carbon_mf(target_au)

    					water_m= target_m*SF.water_mf(target_au)
    					carbon_m=  target_m*SF.carbon_mf(target_au)
    					nitrogen_m= target_m*SF.nitrogen_mf(target_au)
    					accret_tot = water_m+ carbon_m  +nitrogen_m
    
    					mantle_grams =   water_m + carbon_m+ nitrogen_m 
    
    					impactEsc = 0
    					
                        #save initialized values
                        
    			#starts main calculation here
    			if target_i ==  target_i_file:  #change target index to desired embryo
    		
    				Ta =1200
    				target_r = get_radius.r_cm(target_m)  #function to calc radius from SA Jacobson
    				g_accel = cgrav*target_m/(target_r)**2   #cm/s^2
    				SA_planet = 4*np.pi*target_r**2
    
    
    				initial = False
    				if int(project_m)  > 0.1*Mearth: # for embryos
    					print ("sub_embryo")
    
    					##########################deposit/accret gas$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
    					
    					water_a= water_a + project_m*SF.water_mf(project_au)*0.01
    					carbon_a= carbon_a + project_m*SF.carbon_mf(project_au)*0.01
    					nitrogen_a= nitrogen_a+ project_m*SF.nitrogen_mf(project_au)*0.01
    					accret_tot = project_m*SF.water_mf(project_au) + project_m*SF.carbon_mf(project_au)+ project_m*SF.nitrogen_mf(project_au)
    
    					atm_grams= water_a +nitrogen_a + carbon_a 
    
    
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
    
    					#print ("impact_frac=", impact_frac)		
    					impactEsc_save =  impactEsc_save + impactEsc  #save atm loss res
    
    					water_a  = water_a -impactEsc*water_molfrac 
    					carbon_a = carbon_a - impactEsc*carbon_molfrac 
    					nitrogen_a = nitrogen_a - impactEsc*nitrogen_molfrac  
    
    					atm_grams= water_a +nitrogen_a +carbon_a 
                        #########################hydrodynamics escape, energy-limited##############
    					Feuv=29.7*(time+10/1e9)**(-1.23)*(1.23)**(-2)
    					m_elim = (np.pi*Feuv*target_r**3)/(cgrav*target_m)
    
    					water_a  = water_a - m_elim*water_molfrac
    					carbon_a = carbon_a - m_elim*carbon_molfrac
    					nitrogen_a = nitrogen_a - m_elim*nitrogen_molfrac
    
    					atm_grams= water_a +nitrogen_a +carbon_a
    					atm_grams = max(10, atm_grams)
					#############################Henry's Law###################
    					Psurf  = atm_grams*g_accel/SA_planet
    					Psurf_bar = Psurf*1e-6
    					#Tsurf = (Psurf*Vearth)/(water_mol*R_gas)

    					water_m_eq= (henry.kH_t(Ta,'water')*Psurf_bar*water_molfrac)*water_mm*target_m*0.6
    					carbon_m_eq  = (henry.kH_t(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm*target_m*0.6
    					nitrogen_m_eq = (henry.kH_t(Ta,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm*target_m*0.6

    					if water_m_eq > water_m:# and (henry.kH_t(Ta,'water')*Psurf_bar*water_molfrac)*water_mm > water_m_eq/(target_m*0.6):						
    						water_a = water_a -(water_m_eq -water_m)
    						water_m = water_m_eq
    					else:
    						water_a = water_a +(water_m - water_m_eq)
    						water_m = water_m_eq
                            
    					if carbon_m_eq > carbon_m:# and (henry.kH_t(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm > carbon_m_eq/(target_m*0.6):						
    						carbon_a = carbon_a -(carbon_m_eq -carbon_m)
    						carbon_m = carbon_m_eq
    					else:
    						carbon_a = carbon_a +(carbon_m -carbon_m_eq)
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

    					mantle_grams =   water_m + carbon_m + nitrogen_m
                    

    						#print 'target_m=',atm_grams,water_a,nitrogen_a,carbon_a , argon_a,Mcap
    					Psurf  = atm_grams*g_accel/SA_planet
    					Psurf_bar =Psurf*1e-6
    				
    					last_GI_time = time
    					GI_mass = GI_mass + project_m
    

    					#print('escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2))
    
    			
    				else:   #for super-planetesiamls assembled from plantesiamls.
    
    					#print ("planetesimal")
    
    				
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
    						#neon_a= neon_a + Mpl*SF.neon_mf(project_au)
    						#water_a= water_a + Mpl*SF.water_mf(project_au)
    						#carbon_a= carbon_a + Mpl*SF.carbon_mf(project_au)
    						#nitrogen_a= nitrogen_a + Mpl*SF.nitrogen_mf(project_au)
    						#atm_grams= water_a +nitrogen_a +carbon_a + argon_a
    				

    						water_a= water_a + Mpl*SF.water_mf(project_au)
    						carbon_a= carbon_a + Mpl*SF.carbon_mf(project_au)
    						nitrogen_a= nitrogen_a + Mpl*SF.nitrogen_mf(project_au)
    						
    						accret_tot = Mpl*SF.water_mf(project_au) + Mpl*SF.carbon_mf(project_au) +  Mpl*SF.nitrogen_mf(project_au)
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
                                
    						atm_grams = max(10,atm_grams)
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
    						
    						
    						mu =  (carbon_a/atm_grams)*carbon_mm+  (water_a/atm_grams)*water_mm+ (nitrogen_a/atm_grams)*nitrogen_mm      #calculate mean molecular weight 
    						mu = max(1,mu)
    						m_H = 1.67e-26 # hydrogen atomic mass in kg	
    						H = 1e2*(kerg*Ta)/(mu*m_H*g_accel*0.01)    #convert from m to cm 
    						
    						H_earth = H
    						Mcap = 2*np.pi*rho_surf*H_earth**2*target_r  
    						
    
    						r_min =  (3*rho_surf/Rhopl)**(1/3.)*H_earth  #1e5 for earth
    						r_max = (3*(2*np.pi)**(1/2)*rho_surf/4.*Rhopl)**(1/3.)*(H_earth*target_r)**(1/2.) #25e5 for earth
    						if isinstance(r_min, complex) or isinstance(r_max, complex):
    							break
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

    						
    						water_a  = water_a -impactEsc*water_molfrac 
    						carbon_a = carbon_a - impactEsc*carbon_molfrac 
    						nitrogen_a = nitrogen_a - impactEsc*nitrogen_molfrac 
    					
    
    						atm_grams= water_a +nitrogen_a +carbon_a 
    						atm_grams  = max(1,atm_grams)  #avoid zero atm values
    
                                            	#########################hydrodynamics escape, energy-limited##############
    						Feuv=29.7*(time+10/1e9)**(-1.23)*(1.23)**(-2)
    						m_elim = (np.pi*Feuv*target_r**3)/(cgrav*target_m)
    
    						water_a  = water_a - m_elim*water_molfrac
    						carbon_a = carbon_a - m_elim*carbon_molfrac
    						nitrogen_a = nitrogen_a - m_elim*nitrogen_molfrac
            	                                
    						atm_grams= water_a +nitrogen_a +carbon_a 
    
    						atm_grams = max(100, atm_grams)
    						##############Henry's Law####################
    						Psurf  = atm_grams*g_accel/SA_planet
    						Psurf_bar = Psurf*1e-6
    						#Tsurf = (Psurf*Vearth)/(water_mol*R_gas)
						
    						water_m_eq= (henry.kH_t(Ta,'water')*Psurf_bar*water_molfrac)*water_mm*target_m*0.6
    						carbon_m_eq  = (henry.kH_t(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm*target_m*0.6
    						nitrogen_m_eq = (henry.kH_t(Ta,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm*target_m*0.6

    						if water_m_eq > water_m: #and (henry.kH_t(Ta,'water')*Psurf_bar*water_molfrac)*water_mm > water_m_eq/(target_m*0.6):						
    							water_a = water_a -(water_m_eq -water_m)
    							water_m = water_m_eq
    						else:
    							water_a = water_a +(water_m -water_m_eq)
    							water_m = water_m_eq
                            
    						if carbon_m_eq > carbon_m: #and (henry.kH_t(Ta,'carbon')*Psurf_bar*carbon_molfrac)*carbon_mm > carbon_m_eq/(target_m*0.6):						
    							carbon_a = carbon_a -(carbon_m_eq -carbon_m)
    							carbon_m = carbon_m_eq
    						else:
    							carbon_a = carbon_a +(carbon_m -carbon_m_eq)
    							carbon_m = carbon_m_eq
                            
    						if nitrogen_m_eq > nitrogen_m: #and (henry.kH_t(Ta,'nitrogen')*Psurf_bar*nitrogen_molfrac)*nitrogen_mm > nitrogen_m_eq/(target_m*0.6):						
    							nitrogen_a = nitrogen_a -(nitrogen_m_eq - nitrogen_m)
    							nitrogen_m = nitrogen_m_eq
    						else:
    							nitrogen_a = nitrogen_a +(nitrogen_m -nitrogen_m_eq)
    							nitrogen_m = nitrogen_m_eq
												
    						water_a = max(10.,water_a)
    						carbon_a = max(10.,carbon_a)
    						nitrogen_a = max(10.,nitrogen_a)
                        
    						atm_grams= water_a +nitrogen_a +carbon_a 
                        
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
    				

    					#atm_List.append(atm_grams)
    
    					
    					#fracEarth_List.append(target_m/Mearth)
    
    					#print('escape=',impactEsc,'Psurf = ', Psurf,'atm_grams=', atm_grams)
    
    				       # print 'water_frac=', water_grams/atm_grams,'nitro_frac=',nitrogen_grams/atm_grams,'Psurf=', Psurf_bar,'mantle=', mantle_grams
    					#print('escape=','{:0.3e}'.format(impactEsc),'Psurf = ', round(Psurf_bar,3),'mantle_grams=', '{:0.3e}'.format(mantle_grams),'atm_grams=', '{:0.3e}'.format(atm_grams),'Project AU = ', round(project_au,2), 'Ts = ', Ts)
    
    				Mcurr = 0
    				mearth = target_m*Mearth

    				water_tot = water_a +water_m 
    				nitrogen_tot = nitrogen_a + nitrogen_m 
    				carbon_tot = carbon_a + carbon_m 

    	counter  = counter+1
    return (water_tot,carbon_tot,nitrogen_tot)
