#!/usr/bin/env python
Mearth = 5.97e27 

def  r_cm(m_emb):
      #real*8 m,r
	m = (m_emb/Mearth)*1.1855e-4
	if (m<1.18554e-6): #then ! 0.01 Earth masses
		r=1.01854e-3*m**0.333333
	elif(m<2.37108e-4): #then ! 2.0 Earth masses
		r=6.34540e-4*m**0.298657
	elif(m<6.82761e-4): #then ! 5.8 Earth masses
		r=4.31273e-4*m**0.252394
	elif(m<1.08189e-2): #then ! 91 Earth masses
		r=1.55290e-2*m**0.744031
	elif(m<2.2943): #then ! 20000 Earth masses
		r=4.49484e-4*m**(-0.038559)
	else:
		r=2.18062e-4*m**0.832465
     
	return r*1.496e13  # from AU to cm
     
