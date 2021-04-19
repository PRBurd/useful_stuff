
# --------------------------------------------------------
# based on https://arxiv.org/pdf/astro-ph/9905116.pdf
# --------------------------------------------------------
import numpy as np
from scipy.integrate import quad


if __name__ == '__main__':
	print('bitte als library irgendwo einbinden')
else:
	print('successfully loaded')


class cosmology:
	#comology parameters from MOJAVE web page omega_k = 0 -> flat spacetime, redshift z (cosmology)
	def __init__(self,H0=71.0,omega_M = 0.27, omega_Lambda = 0.73, omega_k = 0, z = 1):
		self.H0 = H0
		self.omega_M = omega_M
		self.omega_Lambda = omega_Lambda
		self.omega_k = omega_k
		self.z = z

	#transverse distance between 2 signals
	#@property
	def D_M(self):
		if self.omega_k+self.omega_M+self.omega_Lambda != 1:
			print(r'make sure that $\Omega _k+\Omega _M + \Omega _{\Lambda}$ = 1')
			return 0
		H_0 = (self.H0/3.086)*10**(-19)
		D_h = 299792458/(H_0)
		E = lambda redshift,o_k,o_M,o_Lambda:1/np.sqrt(o_M*(1+redshift)**3+o_k*(1+redshift)**2+o_Lambda)
		D_c = D_h * quad(E,0,self.z,args=(self.omega_k,self.omega_M,self.omega_Lambda))[0]
		if self.omega_k > 0:
			Dc = D_h * (1/np.sqrt(self.omega_k)) * np.sinh(np.sqrt(self.omega_k)*D_c/D_h)
		if self.omega_k == 0:
			Dc = D_c
		if self.omega_k < 0:
			Dc = D_h * (1/np.sqrt(self.omega_k)) * np.sin(np.sqrt(self.omega_k)*D_c/D_h)
		return(Dc) # in m

	#angular diameter distance
	#@property
	def D_A(self):
		return(cosmology.D_M(self) /(1+self.z)) #returns m/rad

	#distance line-of-sight
	#@property
	def D_c(self):
		if self.omega_k+self.omega_M+self.omega_Lambda == 1:
			H_0 = (self.H0/3.086)*10**(-19)
			D_h = 299792458/(H_0)
			E = lambda redshift,o_k,o_M,o_Lambda:1/np.sqrt(o_M*(1+redshift)**3+o_k*(1+redshift)**2+o_Lambda)
			Dc = D_h * quad(E,0,self.z,args=(self.omega_k,self.omega_M,self.omega_Lambda))[0]
			return(Dc) # in m
		else:
			print(r'make sure that $\Omega _k+\Omega _M + \Omega _{\Lambda}$ = 1')
			return(0)


	#lookback time
	#@property
	def t_l(self):
		if self.omega_k+self.omega_M+self.omega_Lambda != 1:
			print(r'make sure that $\Omega _k+\Omega _M + \Omega _{\Lambda}$ = 1')
			return(0)
		H_0 = (self.H0/3.086)*10**(-19)
		t_h = 1/(H_0)
		E = lambda redshift,o_k,o_M,o_Lambda:1/((1+redshift)*np.sqrt(o_M*(1+redshift)**3+o_k*(1+redshift)**2+o_Lambda))
		tl = t_h * quad(E,0,self.z,args=(self.omega_k,self.omega_M,self.omega_Lambda))[0]
		return(tl) # in s

	#--------------------
	#helpful conversions
	#--------------------
	#convert m to Gpc
	def m_to_kpc(self,d):
		return(d /(3.086*10**(16) *10**(3)))
	def m_to_gpc(self,d):
		return(d /(3.086*10**(16) *10**(9)))

	#convert m to Gly
	def m_to_gly(self,d):
		return(d / 9.454254955488E+24)

	#convert Gpc to Gly
	def gpc_to_gly(self,d):
		return(d * 3.2638)

	#convert Gly to Gpc
	def gly_to_gpc(self,d):
		return(d / 3.2638)

	#convert s to Yr
	def s_to_yr(self,t):
		return(t/(60*60*24*365))

	#convert s to GYr
	def s_to_gyr(self,t):
		return(t/(60*60*24*365*10**9))

	def s_to_h(self,t):
		return(t*((self.H0/3.086)*10**(-19)))

	@property
	def kpc_per_arcsec(self):
		return(cosmology.m_to_kpc(self,np.pi/(180*3600) * cosmology.D_A(self)))
	@property
	def kpc_per_arcmin(self):
		return(cosmology.m_to_kpc(self,np.pi/(180*60) * cosmology.D_A(self)))
	@property
	def Gyrs_after_BigBang(self):
		return(cosmology.s_to_gyr(self,1/(((self.H0/3.086)*10**(-19))) - cosmology.t_l(self)))
	
	@property
	def lookback_time_s(self):
		return(cosmology.t_l(self))
	@property
	def lookback_time_Gyr(self):
		return(cosmology.s_to_gyr(self,cosmology.t_l(self)))
	@property 
	def lookback_time_Ht(self):
		return(cosmology.s_to_h(self,cosmology.t_l(self)))
	@property
	def comoving_distance_m(self):
		return(cosmology.D_M(self))
	@property
	def comoving_distance_Gpc(self):
		return(cosmology.m_to_gpc(self,cosmology.D_M(self)))
	@property
	def comoving_distance_Gly(self):
		return(cosmology.s_to_gyr(self,cosmology.D_M(self)))


