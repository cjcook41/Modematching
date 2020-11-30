import numpy as np


class Thermo:
	def __init__(self, x, temp):
		self.T = temp
		self.h = 6.62607363e-34
		self.Na = 6.0221367e23
		self.kb = 1.3806488e-23
		self.c = 2.99792458e8
		self.w = x * (self.c * 100)
	def Helmholtz(self):
		if self.T == 0:
			F_fun = self.w * self.h / 2 
		else: 
			F_fun = self.kb * self.T * np.log(2*np.sinh(self.h * self.w / (2 * self.kb * self.T)))
			#F_fun = np.log(2*np.sinh(self.h * self.w / (2 * self.kb * self.T)))

		return F_fun
	def Enthalpy(self):
		if self.T == 0: 
			H_fun = self.w * self.h / 2
		else:
			H_fun = self.h / 2 * self.w * (np.cosh(self.h * self.w / (2 * self.kb * self.T)) / np.sinh(self.h * self.w / (2 * self.kb * self.T)))

		#H_fun = (self.h * self.w) / (np.exp( (self.h * self.w)/ (self.kb * self.T)) -1)

		return H_fun
	def Entropy(self):
		eps = self.w * self.h #DOS in energy
		#S_fun = self.kb * (eps / (2 * self.kb * self.T) * (np.cosh(eps / (2 * self.kb * self.T)) / np.sinh(eps / (2 * self.kb * self.T))) - np.log(2 * np.sinh(eps / (2 * self.kb * self.T))))  
		R = self.Na * self.kb
		occ = 1 / (np.exp(eps / (self.kb * self.T)) - 1 ) #Planck Distribution // Fultz 2009
		S_fun = (occ + 1) * np.log(occ + 1) - occ * np.log(occ) #Svib for phonons

		return S_fun
	def HeatCap(self):
		Cv_fun = self.kb * np.square((self.h * self.w / (2 * self.kb * self.T))) * np.square(1 / (np.sinh(self.h * self.w / (2 * self.kb * self.T))))
		return Cv_fun
