import numpy as np

def EvaluateDOS(ssfreqs, natoms):
	ThzToCm = 33.35640952
	Na = 6.0221367e23
	kb = 1.3806488e-23
	R = kb*Na
	dim = 3 * natoms

	w_0 = ssfreqs.flatten() * ThzToCm
	ssfreqs = np.reshape(ssfreqs,(dim,-1)) * ThzToCm

	use_imaginary = 0
	if use_imaginary == 0:
		for i in range(np.size(w_0)):
			if w_0[i] < 0:
				w_0[i] = 0
	
	#Create xrange
	xS = np.amin(w_0) - 100#Get Dist +- 100 from freq range in ssfreqs
	xL = np.amax(w_0) + 100
	grid = int(dim * (np.size(ssfreqs) /(dim)) + 200) #100 for xS  + xL
	x = np.linspace(xS,xL,grid)
	
	#find avg dx
	dx = []
	for i in np.arange(np.size(x) - 1):
		tmp = x[i + 1] - x[i]		
		dx = np.append(dx,tmp)
	dx = np.mean(dx)

	#Evaluate PDF in xrange. bw_method is bw; change if necessary
	y = np.zeros(grid)
	bw = 5 #in wavenumbers ; this is consistent with the std / 20 that NyDay 2016 used (Just for all modes)
	for i in np.arange(dim):
		band = ssfreqs[i,:]
		for j in np.arange(np.size(band)):
			mean = np.mean(band)
			g_x = np.exp(-np.square(x - mean) / (2  * np.square(bw)))
			y = y + g_x

	#Normalized Density of States 
	y_0 = y / (np.size(x) * bw * np.sqrt(2 * np.pi))

	#Create plot
	#axes = plt.gca()
	#axes.set_xlim([xS,xL])
	#axes.set_ylim([0,np.amax(y_0)+0.001])

	#Check norm; print A_0 if needed
	A_0 = 0
	for i in range(np.size(x) - 1):
		A_0 = A_0 + y_0[i] * dx

	tL = 500
	tS = 0
	T = np.arange(tS,tL,10)
	FvibArray = []
	HvibArray = []
	SvibArray = []

	for temp in T:
		thermo = Thermo(x,temp)
	
	##FVIB##	
		Fw_0 = thermo.Helmholtz() * y_0
		Fw_0[0] = 0
		F_vib = 0

		for i in range(np.size(x)):
			flag = np.isfinite(Fw_0[i])
			if flag == True:
				F_vib = F_vib + Fw_0[i] * dx

		F_vib = 3 * natoms * Na * F_vib / 1000
		FvibArray = np.append(FvibArray, F_vib)
	
	##ENTHALPY##
		Hw_0 = thermo.Enthalpy() * y_0
		Hw_0[0] = 0
		H_vib = 0
	
		for i in range(np.size(x)):
			flag = np.isfinite(Hw_0[i])
			if flag == True:
				H_vib = H_vib + Hw_0[i] * dx
	
		H_vib = 3 * Na * natoms * H_vib / 1000
		HvibArray = np.append(HvibArray, H_vib)
	
	##ENTROPY##
		Sw_0 = thermo.Entropy() * y_0
		Sw_0[0] = 0
		Svib = 0
	
		for i in range(np.size(x)):
			flag = np.isfinite(Sw_0[i])
			if flag == True:
				Svib = Svib + Sw_0[i] * dx

		Svib = 3 * R * Svib * natoms # // Fultz 2009
		SvibArray = np.append(SvibArray, Svib)

	#AllData = pd.DataFrame({'Temperature': np.array(T),
                                    #'Hvib': np.array(HvibArray),
                                    #'Svib': np.array(SvibArray),
                                    #'Fvib': np.array(FvibArray)})


	return FvibArray, HvibArray, SvibArray, T

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
		return F_fun

	def Enthalpy(self):
		if self.T == 0:
			H_fun = self.w * self.h / 2
		else:
			H_fun = self.h / 2 * self.w * (np.cosh(self.h * self.w / (2 * self.kb * self.T)) / np.sinh(self.h * self.w / (2 * self.kb * self.T)))
		return H_fun

	def Entropy(self):
		eps = self.w * self.h #DOS in energy
		occ = 1 / (np.exp(eps / (self.kb * self.T)) - 1 ) #Planck Distribution // Fultz 2009
		S_fun = (occ + 1) * np.log(occ + 1) - occ * np.log(occ) #Svib for phonons
		return S_fun

	def HeatCap(self):
		Cv_fun = self.kb * np.square((self.h * self.w / (2 * self.kb * self.T))) * np.square(1 / (np.sinh(self.h * self.w / (2 * self.kb * self.T))))
		return Cv_fun
