import numpy as np
from . import GetDOS as dos
from scipy import optimize
import pandas as pd



def QHA(AllFreqs,FreqVols,ev_curve,count,natoms,Press,tL):
	print('STARTING QHA')
	Na = 6.0221367e23
	natoms = sum(natoms)
	dim = 3 * natoms
	AllFreqs = AllFreqs.transpose()

#Curve Data
	Volumes = ev_curve[:,0]
	Energies = ev_curve[:,1]
	Vmin = 0.98 * min(Volumes)
	Vmax = 1.08 * max(Volumes)

	"""
	Freqvol Range
	Larger Fvol_range -> less extrapolation, quadratic fit OK. Not true for narrow ranges
	"""
	Fvol_range = (max(FreqVols) - min(FreqVols)) / min(FreqVols)
	if Fvol_range > 0.1:
		polydeg = 2
	else:
		polydeg = 1

	TArray = []	
	OptVolArray = []
	OptGibbsArray = []
	OptFvibArray = [] 
	OptSvibArray = []
	OptHvibArray = []
	OptElArray = []

	FvibAll = []
	HvibAll = []
	SvibAll = []
	TAll = []

	##Get initial V for minimumization
	VGrid = np.arange(Vmin, Vmax, 0.01)

	for i in range(np.size(AllFreqs,1)):
		freqset = np.array(AllFreqs[:,i]).reshape(-1,dim)
		[Fvib, Hvib, Svib,T] = dos.EvaluateDOS(freqset,natoms,tL)
		print('Done with frequency set', i)
		FvibAll = np.append(FvibAll, Fvib)
		HvibAll = np.append(HvibAll, Hvib)
		SvibAll = np.append(SvibAll, Svib)
		TAll = np.append(TAll, T)
	

	ThermoDim = np.int(np.size(FvibAll)/count)
	FvibAll = np.array(FvibAll).reshape(-1,ThermoDim)
	HvibAll = np.array(HvibAll).reshape(-1,ThermoDim)
	SvibAll = np.array(SvibAll).reshape(-1,ThermoDim)
	TAll = np.array(TAll).reshape(-1,ThermoDim)
	np.savetxt("Fvibs.csv",FvibAll,delimiter=",")
	print('Starting pressure scan...')

	for P in Press:
	#Loop over all T's evaluated in Thermo (Columns in FvibAll)
		for i in range(ThermoDim):
			F_polyfit = np.polyfit(FreqVols, FvibAll[:,i], polydeg)
			H_polyfit = np.polyfit(FreqVols, HvibAll[:,i], polydeg)
			S_polyfit = np.polyfit(FreqVols, SvibAll[:,i], polydeg)
			PolyFvibFunction = np.polyval(F_polyfit,Volumes)
			PV = P * Volumes * 1.0e-24 * Na
			Gibbs = PolyFvibFunction + Energies + PV
			G_polyfit = np.polyfit(Volumes,Gibbs,4)

		#Define initial Gibbs curve \\ G(V) 
			def PolyGFxn(V):
				return np.polyval(G_polyfit,V)
		
			##Get initial V for minimumization
			T_fit = np.polyval(G_polyfit,VGrid) ##Pull these for G vs V curves

			initV = VGrid[np.argmin(T_fit)]
			Extrap = round(0.1 * len(VGrid))

			opt_params = optimize.minimize(PolyGFxn,initV,method='Powell') #Provides initial Vmin and Gmin for Murnaghan Fit // Maybe unnecessary?
			minimum = opt_params.x[0]
			minGibbs = opt_params.fun
			MinIndex = (np.abs(VGrid - minimum)).argmin()

		#Murnaghan EOS \\ Check parentheses
			def Murnaghan(x,c0,c1):
				return x[:,2] + c0*x[:,1] * (1 / (c1*(c1 - 1)) * np.power((x[:,0]/x[:,1]),1-c1) + 1 / c1 * (x[:,0]/x[:,1]) - 1 / (c1-1))

		#Figure out how to combine these fxns... it was an indexing issue
			def DoubleMurn(x,Vmin,Gmin,comp,exp): 
	
				def Murnaghan(x,c0,c1):
					return x[2] + c0*x[1] * (1 / (c1*(c1 - 1)) * np.power((x[0]/x[1]),1-c1) + 1 / c1 * (x[0]/x[1]) - 1 / (c1-1))
				
				Test = np.zeros((1,3))
				Test[:,0] = x
				Test[:,1] = Vmin
				Test[:,2] = Gmin

				if x > Vmin:
					ans = Murnaghan(Test[0,:],exp[0],exp[1])
				else:
					ans = Murnaghan(Test[0,:],comp[0],comp[1])
				return ans

		#Define Compression and Expansion branches for Double Murn
			if P < 0.75:
				VolExp = VGrid[MinIndex - Extrap:-1]
				VolComp = VGrid[0:MinIndex + Extrap]
				GibbsExp = T_fit[MinIndex - Extrap:-1]
				GibbsComp = T_fit[0:MinIndex + Extrap]

			else:
				VolExp = VGrid[MinIndex - int(0.5*Extrap):-1]
				VolComp = VGrid[0:MinIndex + Extrap]
				GibbsExp = T_fit[MinIndex - int(0.5*Extrap):-1]
				GibbsComp = T_fit[0:MinIndex + Extrap]

		#Get optimized parameters for Compression and Expansion branches
			a = [5.0,8.0]
			Comp = np.zeros((len(VolComp),3))
			Comp[:,0] = VolComp
			Comp[:,1] = minimum
			Comp[:,2] = minGibbs
			Exp = np.zeros((len(VolExp),3))
			Exp[:,0] = VolExp
			Exp[:,1] = minimum
			Exp[:,2] = minGibbs

		#Optimize parameters
			Exp_coeffs = optimize.curve_fit(Murnaghan,Exp, GibbsExp,a)[0]
			Comp_coeffs = optimize.curve_fit(Murnaghan,Comp,GibbsComp,a)[0]

			OptimalParams = optimize.minimize(DoubleMurn,minimum,args=(minimum,minGibbs,Comp_coeffs,Exp_coeffs),method='Powell')

			#All of this is for Gibbs curves -- uncomment for spreadsheets
			#Murn = []
			#for vol in VGrid:
			#	y = DoubleMurn(vol,minimum,minGibbs,Comp_coeffs,Exp_coeffs)
			#	t = np.array([vol,y])
			#	Murn.append(t)
			#df = pd.DataFrame(Murn,columns=['Volume','Gibbs'])
			#df.to_csv('Gibbs_curve-{}.csv'.format(TAll[1,i]))

			OptimalVolume = OptimalParams.x[0]
			OptimalGibbs = OptimalParams.fun
			OptimalFvib = np.polyval(F_polyfit,OptimalVolume)
			OptimalSvib = np.polyval(S_polyfit,OptimalVolume)
			OptimalHvib = np.polyval(H_polyfit,OptimalVolume)
			OptimalEl = OptimalGibbs - OptimalFvib # - OptimalPV

			TArray = np.append(TArray,TAll[1,i])
			OptVolArray = np.append(OptVolArray, OptimalVolume)
			OptGibbsArray = np.append(OptGibbsArray, OptimalGibbs)
			OptFvibArray = np.append(OptFvibArray, OptimalFvib)
			OptSvibArray = np.append(OptSvibArray, OptimalSvib)
			OptHvibArray = np.append(OptHvibArray, OptimalHvib)
			OptElArray = np.append(OptElArray, OptimalEl)

		AllData = pd.DataFrame({'Temperature': np.array(TArray),
							'Volume': np.array(OptVolArray),
							'Gibbs': np.array(OptGibbsArray),
							'Helmholtz': np.array(OptFvibArray),
							'Entropy': np.array(OptSvibArray),
							'Enthalpy': np.array(OptHvibArray),
							'Electronic_E' : np.array(OptElArray)})
	return AllData	
	
def PhaseTrans(AllData1,AllData2,Press):
	
	NumPress = np.size(Press)
	NumTemps = int(np.size(AllData1['Temperature'][:]) / np.size(Press))

	T_Trans_Array = []
	S_Trans_Array = []
	H_Trans_Array = []

	for i in range(NumPress):

		GibbsDiff = AllData2['Gibbs'][i * NumTemps: (i + 1) * NumTemps] - AllData1['Gibbs'][i * NumTemps: (i + 1) * NumTemps]
		Temps = AllData2['Temperature'][i * NumTemps: (i + 1) * NumTemps]

		T_fit = np.polyfit(Temps,GibbsDiff,3)
		TempTrans = np.roots(T_fit)[2]

		Enthalpy1 = AllData1['Electronic_E'][i * NumTemps: (i + 1) * NumTemps] + AllData1['Enthalpy'][i * NumTemps: (i + 1) * NumTemps]
		Poly1_Sfit = np.polyfit(Temps, AllData1['Entropy'][i * NumTemps: (i + 1) * NumTemps], 2)
		Poly1_Hfit = np.polyfit(Temps, Enthalpy1, 3)

		TransitionEntropy1 = np.polyval(Poly1_Sfit,TempTrans)
		TransitionEnthalpy1 = np.polyval(Poly1_Hfit,TempTrans)

		Enthalpy2 = AllData2['Electronic_E'][i * NumTemps: (i + 1) * NumTemps] + AllData2['Enthalpy'][i * NumTemps: (i + 1) * NumTemps]
		Poly2_Sfit = np.polyfit(Temps, AllData2['Entropy'][i * NumTemps: (i + 1) * NumTemps], 2)
		Poly2_Hfit = np.polyfit(Temps, Enthalpy2,3)

		TransitionEntropy2 = np.polyval(Poly2_Sfit,TempTrans)
		TransitionEnthalpy2 = np.polyval(Poly2_Hfit,TempTrans)
	
		dS = TransitionEntropy2 - TransitionEntropy1
		dH = TransitionEnthalpy2 - TransitionEnthalpy1

		T_Trans_Array = np.append(T_Trans_Array, TempTrans.real)
		S_Trans_Array = np.append(S_Trans_Array, dS.real)
		H_Trans_Array = np.append(H_Trans_Array, dH.real)

	Phase_Trans_Data = pd.DataFrame({'T_Trans': np.array(T_Trans_Array),
								 'H_Trans': np.array(H_Trans_Array),
								 'S_Trans': np.array(S_Trans_Array),
								 'Pressure': Press })

	print(Phase_Trans_Data)
	return Phase_Trans_Data



