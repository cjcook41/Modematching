import numpy as np
from . import GetDOS as dos
from scipy import optimize
import pandas as pd

def QHA(AllFreqs,FreqVols,ev_curve,count,natoms,Press):
	print('STARTING QHA')
	Na = 6.0221367e23
	natoms = sum(natoms)
	dim = 3 * natoms
	AllFreqs = AllFreqs.transpose()

#Curve Data
	Volumes = ev_curve[:,0]
	Energies = ev_curve[:,1]
	Vmin = min(Volumes)
	Vmax = max(Volumes)

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

	for i in range(np.size(AllFreqs,1)):
		freqset = np.array(AllFreqs[:,i]).reshape(-1,dim)
		[Fvib, Hvib, Svib,T] = dos.EvaluateDOS(freqset,natoms)
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
	#Loop over all T's evaluatied in Thermo (Columns in FvibAll)
		for i in range(ThermoDim):
			F_polyfit = np.polyfit(FreqVols, FvibAll[:,i], 2)
			H_polyfit = np.polyfit(FreqVols, HvibAll[:,i], 2)
			S_polyfit = np.polyfit(FreqVols, SvibAll[:,i], 2)
			PolyFvibFunction = np.polyval(F_polyfit,Volumes)
			PV = P * Volumes * 1.0e-24 * Na
			Gibbs = PolyFvibFunction + Energies + PV
			G_polyfit = np.polyfit(Volumes,Gibbs,4)

		# EP CHECKER Don't allow the minimum to be on the endpoints
			MinIndex = np.argmin(Gibbs)
			if Volumes[MinIndex] == Volumes[-1]:
				MinIndex = MinIndex - 1
			elif Volumes[MinIndex] == Volumes[0]:
				MinIndex = MinIndex + 1

		#Define initial Gibbs curve \\ G(V) 
			def PolyGFxn(V):
				return np.polyval(G_polyfit,V)
		
			##Get initial V for minimumization
			VGrid = np.arange(Vmin,Vmax,0.01)
			T_fit = np.polyval(G_polyfit,VGrid)
			initV = VGrid[np.argmin(T_fit)]

			opt_params = optimize.minimize(PolyGFxn,initV) #Provides initial Vmin and Gmin for Murnaghan Fit
			minimum = opt_params.x[0]
			minGibbs = opt_params.fun

		#Murnaghan EOS \\ Check parentheses
			def Murnaghan(x,c0,c1):
				return x[:,2] + c0*x[:,1] * (1 / (c1*(c1 - 1)) * np.power((x[:,0]/x[:,1]),1-c1) + 1 / c1 * (x[:,0]/x[:,1]) - 1 / (c1-1))

		#Figure out how to combine these fxns... it was an indexing issue
			def DoubleMurn(x,Vmin,Gmin,comp,exp): 
	
				def Murnaghan(x,c0,c1):
					return  x[2] + c0*x[1] * (1 / (c1*(c1 - 1)) * np.power((x[0]/x[1]),1-c1) + 1 / c1 * (x[0]/x[1]) - 1 / (c1-1))
				
				Test = np.zeros((1,3))
				Test[:,0] = x
				Test[:,1] = Vmin
				Test[:,2] = Gmin
				out = []

				if x > Vmin:
					ans = Murnaghan(Test[0,:],exp[0],exp[1])
					out = np.append(out,ans)
				else:
					ans = Murnaghan(Test[0,:],comp[0],comp[1])
					out = np.append(out,ans)
				return out

		#Define Compression and Expansion branches for Double Murn
			if P < 0.75:
				VolComp = Volumes[MinIndex - 1:-1]
				VolExp = Volumes[0:MinIndex  + 6]
				GibbsComp = Gibbs[MinIndex - 1:-1]
				GibbsExp = Gibbs[0:MinIndex  + 6]

			else:
				VolComp = Volumes[MinIndex - 3:-1]
				VolExp = Volumes[0:MinIndex  + 1]
				GibbsComp = Gibbs[MinIndex - 3:-1]
				GibbsExp = Gibbs[0:MinIndex  + 1]
		
		#Sets up weights for the nonlinear fit applied to Murn (Unnecessary?) BREAKS CODE; OVERWRITES VOLS
		#Wcomp = VolComp
		#for i in range(len(VolComp)):
		#	if abs(Wcomp[i] - Gmin_Volume) == 0:
		#		Wcomp[i] = 1.0
		#	elif abs(Wcomp[i] - Gmin_Volume) < 50:
		#		Wcomp[i] = 0.25
			#else:
				#Wcomp[i]
		#Wexp = VolExp
		#for i in range(len(VolExp)):
		#	if abs(Wexp[i] - Gmin_Volume) == 0:
		#		Wexp[i] = 1.0
		#	elif abs(Wexp[i] - Gmin_Volume) < 50:
		#		Wexp[i] = 0.25
		#	else:
		#		Wexp[i] = 0.1

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

			#OptimalVolume = optimize.brent(DoubleMurn,args=(minimum,minGibbs,Comp_coeffs,Exp_coeffs),brack=(Vmin,Vmax))
			#OptimalGibbs = DoubleMurn(OptimalVolume,minimum,minGibbs,Comp_coeffs,Exp_coeffs)
			OptimalParams = optimize.minimize(DoubleMurn,minimum,args=(minimum,minGibbs,Comp_coeffs,Exp_coeffs))
			OptimalVolume = OptimalParams.x[0]
			OptimalGibbs = OptimalParams.fun
			OptimalFvib = np.polyval(F_polyfit,OptimalVolume)
			OptimalSvib = np.polyval(S_polyfit,OptimalVolume)
			OptimalHvib = np.polyval(H_polyfit,OptimalVolume)
			#OptimalPV = P * OptimalVolume * 1E-24 * Na
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
		

		T_Trans_Array = np.append(T_Trans_Array, TempTrans)
		S_Trans_Array = np.append(S_Trans_Array, dS)
		H_Trans_Array = np.append(H_Trans_Array, dH)

	Phase_Trans_Data = pd.DataFrame({'T_Trans': np.array(T_Trans_Array),
								 'H_Trans': np.array(H_Trans_Array),
								 'S_Trans': np.array(S_Trans_Array),
								 'Pressure': Press })

	print(Phase_Trans_Data)
	return Phase_Trans_Data



