from . import ReadData as data
from . import JobControl 
import numpy as np
import sys

def RunJob():


	print('Modematch version 0.0.7')
	sys.setrecursionlimit(3000) ## DEFAULT RECURSION LIM HIT WITH STABLE MARRIAGE MATCHING ALG


	job_description = data.ReadIFile('infile') #Chnge this when im done(?)
	job_params = job_description.get_data()

	P1atom_IDs = job_params[0]
	P2atom_IDs = job_params[1]
	struct_IDs = job_params[2]
	acoustic_ID = job_params[3]
	p1freq_vols = job_params[4]
	p2freq_vols = job_params[5]
	Pmax = job_params[6]

	Press = np.arange(0,Pmax+0.01,0.01)

	if P2atom_IDs == 0:
		PhaseTrans = False
		print('Only 1 structure, no phase trans...')
		natoms = sum(P1atom_IDs)
		#dim = 3 * natoms
		poly1 = struct_IDs[0].rstrip()

	else:
		PhaseTrans = True
		print('Multiple structures. Predicting Phase Transition...')
		P1atoms = sum(P1atom_IDs)
		P2atoms = sum(P2atom_IDs)
		#P1dim = 3 * P1atoms
		#P2dim = 3 * P2atoms
		poly1 = struct_IDs[0].rstrip()
		poly2 = struct_IDs[1].rstrip()

	#Number of Fvibs to calculate for each structure
	P1Count = int(len(p1freq_vols))
	P2Count = int(len(p2freq_vols))

	if PhaseTrans == False:
		P1Control = JobControl.Control(poly1, P1Count, P1atom_IDs)
		P1Control.MakePaths()

		if acoustic_ID[0].lower() == 'true' and acoustic_ID[1].lower() == 'true':
			AllConstants = P1Control.SolveECs()
			NewAcoustics = P1Control.Dispersion(AllConstants)
			ShiftedFreqs = P1Control.ShiftPlusECs(NewAcoustics)

		elif acoustic_ID[0].lower() == 'false' and acoustic_ID[1].lower() == 'false': ##No EC acoustics at all
			ShiftedFreqs = P1Control.Modematch()[0]

		elif acoustic_ID[0].lower() == 'false' and acoustic_ID[1].lower() == 'true': ##User input EC's
			AllConstants = P1Control.ReadECs()
			NewAcoustics = P1Control.Dispersion(AllConstants)
			ShiftedFreqs = P1Control.ShiftPlusECs(NewAcoustics)

		#Now we have the shifted freqs...
		if P1Count == 1:
			print("One structure, no QHA")
			from . import GetDOS as dos
			import pandas as pd
			[Fvib, Hvib, Svib,T] = dos.EvaluateDOS(ShiftedFreqs,natoms)
			ThermoData = pd.DataFrame({'Temperature':np.array(T),
										'Enthalpy': np.array(Hvib),
										'Entropy': np.array(Svib),
										'Helmholtz': np.array(Fvib)})
		else:
			ThermoData = P1Control.QuasiHarmonic(ShiftedFreqs,p1freq_vols,Press)

		filename = poly1 + '.csv'
		ThermoData.to_csv(filename,index=False)
		print('JOB SUCCESSFUL!!')

	else:
		##FIRST STRUCTURE
		P1Control = JobControl.Control(poly1, P1Count, P1atom_IDs)
		P1Control.MakePaths()
		if acoustic_ID[0].lower() == 'true' and acoustic_ID[1].lower() == 'true':
			AllConstants = P1Control.SolveECs()
			NewAcoustics = P1Control.Dispersion(AllConstants)
			ShiftedFreqs = P1Control.ShiftPlusECs(NewAcoustics)

		elif acoustic_ID[0].lower() == 'false' and acoustic_ID[1].lower() == 'false': ##No EC acoustics at all
			ShiftedFreqs = P1Control.Modematch()[0]

		elif acoustic_ID[0].lower() == 'false' and acoustic_ID[1].lower() == 'true': ##User input EC's
			AllConstants = P1Control.ReadECs()
			NewAcoustics = P1Control.Dispersion(AllConstants)
			ShiftedFreqs = P1Control.ShiftPlusECs(NewAcoustics)
		#Now we have the shifted freqs...
		if P1Count == 1:
			P1CheckQHA = False
			print("One structure, no QHA")
			from . import GetDOS as dos
			import pandas as pd
			[Fvib, Hvib, Svib,T] = dos.EvaluateDOS(ShiftedFreqs,P1atoms)
			P1ThermoData = pd.DataFrame({'Temperature':np.array(T),
										'Enthalpy': np.array(Hvib),
										'Entropy': np.array(Svib),
										'Helmholtz': np.array(Fvib)})
		else:
			P1CheckQHA = True
			P1ThermoData = P1Control.QuasiHarmonic(ShiftedFreqs,p1freq_vols,Press)
		filename = poly1 + '.csv'
		P1ThermoData.to_csv(filename,index=False)
		print('DONE WITH ', poly1)

		##SECOND STRUCTURE
		P2Control = JobControl.Control(poly2, P2Count, P2atom_IDs)
		P2Control.MakePaths()

		#Check to see if we're calculating new acoustics / ECs
		if acoustic_ID[0].lower() == 'true' and acoustic_ID[1].lower() == 'true':
			AllConstants = P2Control.SolveECs()
			NewAcoustics = P2Control.Dispersion(AllConstants)
			ShiftedFreqs = P2Control.ShiftPlusECs(NewAcoustics)

		elif acoustic_ID[0].lower() == 'false' and acoustic_ID[1].lower() == 'false': ##No EC acoustics at all
			ShiftedFreqs = P2Control.Modematch()[0]


		elif acoustic_ID[0].lower() == 'false' and acoustic_ID[1].lower() == 'true': ##User input EC's
			AllConstants = P2Control.ReadECs()
			NewAcoustics = P2Control.Dispersion(AllConstants)
			ShiftedFreqs = P2Control.ShiftPlusECs(NewAcoustics)

		#Now we have the shifted freqs...
		if P2Count == 1:
			P2CheckQHA = False
			print("One structure, no QHA")
			import GetDOS as dos
			import pandas as pd
			[Fvib, Hvib, Svib,T] = dos.EvaluateDOS(ShiftedFreqs,P2atoms)
			P2ThermoData = pd.DataFrame({'Temperature':np.array(T),
										'Enthalpy': np.array(Hvib),
										'Entropy': np.array(Svib),
										'Helmholtz': np.array(Fvib)})
		else:
			P2CheckQHA = True
			P2ThermoData = P2Control.QuasiHarmonic(ShiftedFreqs,p2freq_vols,Press)
		filename = poly2 + '.csv'
		P2ThermoData.to_csv(filename,index=False)
		print('DONE WITH ', poly2)


		if P1CheckQHA and P2CheckQHA == True:
			PhaseTransData = P1Control.PhaseTransControl(P1ThermoData,P2ThermoData,Press)
			filename = 'Phase_Transition'
			PhaseTransData.to_csv(filename,index=False)

		print('JOB SUCCESSFUL!!')
