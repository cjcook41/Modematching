import modematch
import numpy as np
import ElasticConstants as ECs
import ReadData as data
import GetDOS as dos
import qha2
import Christoffel
import os

class Control:
	def __init__(self,x,c,a):
		self.xtal = x
		self.count = c
		self.atoms = a
	def Modematch(self):
		print("Shifting ", self.xtal, " Frequencies...")
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + '/' + self.xtal + '/'
		AllFreqs = []
		for x in range(self.count):
			x = str(x+1)
			SubDir = XtalPath + x
			
			SSFreqFile = SubDir + '/' + self.xtal + x + '.ss.yaml'
			RefFreqFile = SubDir + '/' + self.xtal + x  + '.ref.yaml'
			ShiftFreqFile = SubDir + '/' + self.xtal + x + '.shift.yaml'

			Supercell_Params = data.ReadYaml(SSFreqFile)
			[ss_freqs, mesh] = Supercell_Params.get_SS()

			shifted_freqs = modematch.Match(self.atoms,ss_freqs,RefFreqFile,ShiftFreqFile)

			AllFreqs = np.append(AllFreqs, shifted_freqs)
		
		AllFreqs = np.reshape(AllFreqs,(self.count,-1))
		return AllFreqs
		
	def MakePaths(self):
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + '/' + self.xtal + '/'
		if not os.path.exists(XtalPath):
			os.mkdir(XtalPath)
		for file in os.listdir(CurrentDir):	
			if file.startswith(self.xtal):
				try:
					os.rename(file, XtalPath + file)
				except OSError:
					continue
		FileNames = []
		for x in range(self.count):
			x = str(x+1)
			SubDir = XtalPath + x
			if not os.path.exists(SubDir):
				os.mkdir(SubDir)
			SSFreqFile = XtalPath + self.xtal + x + '.ss.yaml'
			RefFreqFile = XtalPath + self.xtal + x + '.ref.yaml'
			ShiftFreqFile = XtalPath + self.xtal + x + '.shift.yaml'
			StrainFile = XtalPath + self.xtal + x + '.strains'
			StressFile = XtalPath + self.xtal + x + '.stresses'
			LattFile = XtalPath + self.xtal + x + '.lattice'

			os.rename(SSFreqFile, SubDir + '/' + self.xtal + x + '.ss.yaml')
			os.rename(RefFreqFile, SubDir + '/' + self.xtal + x + '.ref.yaml')
			os.rename(ShiftFreqFile, SubDir + '/' + self.xtal + x + '.shift.yaml')
			os.rename(StrainFile, SubDir + '/' + self.xtal + x + '.strains')
			os.rename(StressFile, SubDir + '/' + self.xtal + x + '.stresses')
			os.rename(LattFile, SubDir + '/' + self.xtal + x + '.lattice')


	
	def SolveECs(self):
		print("Solving ", self.xtal, " Elastic Constants...")
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + '/' + self.xtal + '/'
		AllConstants = []
		for x in range(self.count):
			x = str(x + 1)
			SubDir = XtalPath + x
			StrainFile = SubDir + '/' + self.xtal + x + '.strains'
			StressFile = SubDir + '/' + self.xtal + x + '.stresses'
				
			EC_Matrix = ECs.ElasticConstants(StressFile,StrainFile)
			EC_Matrix = EC_Matrix.flatten()
			AllConstants = np.append(AllConstants, EC_Matrix)
	
		AllConstants = np.reshape(AllConstants, (self.count,-1))
		return AllConstants
	
	def Dispersion(self,ECs):
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + '/' + self.xtal + '/'
		AllFreqs = []
		for x in range(self.count):	
			y = str(x + 1)
			SubDir = XtalPath + y
			LattFile = SubDir + '/' + self.xtal + y + '.lattice'
			SSFreqFile = SubDir + '/' + self.xtal + y + '.ss.yaml'

			#InputFreqs = self.AllFreqs[x]
			InputECs = ECs[x]
			NewFreqs = Christoffel.ECAcoustics(self.atoms,SSFreqFile,InputECs,LattFile)
			AllFreqs = np.append(AllFreqs, NewFreqs)
		
		AllFreqs = np.reshape(AllFreqs,(self.count,-1))
		return AllFreqs

	def ShiftPlusECs(self,Freqs):
		print("Shift + EC Correction...")
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + '/' + self.xtal + '/'
		AllFreqs = []
		for x in range(self.count):
			y = str(x + 1)
			SubDir = XtalPath + y
			RefFreqFile = SubDir + '/' + self.xtal + y  + '.ref.yaml'
			ShiftFreqFile = SubDir + '/' + self.xtal + y + '.shift.yaml'

			ss_freqs = Freqs[x,:]
			NewFreqs = modematch.Match(self.atoms,ss_freqs,RefFreqFile,ShiftFreqFile)
			AllFreqs = np.append(AllFreqs, NewFreqs)
		AllFreqs = np.reshape(AllFreqs,(self.count,-1))
		return AllFreqs

	def QuasiHarmonic(self,Freqs,Vols,Press):
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + '/' + self.xtal + '/'
		CurveFile = XtalPath + self.xtal + '.dat'
		EVParams = data.Quasiharmonic(CurveFile)
		ev_data = EVParams.get_evcurve()
		Vols = np.double(Vols)
		if np.size(Vols) == 1:
			AllData = dos.EvaluateDOS(Freqs,self.atoms)
		else:
			AllData = qha2.QHA(Freqs,Vols,ev_data,self.count,self.atoms,Press)
		return AllData
		
	def PhaseTransControl(self,AllData1,AllData2,Press):
		PhaseTransData = qha2.PhaseTrans(AllData1,AllData2,Press)
		return PhaseTransData
		
	
