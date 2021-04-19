import os

import numpy as np
from shutil import copyfile

from . import Christoffel
from . import ElasticConstants as ECs
from . import GetDOS as dos
from . import ReadData as data
from . import modematch
from . import qha2


class Control:
	def __init__(self,x,c,a):
		self.xtal = x
		self.count = c
		self.atoms = a
	def Modematch(self):
		print("Shifting ", self.xtal, " Frequencies...")
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + os.path.sep + self.xtal + os.path.sep
		AllFreqs = []
		for x in range(self.count):
			x = str(x+1)
			SubDir = XtalPath + x
			
			SSFreqFile = SubDir + os.path.sep + self.xtal + x + '.ss.yaml'
			RefFreqFile = SubDir + os.path.sep + self.xtal + x + '.ref.yaml'
			ShiftFreqFile = SubDir + os.path.sep + self.xtal + x + '.shift.yaml'

			Supercell_Params = data.ReadYaml(SSFreqFile)
			[ss_freqs, mesh] = Supercell_Params.get_SS()
			shifted_freqs  = modematch.Match(self.atoms,ss_freqs,RefFreqFile,ShiftFreqFile)
			np.savetxt(self.xtal + x + 'freqs.csv', shifted_freqs)
			AllFreqs = np.append(AllFreqs, shifted_freqs)
		
		AllFreqs = np.reshape(AllFreqs,(self.count,-1))
		return AllFreqs, mesh
		
	def MakePaths(self):
		CurrentDir = os.getcwd()
		RawfileDir = CurrentDir + os.path.sep + 'datafiles' + os.path.sep
		XtalPath = CurrentDir + os.path.sep + self.xtal + os.path.sep
		if not os.path.exists(XtalPath):
			os.mkdir(XtalPath)
		for file in os.listdir(RawfileDir):
			if file.startswith(self.xtal):
				try:
					source = RawfileDir + file
					target = XtalPath + file
					copyfile(source, target)
				except OSError as e:
					print('No File %s' % e)
					continue

		for x in range(self.count):
			x = str(x+1)
			SubDir = XtalPath + x
			if not os.path.exists(SubDir):
				os.mkdir(SubDir)
			SSFreqFile = XtalPath + self.xtal + x + '.ss.yaml'
			RefFreqFile = XtalPath + self.xtal + x + '.ref.yaml'
			ShiftFreqFile = XtalPath + self.xtal + x + '.shift.yaml'

			SSFreqFileOut = SubDir + os.path.sep + self.xtal + x + '.ss.yaml'
			RefFreqFileOut = SubDir + os.path.sep + self.xtal + x + '.ref.yaml'
			ShiftFreqFileOut = SubDir + os.path.sep + self.xtal + x + '.shift.yaml'

			try:
				os.rename(SSFreqFile, SSFreqFileOut)
			except FileExistsError:
				os.remove(SSFreqFileOut)
				os.rename(SSFreqFile, SSFreqFileOut)

			try:
				os.rename(RefFreqFile, RefFreqFileOut)
			except FileExistsError:
				os.remove(RefFreqFileOut)
				os.rename(RefFreqFile, RefFreqFileOut)

			try:
				os.rename(ShiftFreqFile, ShiftFreqFileOut)
			except FileExistsError:
				os.remove(ShiftFreqFileOut)
				os.rename(ShiftFreqFile, ShiftFreqFileOut)

	def SolveECs(self):
		print("Solving ", self.xtal, " Elastic Constants...")
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + os.path.sep + self.xtal + os.path.sep
		AllConstants = []
		for x in range(self.count):
			x = str(x + 1)
			SubDir = XtalPath + x

			StrainFile = XtalPath + self.xtal + '.strains'
			StressFile = XtalPath + self.xtal + x + '.stresses'

			#StrainFileOut = SubDir + os.path.sep + self.xtal + '.strains'
			StressFileOut = SubDir + os.path.sep + self.xtal + x + '.stresses'

			#try:
			#	os.rename(StrainFile,StrainFileOut)
			#except FileExistsError:
			#	os.remove(StrainFileOut)
			#	os.rename(StrainFile,StrainFileOut)
			try:
				os.rename(StressFile,StressFileOut)
			except FileExistsError:
				os.remove(StressFileOut)
				os.rename(StressFile,StressFileOut)

			EC_Matrix = ECs.ElasticConstants(StressFileOut,StrainFile)
			ECFile = XtalPath + self.xtal + x + '.ecs'
			np.savetxt(ECFile, EC_Matrix, delimiter='\t')
			EC_Matrix = EC_Matrix.flatten()
			AllConstants = np.append(AllConstants, EC_Matrix)
	
		AllConstants = np.reshape(AllConstants, (self.count,-1))
		return AllConstants

	def ReadECs(self):
		print("Found Elastic constants, solving dispersion...")
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + os.path.sep + self.xtal + os.path.sep
		AllConstants = []

		for x in range(self.count):
			x = str(x + 1)

			ECFile = XtalPath + self.xtal + x + '.ecs'
			ECFileOut = XtalPath + os.path.sep + x + os.path.sep + self.xtal + x + '.ecs'

			try:
				copyfile(ECFile,ECFileOut)
				os.remove(ECFile)
			except OSError as e:
				print('Cannot Find file %s' % e)
				continue

			file = open(ECFileOut, "r")
			ECs = []
			for line in file:
				ECs = np.append(ECs, np.double(line.split()))
			file.close()
			AllConstants = np.append(AllConstants,ECs)
		AllConstants = AllConstants.reshape([self.count,-1])
		return AllConstants



	def Dispersion(self,ECs):
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + os.path.sep + self.xtal + os.path.sep
		AllFreqs = []
		for x in range(self.count):	
			y = str(x + 1)
			SubDir = XtalPath + y

			SSFreqFile = SubDir + os.path.sep + self.xtal + y + '.ss.yaml'
			ShiftedFreqFile = SubDir + os.path.sep + self.xtal + y + '.freqs'
			Supercell_Params = data.ReadYaml(SSFreqFile)
			lattice = Supercell_Params.get_lattice()

			InputECs = ECs[x]
			NewFreqs = Christoffel.ECAcoustics(self.atoms,SSFreqFile,InputECs,lattice)

			np.savetxt(ShiftedFreqFile, NewFreqs)
			AllFreqs = np.append(AllFreqs, NewFreqs)
		
		AllFreqs = np.reshape(AllFreqs,(self.count,-1))
		return AllFreqs

	def ShiftPlusECs(self,Freqs):
		print("Shift + EC Correction...")
		CurrentDir = os.getcwd()  #+ os.path.sep + 'datafiles'
		XtalPath = CurrentDir + os.path.sep + self.xtal + os.path.sep
		AllFreqs = []
		for x in range(self.count):
			y = str(x + 1)
			SubDir = XtalPath + y
			RefFreqFile = SubDir + os.path.sep + self.xtal + y  + '.ref.yaml'
			ShiftFreqFile = SubDir + os.path.sep + self.xtal + y + '.shift.yaml'
			ShiftedFreqFile = SubDir + os.path.sep + self.xtal + y + '.freqs'
			ss_freqs = Freqs[x,:]
			NewFreqs = modematch.Match(self.atoms,ss_freqs,RefFreqFile,ShiftFreqFile)

			np.savetxt(ShiftedFreqFile, NewFreqs)
			AllFreqs = np.append(AllFreqs, NewFreqs)

		AllFreqs = np.reshape(AllFreqs,(self.count,-1))
		return AllFreqs

	def QuasiHarmonic(self,Freqs,Vols,Press):
		CurrentDir = os.getcwd()
		XtalPath = CurrentDir + os.path.sep + self.xtal + os.path.sep
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
		
	
