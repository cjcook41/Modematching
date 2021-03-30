import numpy as np
import re
from itertools import islice


class ReadIFile:
	def __init__(self,f):
		self.infile = f
	def get_data(self):
		file = open(self.infile,"r")
		regexStart = re.compile('\$')
		regexEnd = re.compile('\*')
		P1atom_list = []
		P2atom_list = []
		xtal_list = []
		p1vol_list = []
		p2vol_list = []
		acoustic_list  =[]
		for line in file:
			if(regexStart.search(line) != None):
				Label = (''.join(e for e in line if e.isalnum())).lower()
			elif(regexEnd.search(line) != None):
				continue
			elif(line == "\n"):
				Label = 'Blank'
			else:
				if Label == 'P1Atoms'.lower():
					P1atom_list = np.append(P1atom_list,line)
				elif Label == 'P2Atoms'.lower():
					P2atom_list = np.append(P2atom_list,line)
				elif Label == 'Crystals'.lower():
					xtal_list = np.append(xtal_list,line)
				elif Label == 'Acoustics'.lower():
					acoustic_list = np.append(acoustic_list,line)
				elif Label == 'P1Volumes'.lower():
					p1vol_list = np.append(p1vol_list,line)
				elif Label == 'P2Volumes'.lower():
					p2vol_list = np.append(p2vol_list,line)
				elif Label == 'Pmax'.lower():
					Pmax = float(line)
				elif Label == 'Blank':
					continue
		file.close()
		def AtomIDS(atom_list):
			atom_array = []
			for x in range(np.size(atom_list)):
				test = atom_list[x].split()
				if test[0] == 'C':
					num_c = int(test[1])
				if test[0] == 'H':			
					num_h = int(test[1])
				if test[0] == 'O':
					num_o = int(test[1])
				if test[0] == 'N':
					num_n = int(test[1])
				if test[0] == 'S':
					num_s = int(test[1])
			atom_array = [num_c, num_h, num_o, num_n, num_s]
			return atom_array
		
		def AcousticJobID(acoustic_list):
			for x in range(np.size(acoustic_list)):
				test = acoustic_list[x].split()
				if test[0] == 'EC':
					EC_ID = test[1]
				if test[0] == 'Christoffel':
					Chris_ID = test[1]
			acoustic_job = [EC_ID, Chris_ID]
			return acoustic_job
		
		P1atom_array = AtomIDS(P1atom_list)
		if len(P2atom_list) == 0:
			P2atom_array = 0
		else:
			P2atom_array = AtomIDS(P2atom_list)	
		acoustic_job = AcousticJobID(acoustic_list)
	
		return P1atom_array, P2atom_array, xtal_list, acoustic_job, p1vol_list, p2vol_list, Pmax

class ReadYaml:
	def __init__(self, f):
		self.infile = f
	def get_lattice(self):
		file = open(self.infile, "r")
		LattTag = 'lattice:'
		for line in file:
			if line.rstrip()== LattTag:
				latt = (list(islice(file,0,3)))
				break
		lattice = []
		for i in range(np.size(latt)):
			tmp = (latt[i].replace('- [','').replace(',','').replace(']','')).split()
			tmp = np.double(tmp[0:3])
			lattice = np.append(lattice,tmp)
		file.close()
		return lattice

	def get_SS(self):
		file = open(self.infile, "r")
		FreqTag = 'frequency:'
		KptTag = '- q-position:'
		Frequencies = []
		Kpoints = []
		for line in file:
			line = list(line)
			line = ''.join(line)
			try:
				(line.index(FreqTag))
			except ValueError:
				pass
			else: 
				line = (line.replace(FreqTag,""))
				Frequencies.append(line)
			try:
				(line.index(KptTag))
			except ValueError:
				pass
			else:
				line = (line.replace(KptTag,""))
				line =line.replace('[','').replace(',','').replace(']','')
				line = line.split(' ')
				line = [ele for ele in line if ele][0:3]
				Kpoints.append(line)
		Kpoints = np.array(Kpoints)
		Kpoints = np.double(Kpoints)
		Frequencies = np.array(Frequencies)
		Frequencies = ''.join(Frequencies)
		Frequencies = np.fromstring(Frequencies, dtype = float,sep = '\n')
		file.close()
		return Frequencies, Kpoints
	
	def get_UC(self):
		file = open(self.infile, "r")
		FreqTag = 'frequency:'
		AtomTag = '# atom'
		Frequencies = []
		Eigenvectors = []
		for line in file:
			line = list(line)
			line = ''.join(line)
			try:
				(line.index(FreqTag))
			except ValueError:
				pass
			else: 
				line = (line.replace(FreqTag,""))
				Frequencies.append(line)
			try:
				(line.index(AtomTag))
			except ValueError:
				pass
			else:
				vector = []
				evec = list(islice(file,3))
				evec = ''.join(evec).split('\n')
				for x in evec:
					tmp = x.replace('- [','').replace(',','').replace(']','')
					tmp = tmp.split(' ')
					tmp = [ele for ele in tmp if ele]
					if not tmp:
						continue
					else:
						evec2 = tmp[0]
					vector.append(evec2)
				test = np.array(vector)
				Eigenvectors = np.append(Eigenvectors, test)
		Eigenvectors = Eigenvectors.reshape([-1,3])
		Eigenvectors = np.double(Eigenvectors)
		Frequencies = np.array(Frequencies)
		Frequencies = ''.join(Frequencies)
		Frequencies = np.fromstring(Frequencies, dtype = float,sep = '\n')
		file.close()
		return Frequencies, Eigenvectors
	
class ReadElastic:
	def __init__(self, t):
		self.tensor = t
	def get_tensor(self):
		file = open(self.tensor,"r")
		tens_array = []
		for line in file:
			line = line.split()
			tens_array.append(line)
		tens_array = np.array(tens_array)
		tens_array = tens_array.astype(float)
		file.close()
		return tens_array

class ReadLattice:
	def __init__(self,l):
		self.lattice = l
	def get_vecs(self):
		file = open(self.lattice)
		file = list(file)
		vecs = []
		for line in file:
			line = line.split()
			vecs.append(line)
		vecs = np.array(vecs)
		vecs = vecs.astype(float)
		file.close()
		return vecs
class Quasiharmonic:
	def __init__(self, v):
		self.evcurve = v
	def get_evcurve(self):
		file = open(self.evcurve,"r")
		ev_data = list(file)
		ev_data = [x for x in ev_data if "#" not in x]
		ev_data = ''.join(ev_data)
		ev_data = np.fromstring(ev_data, dtype = float, sep = '\n')
		ev_data = np.reshape(ev_data,(-1,2))
		file.close()
		return ev_data

