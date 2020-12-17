from . import ReadData as data
import numpy as np
import warnings 

def get_dots(dim,shift,ref):
	dotprods = np.array([])
	for i in range(0,dim):
		for j in range(0,dim):
			val = np.dot(ref[i],shift[j])
			val = abs(val)
			dotprods = np.append(dotprods,val)
	dotprods = np.reshape(dotprods,(dim,dim))
	dotprods = np.transpose(dotprods)
	return dotprods

def get_match(dim, dotprods):
	match_array = np.array([])
	for i in range(0, dim):
		ind_max = np.argmax(dotprods[:,i])
		val_max = np.max(dotprods[:,i])
		match = (ind_max,val_max)
		match_array = np.append(match_array,match)
	match_array = np.reshape(match_array,(dim,2))
	match_array = match_array[:,0]
	match_array = match_array.astype(int)
	return match_array

def Match(atom_ID,ss_freqs,RefFile,ShiftFile):
	warnings.filterwarnings("ignore")

	#All Atom Info
	natoms = sum(atom_ID)
	dim = natoms * 3
	ss_freqs = np.reshape(ss_freqs,(-1,dim))

	Ref_UCParams = data.ReadYaml(RefFile)
	[ref_freqs, ref_vecs] = Ref_UCParams.get_UC()
	ref_vecs = np.reshape(ref_vecs,(-1,dim))

	Shift_UCParams = data.ReadYaml(ShiftFile)
	[shift_freqs, shift_vecs] = Shift_UCParams.get_UC()
	shift_vecs = np.reshape(shift_vecs,(-1,dim))

	#Get initial dotprod matrix
	dotprods = get_dots(dim,shift_vecs,ref_vecs)
	match_array = get_match(dim,dotprods)

	seen = set()
	uniq = []
	#Loops until no diplicate matches are found
	for i in match_array:
		if i not in seen:
			uniq.append(i)
			seen.add(i)
		else:
			i = (int(i))
			dupe_index = np.where(match_array == i)
			dupe_index = dupe_index[0][0]
			if dotprods[i,i] > dotprods[i,dupe_index]:
				dotprods[i,dupe_index] = 0
			if dotprods[i,i] < dotprods[i,dupe_index]:
				dotprods[i,i] = 0
			else:
				dotprods[i,dupe_index] = 0 
			match_array = get_match(dim,dotprods)

	#Set new order according to shift
	shift_freqs = shift_freqs[match_array]

	#Apply shifts to dftb freqs
	shifts = [ref_freqs - shift_freqs]
	shifts = np.array(shifts).reshape(dim,-1)

	ss_freqs = ss_freqs.transpose() 
	shifted_freqs = ss_freqs[match_array] + shifts
	shifted_freqs = shifted_freqs.flatten()
	

	
	
	return shifted_freqs

