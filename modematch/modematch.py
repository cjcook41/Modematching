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

def get_match(available_matches, dotprods):
	match_array = np.array([])
	for i in available_matches:
		ind_max = np.argmax(dotprods[:,i])
		match_array = np.append(match_array,ind_max)
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
	available_matches = list(range(0,dim))
	match_array = get_match(available_matches,dotprods)

	seen = set()
	uniq = []
	#Loops until no diplicate matches are found
	for i in range(np.size(match_array)):
		num = match_array[i]
		if num not in seen:
			available_matches.remove(num)
			uniq.append(num)
			seen.add(num)
		else:
			olaps = []
			dupe_index = np.where(match_array == num)
			dupe_index = dupe_index[0][0]
			if dotprods[num,dupe_index] > dotprods[num,i]:
				dotprods[num,i] = 0
				for j in available_matches:
					olaps = np.append(olaps,dotprods[dupe_index,j])
				newmatch = np.argmax(olaps)
				newmatch = available_matches[newmatch]
				match_array[i] = newmatch
			if dotprods[num,dupe_index] < dotprods[num,i]:
				dotprods[num,dupe_index] = 0
				for j in available_matches:
					olaps = np.append(olaps,dotprods[dupe_index,j])
				newmatch = np.argmax(olaps)
				newmatch = available_matches[newmatch]
				match_array[dupe_index] = newmatch

	#Set new order according to shift
	shift_freqs = shift_freqs[match_array]

	#Apply shifts to dftb freqs
	shifts = [ref_freqs - shift_freqs]
	shifts = np.array(shifts).reshape(dim,-1)

	ss_freqs = ss_freqs.transpose() 
	shifted_freqs = ss_freqs[match_array] + shifts
	shifted_freqs = shifted_freqs.flatten()
	

	
	
	return shifted_freqs

