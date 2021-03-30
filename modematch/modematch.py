from . import ReadData as data
import numpy as np
import warnings
from matching.games import StableMarriage

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

def Match(atom_ID,ss_freqs,RefFile,ShiftFile):
	warnings.filterwarnings("ignore")

	#All Atom Info
	natoms = sum(atom_ID)
	dim = natoms * 3
	ss_freqs = np.reshape(ss_freqs,(-1,dim))

	Ref_UCParams = data.ReadYaml(RefFile)
	[ref_freqs, ref_vecs] = Ref_UCParams.get_UC()
	ref_vecs = np.reshape(ref_vecs,(-1,dim))
	ref_freqs[0] = 0
	ref_freqs[1] = 0
	ref_freqs[2] = 0

	Shift_UCParams = data.ReadYaml(ShiftFile)
	[shift_freqs, shift_vecs] = Shift_UCParams.get_UC()
	shift_vecs = np.reshape(shift_vecs,(-1,dim))

	#Get initial dotprod matrix
	dots = get_dots(dim,shift_vecs,ref_vecs)

	dft_ranks = [] #Ranking DFT matches with FIXED DFTB modes
	for i in range(dim):
		sort_dots = np.sort(dots[i,:])[::-1]
		seen = set()
		available_matches = list(range(0,dim))
		for j in range(dim):
			val = sort_dots[j]
			arr = dots[i,:]
			tmp = (np.where(arr == val)[0][0])
			if tmp not in seen:
				available_matches.remove(tmp)
				seen.add(tmp)
			else:
				availdots = []
				for k in available_matches:
					availdots.append(dots[i,k])
				tmp = available_matches[availdots.index(max(availdots))]
				available_matches.remove(tmp)
			dft_ranks = np.append(dft_ranks,tmp)

	dftb_ranks = [] #Ranking DFTB matches with FIXED DFT modes
	for i in range(dim):
		sort_dots = np.sort(dots[:,i])[::-1]
		seen = set()
		available_matches = list(range(0,dim))
		for j in range(dim):
			val = sort_dots[j]
			arr = dots[:,i]
			tmp = (np.where(arr == val)[0][0])
			if tmp not in seen:
				available_matches.remove(tmp)
				seen.add(tmp)
			else:
				availdots = []
				for k in available_matches:
					availdots.append(dots[k,i])
				tmp = available_matches[availdots.index(max(availdots))]
				available_matches.remove(tmp)
			dftb_ranks = np.append(dftb_ranks,tmp)

	dft_ranks = dft_ranks.reshape([-1,dim])
	dftb_ranks = dftb_ranks.reshape([-1,dim])

	dftRankDict = dict(((i), dft_ranks[i][:]) for i in range(dim)) #Favored DFT pairing to a DFTB mode -- "Suitor"
	dftbRankDict = dict(((i), dftb_ranks[i][:]) for i in range(dim)) #Favored DFTB pairing to a DFT mode -- "Reviewer"
	
	matching = StableMarriage.create_from_dictionaries(dftRankDict, dftbRankDict)
	AllMatches = matching.solve()
	AllMatches = list(str(AllMatches).replace('{','').replace('}','').split(','))
	matches = []

	for i in AllMatches:
		line = (i.split(':'))
		DFTmode  = int(line[0])
		DFTBmode = int(line[1])
		matches.append(DFTmode)
		matches.append(DFTBmode)

	matches = np.array(matches).reshape([dim,2])

	##In matches array, the first column is the DFT modes, and the second is DFTB. Need to re-sort in order to get the
	##DFTB modes in ascending order for the shift
	match_array = matches[np.argsort(matches[:,1])]
	match_array = match_array[:,0]
	#Set new order according to match
	shift_freqs = shift_freqs[match_array]

	#Apply shifts to dftb freqs
	shifts = [ref_freqs - shift_freqs]
	shifts = np.array(shifts).reshape(dim,-1)

	ss_freqs = ss_freqs.transpose() 
	shifted_freqs = ss_freqs[match_array] + shifts
	shifted_freqs = shifted_freqs.flatten()
	

	
	
	return shifted_freqs

