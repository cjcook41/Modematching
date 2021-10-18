from . import ReadData as data
import numpy as np

def ECAcoustics(atom_IDs,SSFreqFile,ECs,lattice):
		dim = sum(atom_IDs) * 3
		bohr2ang = 0.529177 ## YAML files start in bohr
		ss_params = data.ReadYaml(SSFreqFile)
		mesh = ss_params.get_SS()[1]
		lattice = lattice.reshape([3,3]) * bohr2ang

		##IN ANGSTROM
		a = lattice[0,:]
		b = lattice[1,:]
		c = lattice[2,:]

		v_cross = np.cross(b,c)
		v_dot = np.dot(a,v_cross)
		volume = v_dot * np.power(1E-10,3) 

		#Get the reciprocal Lattice
		recip_lattice = []
		astar = (2*np.pi) / v_dot * np.cross(b,c)
		bstar = (2*np.pi) / v_dot * np.cross(a,c)
		cstar = (2*np.pi) / v_dot * np.cross(a,b)
		recip_lattice.append(astar)
		recip_lattice.append(bstar)
		recip_lattice.append(cstar)
		recip_lattice = np.reshape(np.array(recip_lattice),([3,3])) #2pi / ang

		num_C = atom_IDs[0]
		num_H = atom_IDs[1]
		num_O = atom_IDs[2]
		num_N = atom_IDs[3]
		num_S = atom_IDs[4]

		mass_C = 12.011
		mass_H = 1.008
		mass_O = 15.999
		mass_N = 14.007
		mass_S = 32.06

		##Some constants
		mass_kg = 1.66054e-27
		THz = 1e12
		bar = 1000
		Pa = 100000
		wvnum = 33.35641

		mass_tot = (num_C * mass_C + num_H * mass_H + num_O * mass_O + num_N * mass_N + num_S * mass_S) * mass_kg #In kg
		density = mass_tot / volume
		
		C = np.reshape(ECs, (6,6))
		C = C * bar * Pa # Elastic constant matrix In Pa

		Christoffel_Mat = np.zeros([3,3])

		kpt_sampling_directions = np.array([0.5,0,0,
											0,0.5,0,
											0,0,0.5,
											0.5,0.5,0,
											0.5, 0, 0.5,
											0,0.5,0.5,
											0.5, 0.5, 0.5,
											0.5,-0.5,0,
											0.5,0,-0.5,
											0,0.5,-0.5,
											0.5,0.5,-0.5,
											0.5,-0.5,0.5,
											-0.5,0.5,0.5]).reshape([-1,3])

		Wmax = []
		for i in range(int(np.size(kpt_sampling_directions)/3)):
			d_cos = np.zeros(3)
			recip_sampling = np.matmul(kpt_sampling_directions[i],recip_lattice)

			#Take direction cosines wrt recip lattice system
			d_cos[0] = (recip_sampling[0] / np.linalg.norm(recip_sampling))
			d_cos[1] = (recip_sampling[1] / np.linalg.norm(recip_sampling))
			d_cos[2] = (recip_sampling[2] / np.linalg.norm(recip_sampling))

			#Building the Christoffel matrix, \Gamma_ik
			A_1 = np.square(d_cos[0])*C[0,0] + np.square(d_cos[1])*C[5,5] + np.square(d_cos[2])*C[4,4] + 2*d_cos[1]*d_cos[2]*C[4,5] + 2*d_cos[2]*d_cos[0]*C[0,4] + 2*d_cos[0]*d_cos[1]*C[0,5]
			A_2 = np.square(d_cos[0])*C[5,5] + np.square(d_cos[1])*C[1,1] + np.square(d_cos[2])*C[3,3] + 2*d_cos[1]*d_cos[2]*C[1,3] + 2*d_cos[2]*d_cos[0]*C[3,5] + 2*d_cos[0]*d_cos[1]*C[1,5]
			A_3 = np.square(d_cos[0])*C[4,4] + np.square(d_cos[1])*C[3,3] + np.square(d_cos[2])*C[2,2] + 2*d_cos[1]*d_cos[2]*C[2,3] + 2*d_cos[2]*d_cos[0]*C[2,4] + 2*d_cos[0]*d_cos[1]*C[3,4]
			alpha_23 = np.square(d_cos[0])*C[4,5] + np.square(d_cos[1])*C[1,3] + np.square(d_cos[2])*C[2,3] + 2*d_cos[1]*d_cos[2]*(1/2*(C[1,2] + C[3,3])) + 2*d_cos[2]*d_cos[0]*(1/2*(C[2,5]+C[3,4])) + 2*d_cos[0]*d_cos[1]*(1/2*(C[1,4]+C[3,5]))
			alpha_13 = np.square(d_cos[0])*C[0,4] + np.square(d_cos[1])*C[3,5] + np.square(d_cos[2])*C[2,4] + 2*d_cos[1]*d_cos[2]*(1/2*(C[2,5] + C[3,4])) + 2*d_cos[2]*d_cos[0]*(1/2*(C[0,2]+C[4,4])) + 2*d_cos[0]*d_cos[1]*(1/2*(C[0,3]+C[4,5]))
			alpha_12 = np.square(d_cos[0])*C[0,5] + np.square(d_cos[1])*C[1,5] + np.square(d_cos[2])*C[3,4] + 2*d_cos[1]*d_cos[2]*(1/2*(C[1,4] + C[3,5])) + 2*d_cos[2]*d_cos[0]*(1/2*(C[0,3]+C[4,5])) + 2*d_cos[0]*d_cos[1]*(1/2*(C[0,1]+C[5,5]))
			
			#Christoffel_Mat = np.zeros([3,3])
			Christoffel_Mat[0,0] = A_1
			Christoffel_Mat[0,1] = alpha_12
			Christoffel_Mat[0,2] = alpha_13
			Christoffel_Mat[1,0] = alpha_12
			Christoffel_Mat[1,1] = A_2
			Christoffel_Mat[1,2] = alpha_23
			Christoffel_Mat[2,0] = alpha_13
			Christoffel_Mat[2,1] = alpha_23
			Christoffel_Mat[2,2] = A_3
			
			PhaseVelocity = np.linalg.eig(Christoffel_Mat)[0] # = pv^2
			PhaseVelocity = abs(PhaseVelocity)
			kzb = np.linalg.norm(np.matmul(kpt_sampling_directions[i],recip_lattice)) * 1E10 / (2*np.pi) #kzb In 1 / m

			PhaseVelocity = np.transpose(PhaseVelocity)
			PhaseVelocity = PhaseVelocity / density # -- kg/(m*s^2) / kg/m^3 -> m^2/s^2
			PhaseVelocity = np.sqrt(PhaseVelocity) # -- m*s^-1 -> Velocity
			nan_ind = np.isnan(PhaseVelocity)
			PhaseVelocity[nan_ind] = 0

			Coeff = 2 * PhaseVelocity * kzb / np.pi # -- m/s * 1/m -> s^-1 = Frequency
			Coeff = Coeff / THz * wvnum # -- Convert to Wavenumber
			Wmax = np.append(Wmax,Coeff)

		Wmax = Wmax.reshape([-1,3])
		#print(Wmax)

		#Assign kpt direction in mesh sampling + Smoothing of the bands at kpt boundaries
		all_ids = []
		for i in range(np.size(mesh[:,0]) - 1):
			x1 = mesh[i,:]
			x2 = mesh[i+1,:]
			if np.array_equal(x1,x2):
				x2 = mesh[i-1,:]
			direction = x2 - x1
			dir_id = np.nonzero(direction)[0]
			direction2 = np.zeros(3)
			direction2[dir_id] = 0.5 #Edge of BZ to search for 13 "Base" kpt directions
			if np.linalg.norm(x2) == 0 and np.linalg.norm(x1) == 0:
				all_ids = np.append(all_ids,all_ids[-1])
			
			for j in range(np.size(kpt_sampling_directions[:,0])):
				bool_array = []
				for k in range(3):
					test = np.in1d(kpt_sampling_directions[j,k],direction2[k])
					bool_array = np.append(bool_array,test)
				if bool_array.all() == 1:
					all_ids = np.append(all_ids, j)


		Wmax_directions = []
		for i in range(np.size(all_ids) - 1):
			if i != 1:
				checkBack = all_ids[i-1]
				checkFront = all_ids[i+1]
				if all_ids[i] != checkBack and checkFront:
					all_ids[i] = checkFront
			Wmax_directions = np.append(Wmax_directions,Wmax[int(all_ids[i]),:])
		Wmax_directions = (Wmax_directions.reshape([-1,3]))

		#Nonzeros of directional Wmax
		#branch1_avg = np.average(Wmax_directions[np.nonzero(Wmax_directions[:,0])])
		#branch2_avg = np.average(Wmax_directions[np.nonzero(Wmax_directions[:,1])])
		#branch3_avg = np.average(Wmax_directions[np.nonzero(Wmax_directions[:,2])])

		#Wmax_directions[Wmax_directions[:,0]==0] = branch1_avg
		#Wmax_directions[Wmax_directions[:,1]==0] = branch2_avg
		#Wmax_directions[Wmax_directions[:,2]==0] = branch3_avg

		Wmax2 = np.array((np.average(Wmax_directions[:,0]), np.average(Wmax_directions[:,1]), np.average(Wmax_directions[:,2])))

		#All of the base kpt directions
		#branch1 = Wmax[:,0]
		#branch2 = Wmax[:,1]
		#branch3 = Wmax[:,2]
		#Wmax2 = np.array((np.average(branch1[np.nonzero(branch1)]),np.average(branch2[np.nonzero(branch2)]),np.average(branch3[np.nonzero(branch3)])))
		#Wmax2 = np.array((np.average(branch1),np.average(branch2),np.average(branch3)))

		#Add ID for endpoint
		all_ids = np.append(all_ids, all_ids[-1])
		TotKpts = np.size(all_ids)

		ss_freqs = ss_params.get_SS()[0]
		ss_freqs = np.reshape(ss_freqs,(-1, dim))
		DispFreqs = []
		mesh = ss_params.get_SS()[1]

		#print(Wmax,kpt_sampling_directions)
		for i in range(TotKpts):
			x1 = mesh[i,:]
			dir_id = int(all_ids[i])
			ratio = np.linalg.norm(x1) / np.linalg.norm(kpt_sampling_directions[dir_id,:])
			#if ratio > 1:
			#	ratio = 1
			y =  Wmax2 * abs(np.sin((ratio) * (np.pi / 2))) # -- x1 / norm(ID) both in recip space - cancel out. Day Thesis eq. 4.43
			DispFreqs = np.append(DispFreqs,y)

		DispFreqs = DispFreqs.reshape([-1,3])
		ss_freqs[:,0] = DispFreqs[:,0] / wvnum
		ss_freqs[:,1] = DispFreqs[:,1] / wvnum
		ss_freqs[:,2] = DispFreqs[:,2] / wvnum

		TEST_TAG = SSFreqFile.replace('.ss.yaml','-acoustics.csv')
		with open(TEST_TAG,'w') as f:
			for line in DispFreqs:
				f.write('{:f},{:f},{:f}\n'.format(line[0],line[1],line[2]))


		return ss_freqs

