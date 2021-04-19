from . import ReadData as data
import numpy as np
import cmath

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

		#a_norm = np.linalg.norm(a)
		#b_norm = np.linalg.norm(b)
		#c_norm = np.linalg.norm(c)

		#beta = np.arccos(np.dot(a,c) / (a_norm * c_norm))
		#alpha = np.arccos(np.dot(b,c) / (b_norm * c_norm))
		#gamma = np.arccos(np.dot(a,b) / (a_norm * b_norm))

		v_cross = np.cross(b,c)
		v_dot = np.dot(a,v_cross)
		volume = v_dot * np.power(1E-10,3) 

		#Get the reciprocal Lattice
		recip_lattice = []
		astar = 1 / v_dot * np.cross(b,c)
		bstar = 1 / v_dot * np.cross(a,c)
		cstar = 1 / v_dot * np.cross(a,b)
		recip_lattice.append(astar)
		recip_lattice.append(bstar)
		recip_lattice.append(cstar)
		recip_lattice = np.reshape(np.array(recip_lattice),([3,3]))


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
		THz = 10e12
		bar = 1000
		Pa = 100000
		wvnum = 33.35641

		mass_tot = (num_C * mass_C + num_H * mass_H + num_O * mass_O + num_N * mass_N) * mass_kg #In kg

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
			d_cos = np.matmul(kpt_sampling_directions[i],recip_lattice)/ np.linalg.norm(np.matmul(kpt_sampling_directions[i], recip_lattice)) ##direction cosines with respect to reciprocal lattice
			
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
			kzb = np.linalg.norm(d_cos) / 1E-10 # --- Length in pi / meter -- kzb 


			PhaseVelocity = np.transpose(PhaseVelocity)
			PhaseVelocity = np.sqrt(PhaseVelocity / density)
			PhaseVelocity = np.real(PhaseVelocity)
			nan_ind = np.isnan(PhaseVelocity)
			PhaseVelocity[nan_ind] = 0

			Coeff = 2 * PhaseVelocity * kzb / np.pi
			Coeff = Coeff / THz * wvnum
			Wmax = np.append(Wmax,Coeff)

		Wmax = Wmax.reshape([-1,3])

		#Using Wmax averages... not going this route
		#a_Wmax = Wmax[:,0]
		#b_Wmax = Wmax[:,1]
		#c_Wmax = Wmax[:,2]

		#a_Wmax = a_Wmax[a_Wmax!=0]
		#b_Wmax = b_Wmax[b_Wmax!=0]
		#c_Wmax = c_Wmax[c_Wmax!=0]

		#a_Wmax = np.average(a_Wmax)
		#b_Wmax = np.average(b_Wmax)
		#c_Wmax = np.average(c_Wmax)

		#avg_Wmax = np.array([a_Wmax, b_Wmax, c_Wmax])

		#Assign kpt direction in mesh sampling + Smoothing of the bands at kpt boundaries
		all_ids = []
		for i in range(np.size(mesh[:,0]) - 1):
			x1 = mesh[i,:]
			x2 = mesh[i+1,:]
			direction = x2 - x1
			dir_id = np.nonzero(direction)[0]
			direction2 = x1
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
		Wmax2 = np.array((np.average(Wmax_directions[:,0]), np.average(Wmax_directions[:,1]), np.average(Wmax_directions[:,2])))
		
		#Add ID for endpoint
		all_ids = np.append(all_ids, all_ids[-1])

		TotKpts = np.size(all_ids)

		ss_freqs = ss_params.get_SS()[0]
		ss_freqs = np.reshape(ss_freqs,(-1, dim))
		DispFreqs = []
		mesh = ss_params.get_SS()[1]

		for i in range(TotKpts):
			x1 = mesh[i,:]
			dir_id = int(all_ids[i])
			y =  Wmax2 * abs(np.sin((np.linalg.norm(x1) / np.linalg.norm(kpt_sampling_directions[dir_id,:])) * np.pi / 2))

			DispFreqs = np.append(DispFreqs,y)

		DispFreqs = DispFreqs.reshape([-1,3])
		ss_freqs[:,0] = DispFreqs[:,0] / wvnum
		ss_freqs[:,1] = DispFreqs[:,1] / wvnum
		ss_freqs[:,2] = DispFreqs[:,2] / wvnum

		return ss_freqs

