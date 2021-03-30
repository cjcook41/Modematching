import numpy as np
from . import ReadData as data
import re
from sympy import *

def ElasticConstants(stresses,strains):

	stresses = data.ReadElastic(stresses)
	strains = data.ReadElastic(strains)
	
	stresses = stresses.get_tensor()
	strains = strains.get_tensor()

	stresses = stresses * 2 # Ha to Ry

	C11, C12, C13, C14, C15, C16 = symbols('C11 C12 C13 C14 C15 C16')
	C22, C23, C24, C25, C26 = symbols('C22 C23 C24 C25 C26')
	C33, C34, C35, C36 = symbols('C33 C34 C35 C36')
	C44, C45, C46 = symbols('C44 C45 C46')
	C55, C56 = symbols('C55 C56')
	C66 = symbols('C66')

	EC_matrix = Matrix([[C11, C12, C13, C14, C15, C16],[C12, C22, C23, C24, C25, C26], [C13, C23, C33, C34, C35, C36],[C14, C24, C34, C44, C45, C46],[C15, C25, C35, C45, C55, C56],[C16, C26, C36, C46, C56, C66]])

	e_vec = []
	t_vec = []
	constants = []

#Initialize stress const arrays
	C11_ = []
	C12_ = []
	C13_ = []
	C14_ = []
	C15_ = []
	C16_ = []
	C22_ = []
	C23_ = []
	C24_ = []
	C25_ = []
	C26_ = []
	C33_ = []
	C34_ = []
	C35_ = []
	C36_ = []
	C44_ = []
	C45_ = []
	C46_ = []
	C55_ = []
	C56_ = []
	C66_ = []

	bohr2ang = 0.529117
	ang2m = 1E10
	ry2joule = 2.179872E-18
	pa2kbar = 1E8

	fit_x = []
	fit_y1 = []
	fit_y2 = []
	fit_y3 = []
	fit_y4 = []
	fit_y5 = []
	fit_y6 = []


#Get range of stress/strains to fit over

	e_range = np.arange(np.min(strains) - 0.0005,np.max(strains) + 0.0005 ,0.0005)
	t_range1 = []
	t_range2 = []
	t_range3 = []
	t_range4 = []
	t_range5 = []
	t_range6 = []

	new_strain = []
	new_stress = []

#This will loop to create unique fits based on the strain being applied
	for i in range(np.size(strains,0)):
		tmp = stresses[i,:]
		t_vec = np.append(t_vec,tmp)
		tmp2 = strains[i,:]
		e_vec = np.append(e_vec,tmp2)
		if np.mod(i+1,3) == 0:
			e_vec = np.reshape(e_vec,(3,3))
			t_vec = np.reshape(t_vec,(3,3))
			tag1 = sum(np.concatenate(np.nonzero(e_vec)))
			tag2 = np.size(np.nonzero(e_vec))
			#STRAIN 1
			if tag1 == 0 and tag2 == 2:
				fit_x = np.append(fit_x, e_vec[0,0])
				fit_y1 = np.append(fit_y1, t_vec[0,0]) #C11
				fit_y2 = np.append(fit_y2, t_vec[1,1]) #C21
				fit_y3 = np.append(fit_y3, t_vec[2,2]) #C31
				fit_y4 = np.append(fit_y4, t_vec[2,1]) #C41
				fit_y5 = np.append(fit_y5, t_vec[0,2]) #C51
				fit_y6 = np.append(fit_y6, t_vec[0,1]) #C61
				if len(fit_x) == 4:
					p1 = np.polyfit(fit_x, fit_y1,2)
					p2 = np.polyfit(fit_x, fit_y2,2)
					p3 = np.polyfit(fit_x, fit_y3,2)
					p4 = np.polyfit(fit_x, fit_y4,2)
					p5 = np.polyfit(fit_x, fit_y5,2)
					p6 = np.polyfit(fit_x, fit_y6,2)
					t_range1 = np.append(t_range1,np.polyval(p1,e_range))
					t_range2 = np.append(t_range2,np.polyval(p2,e_range))
					t_range3 = np.append(t_range3,np.polyval(p3,e_range))
					t_range4 = np.append(t_range4,np.polyval(p4,e_range))
					t_range5 = np.append(t_range5,np.polyval(p5,e_range))
					t_range6 = np.append(t_range6,np.polyval(p6,e_range))
					for i in range(np.size(e_range,0)):
						e_vec[0,0] = e_range[i]
						t_vec[0,0] = t_range1[i]
						t_vec[1,1] = t_range2[i]
						t_vec[2,2] = t_range3[i]
						t_vec[2,1] = t_range4[i]
						t_vec[1,2] = t_range4[i]
						t_vec[0,2] = t_range5[i]
						t_vec[2,0] = t_range5[i]
						t_vec[0,1] = t_range6[i]
						t_vec[1,0] = t_range6[i]
						new_stress = np.append(new_stress,t_vec)
						new_strain = np.append(new_strain, e_vec)
					fit_x = []
					fit_y1 = []
					fit_y2 = []
					fit_y3 = []
					fit_y4 = []	
					fit_y5 = []
					fit_y6 = []
					e_vec = []
					t_vec = []
					t_range1 = []
					t_range2 = []
					t_range3 = []
					t_range4 = []
					t_range5 = []
					t_range6 = []

			#STRAIN 2
			if tag1 == 2 and tag2 == 2:
				fit_x = np.append(fit_x, e_vec[1,1])
				fit_y2 = np.append(fit_y2, t_vec[1,1]) #C22
				fit_y3 = np.append(fit_y3, t_vec[2,2]) #C32
				fit_y4 = np.append(fit_y4, t_vec[2,1]) #C42
				fit_y5 = np.append(fit_y5, t_vec[0,2]) #C52
				fit_y6 = np.append(fit_y6, t_vec[0,1]) #C62
				if len(fit_x) == 4:
					p2 = np.polyfit(fit_x, fit_y2,2)
					p3 = np.polyfit(fit_x, fit_y3,2)
					p4 = np.polyfit(fit_x, fit_y4,2)
					p5 = np.polyfit(fit_x, fit_y5,2)
					p6 = np.polyfit(fit_x, fit_y6,2)
					t_range2 = np.append(t_range2,np.polyval(p2,e_range))
					t_range3 = np.append(t_range3,np.polyval(p3,e_range))
					t_range4 = np.append(t_range4,np.polyval(p4,e_range))
					t_range5 = np.append(t_range5,np.polyval(p5,e_range))
					t_range6 = np.append(t_range6,np.polyval(p6,e_range))
					for i in range(np.size(e_range,0)):
						e_vec[1,1] = e_range[i]
						t_vec[1,1] = t_range2[i]
						t_vec[2,2] = t_range3[i]
						t_vec[2,1] = t_range4[i]
						t_vec[1,2] = t_range4[i]
						t_vec[0,2] = t_range5[i]
						t_vec[2,0] = t_range5[i]
						t_vec[0,1] = t_range6[i]
						t_vec[1,0] = t_range6[i]
						new_stress = np.append(new_stress,t_vec)
						new_strain = np.append(new_strain, e_vec)
					fit_x = []
					fit_y2 = []
					fit_y3 = []
					fit_y4 = []
					fit_y5 = []
					fit_y6 = []
					e_vec = []
					t_vec = []
					t_range2 = []
					t_range3 = []
					t_range4 = []
					t_range5 = []
					t_range6 = []
		#STRAIN 3
			if tag1 == 4 and tag2 == 2:
				fit_x = np.append(fit_x, e_vec[2,2])
				fit_y3 = np.append(fit_y3, t_vec[2,2]) #C33
				fit_y4 = np.append(fit_y4, t_vec[2,1]) #C43
				fit_y5 = np.append(fit_y5, t_vec[0,2]) #C53
				fit_y6 = np.append(fit_y6, t_vec[0,1]) #C63
				if len(fit_x) == 4:
					p3 = np.polyfit(fit_x, fit_y3,2)
					p4 = np.polyfit(fit_x, fit_y4,2)
					p5 = np.polyfit(fit_x, fit_y5,2)
					p6 = np.polyfit(fit_x, fit_y6,2)
					t_range3 = np.append(t_range3,np.polyval(p3,e_range))
					t_range4 = np.append(t_range4,np.polyval(p4,e_range))
					t_range5 = np.append(t_range5,np.polyval(p5,e_range))
					t_range6 = np.append(t_range6,np.polyval(p6,e_range))
					for i in range(np.size(e_range,0)):
						e_vec[2,2] = e_range[i]
						t_vec[2,2] = t_range3[i]
						t_vec[2,1] = t_range4[i]
						t_vec[1,2] = t_range4[i]
						t_vec[0,2] = t_range5[i]
						t_vec[2,0] = t_range5[i]
						t_vec[0,1] = t_range6[i]
						t_vec[1,0] = t_range6[i]
						new_stress = np.append(new_stress,t_vec)
						new_strain = np.append(new_strain, e_vec)
					fit_x = []
					fit_y3 = []
					fit_y4 = []	
					fit_y5 = []
					fit_y6 = []
					e_vec = []
					t_vec = []
					t_range3 = []
					t_range4 = []
					t_range5 = []
					t_range6 = []
			#STRAIN 4
			if tag1 == 2 and tag2 == 4:
				fit_x = np.append(fit_x, e_vec[0,1])
				fit_y4 = np.append(fit_y4, t_vec[0,1]) #C66
				if len(fit_x) == 4:
					p4 = np.polyfit(fit_x, fit_y4,2)
					t_range4 = np.append(t_range4,np.polyval(p4,e_range))
					for i in range(np.size(e_range,0)):
						e_vec[0,1] = e_range[i]
						e_vec[1,0] = e_range[i]
						t_vec[0,1] = t_range4[i]
						t_vec[1,0] = t_range4[i]
						new_stress = np.append(new_stress,t_vec)
						new_strain = np.append(new_strain, e_vec)
					fit_x = []
					fit_y4 = []	
					e_vec = []
					t_vec = []
					t_range4 = []
	
			#STRAIN 5
			if tag1 == 4 and tag2 == 4:
				fit_x = np.append(fit_x, e_vec[0,2])
				fit_y4 = np.append(fit_y4, t_vec[0,1]) #C65
				fit_y5 = np.append(fit_y5, t_vec[0,2]) #C55
				if len(fit_x) == 4:
					p4 = np.polyfit(fit_x, fit_y4,2)
					p5 = np.polyfit(fit_x, fit_y5,2)
					t_range4 = np.append(t_range4,np.polyval(p4,e_range))
					t_range5 = np.append(t_range5,np.polyval(p5,e_range))
					for i in range(np.size(e_range,0)):
						e_vec[0,2] = e_range[i]
						e_vec[2,0] = e_range[i]
						t_vec[0,1] = t_range4[i]
						t_vec[1,0] = t_range4[i]
						t_vec[0,2] = t_range5[i]
						t_vec[2,0] = t_range5[i]
						new_stress = np.append(new_stress,t_vec)
						new_strain = np.append(new_strain, e_vec)
					fit_x = []
					fit_y4 = []	
					fit_y5 = []
					e_vec = []
					t_vec = []
					t_range4 = []
					t_range5 = []
	


			#STRAIN 6
			if tag1 == 6 and tag2 ==4:	
				fit_x = np.append(fit_x, e_vec[1,2])
				fit_y4 = np.append(fit_y4, t_vec[0,1])
				fit_y5 = np.append(fit_y5, t_vec[0,2])
				fit_y6 = np.append(fit_y6, t_vec[1,2])  
				if len(fit_x) == 4:
					p4 = np.polyfit(fit_x, fit_y4,2)
					p5 = np.polyfit(fit_x, fit_y5,2)
					p6 = np.polyfit(fit_x, fit_y6,2)
					t_range4 = np.append(t_range4,np.polyval(p4,e_range))
					t_range5 = np.append(t_range5,np.polyval(p5,e_range))
					t_range6 = np.append(t_range6,np.polyval(p6,e_range))
					for i in range(np.size(e_range,0)):
						e_vec[1,2] = e_range[i]
						e_vec[2,1] = e_range[i]
						t_vec[0,1] = t_range4[i]
						t_vec[1,0] = t_range4[i]
						t_vec[0,2] = t_range5[i]
						t_vec[2,0] = t_range5[i]
						t_vec[1,2] = t_range6[i]
						t_vec[2,1] = t_range6[i]
						new_stress = np.append(new_stress,t_vec)
						new_strain = np.append(new_strain, e_vec)
					fit_x = []
					fit_y4 = []
					fit_y5 = []
					fit_y6 = []
					e_vec = []
					t_vec = []
					t_range4 = []
					t_range5 = []
					t_range6 = []
			t_vec = []
			e_vec = []
	
	new_strain = np.reshape(new_strain,(-1,3))
	new_stress = np.reshape(new_stress,(-1,3))

	new_stress = (new_stress) / np.power(bohr2ang,3) * np.power(ang2m,3) * ry2joule / pa2kbar

#Zero out small strains(Should be zero!!)
	low = abs(new_strain) < 1E-10
	new_strain[low] = 0

	for i in range(np.size(new_strain,0)):
		tmp = new_stress[i,:]
		t_vec = np.append(t_vec,tmp)
		tmp2 = new_strain[i,:]
		e_vec = np.append(e_vec,tmp2)
		if np.mod(i+1,3) == 0:
		#Test for which strain is being applied
		#important : determines which C values to evaluate
		#STRAINS 1 - 6   ->> Triclinic   ; all 6 applied
		#STRAIN  7       ->> Cubic       ; 2 applied (3 and 4+5+6)

		#strain No  tag1 tag2
		# 1	 0   2
		# 2	 2   2
		# 3	 4   2
		# 4	 2   4
		# 5	 4   4
		# 6	 6   4
		# 7	 12  12

			e_vec = np.reshape(e_vec,(3,3))
			t_vec = np.reshape(t_vec,(3,3))
			tag1 = sum(np.concatenate(np.nonzero(e_vec)))
			tag2 = np.size(np.nonzero(e_vec))
		
		#REMOVE ZEROS -> CHANGE THIS.
		#Necessary because of tagging system, manually dropping a point
			if tag1 == 0 and tag2 == 0:
				e_vec = []
				t_vec = []	
	
		#STRAIN 1#
			if tag1 == 0 and tag2 == 2:
				t_vec = [t_vec[0][0], t_vec[1][1], t_vec[2][2], t_vec[2][1], t_vec[0][2], t_vec[0][1]]
				e_vec = [e_vec[0][0], e_vec[1][1], e_vec[2][2],  2 * e_vec[2][1], 2 * e_vec[0][2], 2 * e_vec[0][1]]
				t_vec = Matrix(t_vec)
				e_vec = Matrix(e_vec)
				eqs = EC_matrix.multiply(e_vec)
				eqs = eqs - t_vec 
				ans = solve([eqs[0] , eqs[1], eqs[2], eqs[3], eqs[4], eqs[5]], [C11, C12, C13, C14, C15, C16])
				C11_.append(ans[C11])
				C12_.append(ans[C12])
				C13_.append(ans[C13])
				C14_.append(ans[C14])
				C15_.append(ans[C15])
				C16_.append(ans[C16])
				t_vec = []
				e_vec = []
		
		#STRAIN 2#
			if tag1 == 2 and tag2 == 2:
				t_vec = [t_vec[0][0], t_vec[1][1], t_vec[2][2], t_vec[2][1], t_vec[0][2], t_vec[0][1]]
				e_vec = [e_vec[0][0], e_vec[1][1], e_vec[2][2], 2 * e_vec[2][1], 2 * e_vec[0][2], 2 * e_vec[0][1]]
				t_vec = Matrix(t_vec)
				e_vec = Matrix(e_vec)
				eqs = EC_matrix.multiply(e_vec)
				eqs = eqs - t_vec 
				ans = solve([eqs[0], eqs[1], eqs[2], eqs[3], eqs[4], eqs[5]], [C12, C22, C23, C24, C25, C26])
				C12_.append(ans[C12])
				C22_.append(ans[C22])
				C23_.append(ans[C23])
				C24_.append(ans[C24])
				C25_.append(ans[C25])
				C26_.append(ans[C26])
				t_vec = []
				e_vec = []

			#STRAIN 3#
			if tag1 == 4 and tag2 == 2:
				t_vec = [t_vec[0][0], t_vec[1][1], t_vec[2][2], t_vec[2][1], t_vec[0][2], t_vec[0][1]]
				e_vec = [e_vec[0][0], e_vec[1][1], e_vec[2][2], 2 * e_vec[2][1], 2 * e_vec[0][2], 2 * e_vec[0][1]]
				t_vec = Matrix(t_vec)
				e_vec = Matrix(e_vec)
				eqs = EC_matrix.multiply(e_vec)
				eqs = eqs - t_vec 
				ans = solve([eqs[0], eqs[1], eqs[2], eqs[3], eqs[4], eqs[5]], [C13, C23, C33, C34, C35, C36])
				C13_.append(ans[C13])
				C23_.append(ans[C23])
				C33_.append(ans[C33])
				C34_.append(ans[C34])
				C35_.append(ans[C35])
				C36_.append(ans[C36])
				t_vec = []
				e_vec = []
		
		#STRAIN 4#
			if tag1 == 2 and tag2 ==4:
				t_vec = [t_vec[0][0], t_vec[1][1], t_vec[2][2], t_vec[2][1], t_vec[0][2], t_vec[0][1]]
				e_vec = [e_vec[0][0], e_vec[1][1], e_vec[2][2], 2 * e_vec[2][1], 2 * e_vec[0][2], 2 * e_vec[0][1]]
				t_vec = Matrix(t_vec)
				e_vec = Matrix(e_vec)
				eqs = EC_matrix.multiply(e_vec)
				eqs = eqs - t_vec
				ans = solve([eqs[0], eqs[1], eqs[2], eqs[3], eqs[4], eqs[5]], [C16, C26, C36, C46, C56, C66])
				C16_.append(ans[C16])
				C26_.append(ans[C26])
				C36_.append(ans[C36])
				C46_.append(ans[C46])
				C56_.append(ans[C56])
				C66_.append(ans[C66])
				t_vec = []
				e_vec = []

			#STRAIN 5#
			if tag1 == 4 and tag2 ==4:
				t_vec = [t_vec[0][0], t_vec[1][1], t_vec[2][2], t_vec[2][1], t_vec[0][2], t_vec[0][1]]
				e_vec = [e_vec[0][0], e_vec[1][1], e_vec[2][2], 2 * e_vec[2][1], 2 * e_vec[0][2], 2 * e_vec[0][1]]
				t_vec = Matrix(t_vec)
				e_vec = Matrix(e_vec)
				eqs = EC_matrix.multiply(e_vec)
				eqs = eqs - t_vec
				ans = solve([eqs[0], eqs[1], eqs[2], eqs[3], eqs[4], eqs[5]], [C15, C25, C35, C45, C55, C56])
				C15_.append(ans[C15])
				C25_.append(ans[C25])
				C35_.append(ans[C35])
				C45_.append(ans[C45])
				C55_.append(ans[C55])
				C56_.append(ans[C56])
				t_vec = []
				e_vec = []
		
			#STRAIN 6#
			if tag1 == 6 and tag2 ==4:
				t_vec = [t_vec[0][0], t_vec[1][1], t_vec[2][2], t_vec[2][1], t_vec[0][2], t_vec[0][1]]
				e_vec = [e_vec[0][0], e_vec[1][1], e_vec[2][2], 2 * e_vec[2][1], 2 * e_vec[0][2], 2 * e_vec[0][1]]
				t_vec = Matrix(t_vec)
				e_vec = Matrix(e_vec)
				eqs = EC_matrix.multiply(e_vec)
				eqs = eqs - t_vec
				ans = solve([eqs[0], eqs[1], eqs[2], eqs[3], eqs[4], eqs[5]], [C14, C24, C34, C44, C45, C46])
				C14_.append(ans[C14])
				C24_.append(ans[C24])
				C34_.append(ans[C34])
				C44_.append(ans[C44])
				C45_.append(ans[C45])
				C46_.append(ans[C46])
				t_vec = []
				e_vec = []	

			#STRAIN 7#
			#CUBIC CASE - C44 EVAL#
			if tag1 == 12 and tag2 == 12:
		
			#C45
				EC_matrix[4,3] = 0
				EC_matrix[3,4] = 0

			#C46
				EC_matrix[5,3] = 0
				EC_matrix[3,5] = 0

			#C14
				EC_matrix[3,0] = 0
				EC_matrix[0,3] = 0

			#C15
				EC_matrix[0,4] = 0
				EC_matrix[4,0] = 0

			#C16
				EC_matrix[0,5] = 0 
				EC_matrix[5,0] = 0
		
			#C24
				EC_matrix[1,3] = 0
				EC_matrix[3,1] = 0

			#C25
				EC_matrix[1,4] = 0
				EC_matrix[4,1] = 0

			#C26
				EC_matrix[1,5] = 0
				EC_matrix[5,1] = 0

			#C34
				EC_matrix[2,3] = 0
				EC_matrix[3,2] = 0

			#C35
				EC_matrix[2,4] = 0
				EC_matrix[4,2] = 0

			#C56
				EC_matrix[4,5] = 0
				EC_matrix[5,4] = 0

			#C36
				EC_matrix[5,2] = 0
				EC_matrix[2,5] = 0

				t_vec = [0,0,0, t_vec[2][1], t_vec[0][2], t_vec[0][1]]
				e_vec = [e_vec[0][0], e_vec[1][1], e_vec[2][2], 2 * e_vec[2][1], 2 * e_vec[0][2], 2 * e_vec[0][1]]

				t_vec = Matrix(t_vec)
				e_vec = Matrix(e_vec)
				eqs = EC_matrix.multiply(e_vec)
				eqs = eqs - t_vec 
				ans = solve([eqs[0], eqs[1], eqs[2], eqs[3], eqs[4], eqs[5]], [C44,C55,C66])
				t_vec = []
				e_vec = []
			
	

	#Inflection Points -- EQ Elastic Constant at 0 strain
	#Strain1
	C11 = (max(C11_) + min(C11_)) / 2
	C12 = (max(C12_) + min(C12_)) / 2
	C13 = (max(C13_) + min(C13_)) / 2
	C14 = (max(C14_) + min(C14_)) / 2
	C15 = (max(C15_) + min(C15_)) / 2
	C16 = (max(C16_) + min(C16_)) / 2

	#Strain2
	C22 = (max(C22_) + min(C22_)) / 2
	C23 = (max(C23_) + min(C23_)) / 2
	C24 = (max(C24_) + min(C24_)) / 2
	C25 = (max(C25_) + min(C25_)) / 2
	C26 = (max(C26_) + min(C26_)) / 2

	#Strain3
	C33 = (max(C33_) + min(C33_)) / 2
	C34 = (max(C34_) + min(C34_)) / 2
	C35 = (max(C35_) + min(C35_)) / 2
	C36 = (max(C36_) + min(C36_)) / 2

	#Strain4
	C66 = (max(C66_) + min(C66_)) / 2

	#Strain5
	C55 = (max(C55_) + min(C55_)) / 2
	C56 = (max(C56_) + min(C56_)) / 2

	#Strain6
	C44 = (max(C44_) + min(C44_)) / 2
	C45 = (max(C45_) + min(C45_)) / 2
	C46 = (max(C46_) + min(C46_)) / 2


	EC_matrix = np.array([[C11, C12, C13, C14, C15, C16],[C12, C22, C23, C24, C25, C26], [C13, C23, C33, C34, C35, C36],[C14, C24, C34, C44, C45, C46],[C15, C25, C35, C45, C55, C56],[C16, C26, C36, C46, C56, C66]])
	EC_matrix = -(EC_matrix)

	return EC_matrix
