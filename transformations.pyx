from libc.math cimport sin, cos, M_PI, sqrt, acos
import numpy as np
#######################Reference Structures###########################
#cdef double N_crd[3], H1_crd[3], H2_crd[3], H3_crd[3] 
# N_crd[:] = [0.000000000000000, 0.000000000000000, 0.066352583]
# H1_crd[:] = [0.939815876000000, 0.000000000000000, -0.307353191]
# H2_crd[:] = [-0.469907938000000, -0.813904423495926, -0.307353191]
# H3_crd[:] = [-0.469907938000000, 0.813904423495926, -0.307353191]

####################### WRITING TO BINARY FILE ###########################
# from array import array

## Reads ascii file and converts to list of floats
# inpf = open('nh3_pots.txt', 'r')
# data = inpf.readlines()
# data2 = map(lambda x: float(x.strip()), data)

## converts to constrained array of floats and writes to bin file
# output_file = open('nh3_pot_bin', 'wb')
# float_array = array('d', data2)
# float_array.tofile(output_file)
# output_file.close()

## reads bin file and converts to constrained float array
# input_file = open('nh3_pot_bin', 'rb')
# float_array = array('d')
# float_array.fromstring(input_file.read())
# input_file.close()
# print len(float_array)

####################### EXTRACTING A FLATTENED POTENTIAL VECTOR #######################
# outp = open("nh3_pots.txt", "a")
# for i in range(41, 82):
# 	inp = open("../dawesr/fort.9"+str(i), "r")
# 	line = inp.readline()
# 	while line != "":
# 		data = line.split()[6]
# 		outp.write(data + "\n")
# 		line = inp.readline()
# 	inp.close()
# outp.close()

####################### CARTESIAN COORDINATES DATA STRUCTURE #######################
# Memory view into array of form
# [N1, N2, N3, H11, H12, H13, H21, H22, H23, H31, H32, H33, n1, n2, n3, h11, h12, h13, h21, h22, h23, h31, h32, h33]


def index(R_ind, g1_ind, b1_ind, g2_ind, b2_ind, a_ind):
	r = 41
	g1 = 10
	b1 = 14
	g2 = 10
	b2 = 14
	a = 16
	ind = R_ind*(g1*b1*g2*b2*a) + g1_ind*(b1*g2*b2*a) + b1_ind*(g2*b2*a) + g2_ind*(b2*a) + b2_ind*(a) + a_ind
	return ind
# print index(0, 8, 8, 6, 1, 10)

# Takes a profile filename and returns a list of 2 body data containing
# labels, masses, and flattened coordinates lists
def read_profile(filename):
	import numpy as np
	inp = open(filename, 'r')
	lines = inp.readlines()
	inp.close()
	cur_body = None
	bodies_data = [[[],[],[]],[[],[],[]]]
	for item in lines:
		line = item.split()
		if item[0] == "#": pass
		elif line == []: pass
		elif line == ["Body1"]:
			cur_body = 0
		elif line == ["Body2"]:
			cur_body = 1
		elif len(line) == 5:
			bodies_data[cur_body][0].append(line[0])
			bodies_data[cur_body][1].append(line[1])
			bodies_data[cur_body][2].extend(line[2:])
		else: print "Parsing error. Check format of the profile: %s" %(filename)
	return bodies_data

# Input: 1/ double: alpha, 2/ double: beta, 3/ double: gamma`
# Output: 4/ double 2D memvw (len 3x3)
cdef void rotmat(double a, double b, double g, double[:,:] rot_op) except *:
	rot_op[0][0] = cos(a)*cos(b)*cos(g) - sin(a)*sin(g)
	rot_op[0][1] = -cos(a)*cos(b)*sin(g) - sin(a)*cos(g)
	rot_op[0][2] = cos(a)*sin(b)

	rot_op[1][0] = sin(a)*cos(b)*cos(g) + cos(a)*sin(g)
	rot_op[1][1] = -sin(a)*cos(b)*sin(g) + cos(a)*cos(g)
	rot_op[1][2] = sin(a)*sin(b)

	rot_op[2][0] = -sin(b)*cos(g)
	rot_op[2][1] = sin(b)*sin(g)
	rot_op[2][2] = cos(b)
	return

# Input: 1/ double memvw (len n), 2/ double memvw (len 3n)
# Output: 3/ double memvw (len 3)
cdef void centre_of_mass(double[:] masses, double[:] coords, double[:] com_op) except *:
	cdef int i, natom
	cdef double total_mass = 0.0
	cdef double [:] masses_vw = masses
	natom = masses_vw.shape[0]
	if natom < 1: print "Error: Need atleast 1 atom for centre of mass calculation"
	# Perform a looping addition before dividing my total mass of body
	for i in range(natom):
		total_mass += masses[i]
		com_op[0] += masses[i]*coords[i*3 + 0]
		com_op[1] += masses[i]*coords[i*3 + 1]
		com_op[2] += masses[i]*coords[i*3 + 2]
	com_op[0] /= total_mass
	com_op[1] /= total_mass
	com_op[2] /= total_mass 
	return

# Takes Euler angles and the reference structure of atom (in one body) and the centre of mass
# and returns the cartesian coordinates of the atom
# Input: 1/ double, 2/ double, 3/ double memvw (len 3), 4/ double memvw (len 3)
# Output: 5/ double memvw (len 3)
cdef void int_cart(double a, double b, double g, double[:] ref_struct, double[:] com, double[:] cart_op) except *:
	cdef double rot_mat[3][3]
	cdef double [:,:] rotmat_vw = rot_mat
	cdef int i, j
	cdef double tot
	rotmat(a, b, g, rotmat_vw)
	# Perform matrix-vector multiplication; apply rotation matrix to ref position
	for i in range(3):
		tot = 0.0
		for j in range(3):
			tot += rotmat_vw[i][j] * ref_struct[j]
		cart_op[i] = tot + com[i]
	return

cdef void derivatives(double N1, N2, N3, H11, H12, H13, H21, H22, H23,\
  H31, H32, H33, m_N, m_H) except *:
	cdef:
		double t3, cosb_N_deriv, cosb_H_deriv, resid_N, resid_H
		double Ncm1, Ncm2, Ncm3, Ncm1_sq, Ncm2_sq, Ncm3_sq
		double H1cm1, H1cm2, H1cm3
		double Ncm_sq, Ncm, Ncm12_sq, Ncm12, a_arccos_deriv
		double g_c1, g_c2, g_c345_sq, g_c345
		double g_c3, g_c4, g_c5, g_c3_sq, g_c4_sq, g_c5_sq
		double g_c6, g_c7, g_c8, g_c9, g_c10, g_c11, g_c12, g_c13, g_c14

		double cosb_N1, cosb_N2, cosb_N3
		double cosb_H11, cosb_H12, cosb_H13
		double cosb_H21, cosb_H22, cosb_H23
		double cosb_H31, cosb_H32, cosb_H33
		
		double a_N1, a_N2, a_N3
		double a_H11,a_H12, a_H13
		double a_H21,a_H22, a_H23
		double a_H31,a_H32, a_H33		

	t3 = 1 / (m_N + 3 * m_H)
	cosb_N_deriv = 1.0 - t3 * m_N 
	cosb_H_deriv = 1.0 - t3 * m_H
	resid_N = t3 * m_N
	resid_H = t3 * m_H
	Ncm1 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))
	Ncm2 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))
	Ncm3 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))

	H1cm1 = H11 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))
	H1cm2 = H12 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))
	H1cm3 = H13 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))

	Ncm1_sq = Ncm1 * Ncm1
	Ncm2_sq = Ncm2 * Ncm2
	Ncm3_sq = Ncm3 * Ncm3
	Ncm_sq = Ncm1_sq + Ncm2_sq + Ncm3_sq
	Ncm = sqrt(Ncm_sq)

	Ncm12_sq = Ncm1_sq + Ncm2_sq
	Ncm12 = sqrt(Ncm12_sq)
	a_arccos_deriv = (sqrt(1 - Ncm1_sq / Ncm12_sq))

	g_c1 = H1cm2 * resid_H - Ncm2 * resid_H
	g_c2 = Ncm1 * resid_H - H1cm1 * resid_H
	g_c3 = Ncm2 * H1cm3 - Ncm3 * H1cm2
	g_c4 = Ncm3 * H1cm1 - Ncm1 * H1cm3
	g_c5 = Ncm1 * H1cm2 - Ncm2 * H1cm1
	g_c6 = g_c4 * Ncm1 / Ncm12 - g_c3 * a_arccos_deriv
	g_c3_sq = g_c3 * g_c3
	g_c4_sq = g_c4 * g_c4
	g_c5_sq = g_c5 * g_c5
	g_c6_sq = g_c6 * g_c6
	g_c345_sq = g_c3_sq + g_c4_sq + g_c5_sq
	g_c345 = sqrt(g_c345_sq)
	g_c7 = sqrt(1.0 - g_c6_sq / g_c345_sq)
	g_c8 = Ncm1 * resid_N + H1cm1 * cosb_N_deriv
	g_c9 = Ncm2 * resid_N + H1cm2 * cosb_N_deriv
	g_c10 = Ncm3 * resid_N + H1cm3 * cosb_N_deriv
	g_c11 = H1cm1 * resid_H + Ncm1 * cosb_H_deriv
	g_c12 = H1cm2 * resid_H + Ncm2 * cosb_H_deriv
	g_c13 = H1cm3 * resid_H + Ncm3 * cosb_H_deriv
	g_c14 = H1cm3 * resid_H - Ncm3 * resid_H


	## NOTE: Hi1, Hi2, Hi3 are each the same for all i's
	##       except for gamma  
	## Should we use one result 3 times each for efficiency?
	cosb_N1 = -Ncm3 / Ncm / Ncm_sq * Ncm1 * cosb_N_deriv
	cosb_N2 = -Ncm3 / Ncm / Ncm_sq * Ncm2 * cosb_N_deriv
	cosb_N3 = cosb_N_deriv / Ncm - Ncm3_sq / Ncm / Ncm_sq * cosb_N_deriv
	cosb_H11 = Ncm3 / Ncm / Ncm_sq * Ncm1 * resid_H
	cosb_H12 = Ncm3 / Ncm / Ncm_sq * Ncm2 * resid_H
	cosb_H13 = -resid_H / Ncm + Ncm3_sq / Ncm / Ncm_sq * resid_H
	cosb_H21 =  Ncm3 / Ncm / Ncm_sq * Ncm1 * resid_H
	cosb_H22 = Ncm3 / Ncm / Ncm_sq * Ncm2 * resid_H
	cosb_H23 = -resid_H / Ncm  + Ncm3_sq / Ncm / Ncm_sq * resid_H
	cosb_H31 = Ncm3 / Ncm / Ncm_sq * Ncm1 * resid_H
	cosb_H32 = Ncm3 / Ncm / Ncm_sq * Ncm2 * resid_H
	cosb_H33 = -resid_H / Ncm + Ncm3_sq / Ncm / Ncm_sq * resid_H

	a_N1 = -(cosb_N_deriv / Ncm12 - Ncm1_sq / Ncm12 / Ncm12_sq * cosb_N_deriv) / a_arccos_deriv
	a_N2 = Ncm1 / Ncm12 / Ncm12_sq * Ncm2 * cosb_N_deriv / a_arccos_deriv
	a_N3 = 0.0
	a_H11 = (resid_H / Ncm12 - Ncm1_sq / Ncm12 / Ncm12_sq * resid_H) / a_arccos_deriv
	a_H12 = -Ncm1 / Ncm12 / Ncm12_sq * Ncm2 * resid_H / a_arccos_deriv
	a_H13 = 0.0
	a_H21 = (resid_H / Ncm12 - Ncm1_sq / Ncm12 / Ncm12_sq * resid_H) / a_arccos_deriv
	a_H22 = -Ncm1 / Ncm12 / Ncm12_sq * Ncm2 * resid_H / a_arccos_deriv
	a_H23 = 0.0
	a_H31 = (resid_H / Ncm12 - Ncm1_sq / Ncm12 / Ncm12_sq * resid_H) / a_arccos_deriv
	a_H32 = -Ncm1 / Ncm12 / Ncm12_sq * Ncm2 * resid_H / a_arccos_deriv
	a_H33 = 0.0

	g_N1 = -((g_c3 / a_arccos_deriv * (Ncm1_sq * Ncm1  * cosb_N_deriv / (Ncm12_sq * Ncm12_sq) - Ncm1 * cosb_N_deriv / Ncm12_sq)\
		- g_c10 * Ncm1 / Ncm12 + g_c4 * cosb_N_deriv / Ncm12  - g_c4 * Ncm1_sq * cosb_N_deriv / Ncm12 / Ncm12_sq) / g_c345\
		- g_c6 * (-g_c4 * g_c10 + g_c5 * g_c9) / g_c345 / g_c345_sq) / g_c7
	
	g_N2 = -((-g_c10 * a_arccos_deriv\
			- g_c3 * Ncm1_sq * Ncm2 * cosb_N_deriv / a_arccos_deriv / (Ncm12_sq * Ncm12_sq)\
			- g_c4 * Ncm1 * Ncm2 * cosb_N_deriv / Ncm12 / Ncm12_sq) / g_c345\
		- g_c6 * (g_c3 * g_c10 - g_c5 * g_c8) / g_c345 / g_c345_sq) / g_c7

	g_N3 = -((g_c8 * Ncm1 / Ncm12 + g_c9 * a_arccos_deriv) / g_c345 - g_c6 * (g_c4 * g_c8 - g_c3 * g_c9) / g_c345 / g_c345_sq) / g_c7

	g_H11 = -((-g_c3 * (Ncm1 * resid_H / Ncm12_sq - Ncm1_sq * Ncm1 * resid_H / (Ncm12_sq * Ncm12_sq)) / a_arccos_deriv\
		+ g_c13 * Ncm1 / Ncm12\
		- g_c4 * resid_H / Ncm12\
		+ g_c4 * Ncm1_sq * resid_H / Ncm12 / Ncm12_sq) / g_c345\
	- g_c6 * (g_c4 * g_c13 - g_c5 * g_c12) / g_c345 / g_c345_sq) / g_c7
 
	g_H12 = -((g_c13 * a_arccos_deriv\
		+ g_c3 * Ncm1_sq * Ncm2 * resid_H / a_arccos_deriv / (Ncm12_sq * Ncm12_sq)\
		+ g_c4 * Ncm1 * Ncm2 * resid_H / Ncm12_sq / Ncm12) / g_c345\
	- g_c6 * (g_c5 * g_c11 - g_c3 * g_c13) / g_c345 / g_c345_sq) / g_c7

	g_H13 = -((-g_c11 * Ncm1 / Ncm12 - g_c12 * a_arccos_deriv) / g_c345\
		- g_c6 * (g_c3 * g_c12 - g_c4 * g_c11) / g_c345 / g_c345_sq) / g_c7

	g_H21 = -((-g_c3 * (Ncm1 * resid_H / Ncm12_sq - Ncm1_sq * Ncm1 * resid_H / (Ncm12_sq * Ncm12_sq)) / a_arccos_deriv\
		+ g_c14 * Ncm1 / Ncm12\
		- g_c4 * resid_H / Ncm12\
		+ g_c4 * Ncm1_sq * resid_H / Ncm12 / Ncm12_sq) / g_c345\
	- g_c6 * (g_c4 * g_c14 - g_c5 * g_c1) / g_c345 / g_c345_sq) / g_c7

	g_H22 = -((g_c14 * a_arccos_deriv\
		+ g_c3 * Ncm1_sq * Ncm2 * resid_H / a_arccos_deriv / (Ncm12_sq * Ncm12_sq)\
		+ g_c4 * Ncm1 * Ncm2 * resid_H / Ncm12_sq / Ncm12) / g_c345\
	- g_c6 * (-g_c5 * g_c2 - g_c3 * g_c14) / g_c345 / g_c345_sq) / g_c7

	g_H23 = -((g_c2 * Ncm1 / Ncm12 - g_c1 * a_arccos_deriv) / g_c345\
		- g_c6 / g_c345 / g_c345_sq * (g_c3 * g_c1 + g_c4 * g_c2)) / g_c7\

	g_H31 = -((-g_c3 * (Ncm1 * resid_H / Ncm12_sq - Ncm1_sq * Ncm1 * resid_H / (Ncm12_sq * Ncm12_sq)) / a_arccos_deriv\
		+ g_c14 * Ncm1 / Ncm12\
		- g_c4 * resid_H / Ncm12\
		+ g_c4 * Ncm1_sq * resid_H / Ncm12 / Ncm12_sq) / g_c345\
	- g_c6 * (g_c4 * g_c14 - g_c5 * g_c1) / g_c345 / g_c345_sq) / g_c7


	g_H32 = -((g_c14 * a_arccos_deriv\
		+ g_c3 * Ncm1_sq * Ncm2 * resid_H / a_arccos_deriv / (Ncm12_sq * Ncm12_sq)\
		+ g_c4 * Ncm1 * Ncm2 * resid_H / Ncm12_sq / Ncm12) / g_c345\
	- g_c6 * (-g_c5 * g_c2 - g_c3 * g_c14) / g_c345 / g_c345_sq) / g_c7

	g_H33 = -((g_c2 * Ncm1 / Ncm12 - g_c1 * a_arccos_deriv) / g_c345\
		- g_c6 / g_c345 / g_c345_sq * (g_c3 * g_c1 + g_c4 * g_c2)) / g_c7

	print g_H32
	return


#derivatives(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 12.0, 2.0)
derivatives(10.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 20.0, 50.0)
###########Test for rotmat###########
#cdef double rotation[3][3]

# [N1, N2, N3, H11, H12, H13, H21, H22, H23, H31, H32, H33, n1, n2, n3, h11, h12, h13, h21, h22, h23, h31, h32, h33]
# cdef void cart_int(double[:] carts, double m_N, m_H, double [:] ints) except *:
# 	t3 = 1 / (m_N + 3 * m_H)
# 	cosb_N_deriv = 1.0 - t3 * m_N 
# 	cosb_H_deriv = 1.0 - t3 * m_H
# 	resid_N = t3 * m_N
# 	resid_H = t3 * m_H
# 	Ncm1 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))
# 	Ncm2 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))
# 	Ncm3 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))

# 	Ncm1_sq = Ncm1 * Ncm1
# 	Ncm2_sq = Ncm2 * Ncm2
# 	Ncm3_sq = Ncm3 * Ncm3

# 	Ncm_sq = Ncm1_sq + Ncm2_sq + Ncm3_sq
# 	Ncm = sqrt(Ncm_sq)
# 	Ncm12_sq = Ncm1_sq + Ncm2_sq
# 	Ncm12 = sqrt(Ncm12_sq)

# 	H1cm1 = H11 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))
# 	H1cm2 = H12 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))
# 	H1cm3 = H13 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))

# 	g_c3 = Ncm2 * H1cm3 - Ncm3 * H1cm2
# 	g_c4 = Ncm3 * H1cm1 - Ncm1 * H1cm3
# 	g_c5 = Ncm1 * H1cm2 - Ncm2 * H1cm1
# 	g_c3_sq = g_c3 * g_c3
# 	g_c4_sq = g_c4 * g_c4
# 	g_c5_sq = g_c5 * g_c5
# 	g_c345_sq = g_c3_sq + g_c4_sq + g_c5_sq
# 	g_c345 = sqrt(g_c345_sq)
# 	a_arccos_deriv = (sqrt(1 - Ncm1_sq / Ncm12_sq))

# 	ints[0] = acos(Ncm1 / Ncm12)
# 	ints[1] = Ncm3 / Ncm
# 	ints[2] = acos((-g_c3 * a_arccos_deriv + g_c4 * Ncm1 / Ncm12) / g_c345)

cdef void cart_int(double[:] carts, double [:] ints, double m_N, m_H) except *:
	t3 = 1 / (m_N + 3 * m_H)
	cosb_N_deriv = 1.0 - t3 * m_N 
	cosb_H_deriv = 1.0 - t3 * m_H
	resid_N = t3 * m_N
	resid_H = t3 * m_H
	Ncm1 = carts[0] - t3 * (m_N * carts[0] + m_H * (carts[3] + carts[6] + carts[9]))
	Ncm2 = carts[1] - t3 * (m_N * carts[1] + m_H * (carts[4] + carts[7] + carts[10]))
	Ncm3 = carts[2] - t3 * (m_N * carts[2] + m_H * (carts[5] + carts[8] + carts[11]))

	Ncm1_sq = Ncm1 * Ncm1
	Ncm2_sq = Ncm2 * Ncm2
	Ncm3_sq = Ncm3 * Ncm3

	Ncm_sq = Ncm1_sq + Ncm2_sq + Ncm3_sq
	Ncm = sqrt(Ncm_sq)
	Ncm12_sq = Ncm1_sq + Ncm2_sq
	Ncm12 = sqrt(Ncm12_sq)

	H1cm1 = carts[3] - t3 * (m_N * carts[0] + m_H * (carts[3] + carts[6] + carts[9]))
	H1cm2 = carts[4] - t3 * (m_N * carts[1] + m_H * (carts[4] + carts[7] + carts[10]))
	H1cm3 = carts[5] - t3 * (m_N * carts[2] + m_H * (carts[5] + carts[8] + carts[11]))

	g_c3 = Ncm2 * H1cm3 - Ncm3 * H1cm2
	g_c4 = Ncm3 * H1cm1 - Ncm1 * H1cm3
	g_c5 = Ncm1 * H1cm2 - Ncm2 * H1cm1
	g_c3_sq = g_c3 * g_c3
	g_c4_sq = g_c4 * g_c4
	g_c5_sq = g_c5 * g_c5
	g_c345_sq = g_c3_sq + g_c4_sq + g_c5_sq
	g_c345 = sqrt(g_c345_sq)
	a_arccos_deriv = (sqrt(1 - Ncm1_sq / Ncm12_sq))

	ints[0] = acos(Ncm1 / Ncm12)
	ints[1] = Ncm3 / Ncm
	ints[2] = acos((-g_c3 * a_arccos_deriv + g_c4 * Ncm1 / Ncm12) / g_c345)
	return

cdef double interpolator(double[:] carts, double m_N, m_H) except *:
	cdef double ints[3]
	cdef double [:] ints_vw = ints

	## Need to use centre_of_mass to find R (z value of com)
	## subtract R from z direction values (n*3+2 index)
	## recentre monomer 2 and find internal angles
	## interpolate from data characteristics

	a = carts[:3] + carts[4:7]
	for i in a: print i
	# cart_int(carts, ints, m_N, m_H)



testints = np.zeros(3)
cdef double [:] testcarts_vw = np.array([10.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
cdef double [:] testints_vw = testints
cart_int(testcarts_vw, testints_vw, 20.0, 50.0)
interpolator(testcarts_vw, 20.0, 50.0)
print testints





###########Test for centre_of_mass###########
# cdef enum: t_atom2 = 4
# cdef double mat[t_atom2]
# cdef:
# 	double t_masses[4]
# 	double [:] t_masses_vw = t_masses
# 	double t_coords[12]
# 	double [:] t_coords_vw = t_coords
# 	double t_com[3]
# 	double [:] t_com_vw = t_com
# t_masses[:] = [14.0, 1.0, 1.0, 1.0]
# t_coords[:] = [0.000000000000000, 0.000000000000000, 0.066352583,
# 0.939815876000000, 0.000000000000000, -0.307353191,
# -0.469907938000000, -0.813904423495926, -0.307353191,
# -0.469907938000000, 0.813904423495926, -0.307353191]
# centre_of_mass(t_masses_vw, t_coords_vw, t_com_vw)
# print t_com_vw[0], t_com_vw[1], t_com_vw[2]


###########Test for int_cart###########
# cdef double [:] t_ref = np.array([0.000000000000000, 0.000000000000000, 0.066352583])
# cdef double [:] t_com = np.array([0.0, 0.0, 10.0])
# cdef double t_cart[3]
# cdef double [:] t_cart_vw = t_cart
# int_cart(0.0, M_PI/2.0, 0.0, t_ref, t_com, t_cart_vw)
# print [t_cart_vw[0], t_cart_vw[1], t_cart_vw[2]]



