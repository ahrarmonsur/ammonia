from libc.math cimport sin, cos, M_PI, sqrt
import numpy as np
#######################Reference Structures###########################
#cdef double N_crd[3], H1_crd[3], H2_crd[3], H3_crd[3] 
# N_crd[:] = [0.000000000000000, 0.000000000000000, 0.066352583]
# H1_crd[:] = [0.939815876000000, 0.000000000000000, -0.307353191]
# H2_crd[:] = [-0.469907938000000, -0.813904423495926, -0.307353191]
# H3_crd[:] = [-0.469907938000000, 0.813904423495926, -0.307353191]

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
		double t3, resid_N, resid_H
		double rcm1, rcm2, rcm3, rcm1_sq, rcm2_sq, rcm3_sq
		double rcm_sq, rcm, rcm12_sq, rcm12

		double cosb_N1, cosb_N2, cosb_N3
		double cosb_H11, cosb_H12, cosb_H13
		double cosb_H21, cosb_H22, cosb_H23
		double cosb_H31, cosb_H32, cosb_H33
		
		double a_N1, a_N2, a_N3
		double a_H11,a_H12, a_H13
		double a_H21,a_H22, a_H23
		double a_H31,a_H32, a_H33		

	t3 = 1 / (m_N + 3 * m_H);
	resid_N = -t3 * m_N + 1;
	resid_H = t3 * m_H
	rcm1 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
	rcm2 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
	rcm3 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
	rcm1_sq = rcm1 * rcm1;
	rcm2_sq = rcm2 * rcm2;
	rcm3_sq = rcm3 * rcm3;
	rcm_sq = rcm1_sq + rcm2_sq + rcm3_sq;
	rcm = sqrt(rcm_sq);

	rcm12_sq = rcm1_sq + rcm2_sq;
	rcm12 = sqrt(rcm12_sq);
	a_arccos_deriv = (sqrt(1 - rcm1_sq / rcm12_sq));




	## NOTE: Hi1, Hi2, Hi3 are each the same for all i's
	## Should we use one result 3 times each for efficiency?
	cosb_N1 = -rcm3 / rcm / rcm_sq * rcm1 * resid_N;
	cosb_N2 = -rcm3 / rcm / rcm_sq * rcm2 * resid_N;
	cosb_N3 = resid_N / rcm - rcm3_sq / rcm / rcm_sq * resid_N;
	cosb_H11 = rcm3 / rcm / rcm_sq * rcm1 * resid_H;
	cosb_H12 = rcm3 / rcm / rcm_sq * rcm2 * resid_H;
	cosb_H13 = -resid_H / rcm + rcm3_sq / rcm / rcm_sq * resid_H;
	cosb_H21 =  rcm3 / rcm / rcm_sq * rcm1 * resid_H;
	cosb_H22 = rcm3 / rcm / rcm_sq * rcm2 * resid_H;
	cosb_H23 = -resid_H / rcm  + rcm3_sq / rcm / rcm_sq * resid_H;
	cosb_H31 = rcm3 / rcm / rcm_sq * rcm1 * resid_H;
	cosb_H32 = rcm3 / rcm / rcm_sq * rcm2 * resid_H;
	cosb_H33 = -resid_H / rcm + rcm3_sq / rcm / rcm_sq * resid_H;

	a_N1 = -(resid_N / rcm12 - rcm1_sq / rcm12 / rcm12_sq * resid_N) / a_arccos_deriv;
	a_N2 = rcm1 / rcm12 / rcm12_sq * rcm2 * resid_N / a_arccos_deriv;
	a_N3 = 0.0


	print a_N2, a_N1

#derivatives(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 12.0, 2.0)
derivatives(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 20.0, 50.0)
###########Test for rotmat###########
#cdef double rotation[3][3]

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



