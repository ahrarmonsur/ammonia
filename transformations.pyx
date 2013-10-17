from libc.math cimport sin, cos, M_PI
cdef double rotation[3][3]

#######################Reference Structures###########################
cdef double N_crd[3], H1_crd[3], H2_crd[3], H3_crd[3] 
N_crd[:] = [0.000000000000000, 0.000000000000000, 0.066352583]
H1_crd[:] = [0.939815876000000, 0.000000000000000, -0.307353191]
H2_crd[:] = [-0.469907938000000, -0.813904423495926, -0.307353191]
H3_crd[:] = [-0.469907938000000, 0.813904423495926, -0.307353191]

cdef void rotmat(double a, double b, double g, double rot[3][3]) except *:
	rot[0][0] = cos(a)*cos(b)*cos(g) - sin(a)*sin(g)
	rot[0][1] = -cos(a)*cos(b)*sin(g) - sin(a)*cos(g)
	rot[0][2] = cos(a)*sin(b)

	rot[1][0] = sin(a)*cos(b)*cos(g) + cos(a)*sin(g)
	rot[1][1] = -sin(a)*cos(b)*sin(g) + cos(a)*cos(g)
	rot[1][2] = sin(a)*sin(b)

	rot[2][0] = -sin(b)*cos(g)
	rot[2][1] = sin(b)*sin(g)
	rot[2][2] = cos(b)
	return

def read_profile(filename):
	inp = open(filename, 'r')
	lines = inp.readlines()
	

# cdef double * rotmat(double a, double b, double g) except *:
# 	cdef static double rot[3][3]
# 	rot[0][0] = cos(a)*cos(b)*cos(g) - sin(a)*sin(g)
# 	rot[0][1] = -cos(a)*cos(b)*sin(g) - sin(a)*cos(g)
# 	rot[0][2] = cos(a)*sin(b)

# 	rot[1][0] = sin(a)*cos(b)*cos(g) + cos(a)*sin(g)
# 	rot[1][1] = -sin(a)*cos(b)*sin(g) + cos(a)*cos(g)
# 	rot[1][2] = sin(a)*sin(b)

# 	rot[2][0] = -sin(b)*cos(g)
# 	rot[2][1] = sin(b)*sin(g)
# 	rot[2][2] = cos(b)
# 	return rot


rotmat(M_PI/2.0, 0.0, 0.0, rotation)

