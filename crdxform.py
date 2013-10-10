#######################Reference Structures###########################
import numpy as np
N_crd = np.matrix([[0.000000000000000],[0.000000000000000],[0.066352583]])
H1_crd = np.matrix([[0.939815876000000],[0.000000000000000],[-0.307353191]])
H2_crd = np.matrix([[-0.469907938000000],[-0.813904423495926],[-0.307353191]])
H3_crd = np.matrix([[-0.469907938000000],[0.813904423495926],[-0.307353191]])

def data_init(filename):
	test = open(filename, "r")
	lines = [item.split() for item in test.readlines()]
	data = []
	for coords in lines:
		data.append([float(val) for val in coords])
	return data

def rotmat(a, b, g):
	import numpy as np
	R_a = np.matrix([[np.cos(a), -np.sin(a), 0],
				  [np.sin(a), np.cos(a), 0],
				  [0, 0, 1]])
	b = np.arccos(b)
	R_b = np.matrix([[np.cos(b), 0, np.sin(b)],
				  [0, 1, 0],
				  [-np.sin(b), 0, np.cos(b)]])
	R_g = np.matrix([[np.cos(g), -np.sin(g), 0],
				  [np.sin(g), np.cos(g), 0],
				  [0, 0, 1]])
	# m1 = np.matrix([[np.cos(a)*np.cos(b)*np.cos(g) - np.sin(a)*np.sin(g), np.sin(a)*np.cos(b)*np.cos(g) + np.cos(a)*np.sin(g), -np.sin(b)*np.cos(g)],
	# 			 [-np.cos(a)*np.cos(b)*np.sin(g) - np.sin(a)*np.cos(g), -np.sin(a)*np.cos(b)*np.sin(g) + np.cos(a)*np.cos(g), np.sin(b)*np.sin(g)],
	# 			 [np.cos(a)*np.sin(b), np.sin(a)*np.sin(b), np.cos(b)]])
	m = R_a*R_b*R_g
	return m

def transform(incrd):
	import numpy as np

	R, g1, b1, g2, b2, a, pot = incrd
	rotM1 = rotmat(a/2.0, b1, g1)
	rotM2 = rotmat(-a/2.0, b2, g2)
	cm1 = np.matrix([[0.0], [0.0], [0.0]])
	cm2 = np.matrix([[0.0], [0.0], [R]])

	N1 = np.array(rotM1*N_crd + cm1).flatten().tolist()
	H11 = np.array(rotM1*H1_crd + cm1).flatten().tolist()
	H12 = np.array(rotM1*H2_crd + cm1).flatten().tolist()
	H13 = np.array(rotM1*H3_crd + cm1).flatten().tolist()

	N2 = np.array(rotM1*N_crd + cm2).flatten().tolist()
	H21 = np.array(rotM1*H1_crd + cm2).flatten().tolist()
	H22 = np.array(rotM1*H2_crd + cm2).flatten().tolist()
	H23 = np.array(rotM1*H3_crd + cm2).flatten().tolist()

	return [N1, H11, H12, H13, N2, H21, H22, H23, pot]

def int_cart(filename):
	import numpy as np
	data = [transform(coord) for coord in data_init(filename)]
	return data


# print transform([2.5, 0.0, -0.986283808696811, 0.0, -0.986283808696811, 0.0, -113.09755364040774])


# data = data_init("testdata.txt")
# for i in range(10):
#  	print data[i]
print int_cart("testdata.txt")

# print rotmat(np.pi/2.0,np.pi/4.0,np.pi/8.0)
