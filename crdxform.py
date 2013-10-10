from numpy import matrix, dot, cos, sin, pi
def rotmat(a, b, g):
	R_a = matrix([[cos(a), -sin(a), 0],
				  [sin(a), cos(a), 0],
				  [0, 0, 1]])
	R_b = matrix([[cos(b), 0, sin(b)],
				  [0, 1, 0],
				  [-sin(b), 0, cos(b)]])
	R_g = matrix([[cos(g), -sin(g), 0],
				  [sin(g), cos(g), 0],
				  [0, 0, 1]])
	# m1 = matrix([[cos(a)*cos(b)*cos(g) - sin(a)*sin(g), sin(a)*cos(b)*cos(g) + cos(a)*sin(g), -sin(b)*cos(g)],
	# 			 [-cos(a)*cos(b)*sin(g) - sin(a)*cos(g), -sin(a)*cos(b)*sin(g) + cos(a)*cos(g), sin(b)*sin(g)],
	# 			 [cos(a)*sin(b), sin(a)*sin(b), cos(b)]])
	m = R_a*R_b*R_g
	return m



# print rotmat(pi/2.0,pi/4.0,pi/8.0)
