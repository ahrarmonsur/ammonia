from numpy import matrix, dot, cos, sin, pi
def rotmat(a1, b1, g1):
	m1 = matrix([[cos(a1)*cos(b1)*cos(g1) - sin(a1)*sin(g1), sin(a1)*cos(b1)*cos(g1) + cos(a1)*sin(g1), -sin(b1)*cos(g1)],
				 [-cos(a1)*cos(b1)*sin(g1) - sin(a1)*cos(g1), -sin(a1)*cos(b1)*sin(g1) + cos(a1)*cos(g1), sin(b1)*sin(g1)],
				 [cos(a1)*sin(b1), sin(a1)*sin(b1), cos(b1)]])
	vec = matrix([[1],[0],[0]])
	return m1*vec

print rotmat(-pi/2.0,0.0,0.0)
#print dir(np)