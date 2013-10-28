import numpy as np


# def index(R_ind, g1_ind, b1_ind, g2_ind, b2_ind, a_ind):
# 	r = 41
# 	g1 = 10
# 	b1 = 14
# 	g2 = 10
# 	b2 = 14
# 	a = 16
# 	ind = R_ind*(g1*b1*g2*b2*a) + g1_ind*(b1*g2*b2*a) + b1_ind*(g2*b2*a) + g2_ind*(b2*a) + b2_ind*(a) + a_ind
# 	return ind
# print index(0, 8, 8, 6, 1, 10)

############################# FINDING CHARACTERISTICS OF INTERNAL COORDINATES #############################
# b1 = []; g1 = []; b2 = []; g2 = []; a = []; b_angle =[]
# inp = open("../dawesr/fort.941", "r")
# line = inp.readline()
# while line != "":
# 	data = line.split()
# 	if data[1] not in g1: g1.append(data[1])
# 	if data[2] not in b1: b1.append(data[2])
# 	if data[3] not in g2: g2.append(data[3])
# 	if data[4] not in b2: b2.append(data[4])
# 	if data[5] not in a: a.append(data[5])
# 	line = inp.readline()
# print len(a)
# print len(g1)
# print len(b1)
# print len(g2)
# print len(b2)
# print 16*14*10*14*10
#for i in range(len(a)-1): print float(a[i+1]) - float(a[i])
# for item in b1: b_angle.append(np.arc	)


############################# EXTRACTING A FLATTENED POTENTIAL VECTOR #############################
# outp = open("nh3_pots.txt", "a")
# for i in range(41, 82):
# 	inp = open("../dawesr/fort.9"+str(i), "r")
# 	line = inp.readline()
# 	while line != "":
# 		data = line.split()[6]
# 		outp.write(data+"\n")
# 		line = inp.readline()
# 	inp.close()
# outp.close()


# def read_profile(filename):
# 	inp = open(filename, 'r')
# 	lines = inp.readlines()
# 	cur_body = None
# 	labels = []
# 	masses = []
# 	coords = [[],[]]
# 	for item in lines:
# 		line = item.split()
# 		if item[0] == "#": pass
# 		elif line == []: pass
# 		elif line == ["Body1"]:
# 			cur_body = 0
# 		elif line == ["Body2"]:
# 			cur_body = 1
# 		elif len(line) == 5:
# 			labels.append(line[0])
# 			masses.append(line[1])
# 			coords[cur_body].extend(line[2:])
# 		else: print "Parsing error. Check format of the profile: %s" %(filename)
# 	return labels, masses, coords


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
a =  read_profile("testprofile.txt")

import matplotlib.pyplot as plt
help(plt.legend)