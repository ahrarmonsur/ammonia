# index = [[["N1", "N2", "N3"],
# 		  ["H11", "H12", "H13"],
# 		  ["H21", "H22", "H23"],
# 		  ["H31", "H32", "H33"]],
# 		 [["n1", "n2", "n3"],
# 		  ["h11", "h12", "h13"],
# 		  ["h21", "h22", "h23"],
# 		  ["h31", "h32", "h33"]]]
# print index[1][2][2]
# def fun(x):
# 	return 7*x

# dic = [fun,3,2]
# print dic[0](3)
from array import array

# inpf = open('nh3_pots.txt', 'r')
# data = inpf.readlines()
# data2 = map(lambda x: float(x.strip()), data)

# output_file = open('nh3_pot_bin', 'wb')
# float_array = array('d', data2)
# float_array.tofile(output_file)
# output_file.close()


input_file = open('nh3_pot_bin', 'rb')
float_array = array('d')
float_array.fromstring(input_file.read())
input_file.close()
print len(float_array)