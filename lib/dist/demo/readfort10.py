import numpy as np
def readFort10(path):
	file = open(path, "r")
	splitlines = (file.readlines())
	print(splitlines)
	e1 = []
	e2 = []
	e3 = []
	for line in splitlines:
		print(line)
		row = (str.split(line))
		row_double = np.array(row)
		e1.append(row[46])
		e2.append(row[47])
		e3.append(row[59])
	#tasstr = (str.split(splitlines))
	#tasstr = np.array(tasstr[3:])
	#tas = tasstr.astype(np.float)
	#tasm = tas.reshape(6,6)
	#e1 = 1
	#e2 = 2
	
	return [np.array(e1), np.array(e2), np.array(e3)]


path = '/home/tobias/codes/SixTrackTobias/test/orbit6d-element-quadrupole/fort.10'
print(readFort10(path))