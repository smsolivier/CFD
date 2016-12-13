import numpy as np 
import texTools as tex 

case = 'N339' 

fields = ['U', 'k', 'C'] 

table = tex.table()

for i in range(len(fields)):

	dr = 'InputUncertainty/' + case + '-' + fields[i] + '.csv'

	f = open(dr, 'r')

	coefs = [] # store list of coeffs 
	u = [] 

	next(f) # skip header 
	for line in f:
		if line.startswith('Overall'):
			break 

		val = line.split(',')

		coefs.append(float(val[1]))
		u.append(float(val[2]))

	f.close()

	table.addLine(
		fields[i],
		tex.utils.writeNumber(coefs[0]),
		tex.utils.writeNumber(coefs[1]),
		tex.utils.writeNumber(coefs[2]),
		tex.utils.writeNumber(coefs[3]),
		tex.utils.writeNumber(coefs[4]),
		tex.utils.writeNumber(coefs[5]),
		tex.utils.writeNumber(coefs[6])
		)

table.save('report/coefs.tex')