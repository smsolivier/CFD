class macro:
	''' Makes Latex macros ''' 

	nameList = [] # stores names of macros 
	valList = [] # stores values 

	def define(self, name, number, fmt='{:.2f}'):
		''' Compiles lists of latex macros ''' 

		self.nameList.append(name)
		self.valList.append(fmt.format(number))

	def save(self, outName='macros.tex'):
		''' prints list of macros and values to terminal 
			writes macros to file outName 
		''' 
		out = open(outName, 'w')

		print(outName+':')
		for i in range(len(self.nameList)):

			# print to terminal
			print('\t\\' + self.nameList[i], '=', self.valList[i])

			# write to file 
			out.write('\\newcommand{\\' + self.nameList[i] + '}{\\num{' 
				+ self.valList[i] + '}\\xspace} \n')

		out.close()
		
if __name__ == "__main__":
	m = macro() 
	m.define('test1', 5)
	m.define('test2', 10)

	m.save()
