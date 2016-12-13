class table:
	''' Makes a Latex formatted table ''' 

	str = '' # stores table 

	def addLine(self, *args):
		''' stores list of values into str with & between and \\ \n at the end ''' 
		
		n = len(args)

		for i in range(n - 1):
			self.str += args[i] + ' & ' # add args to str

		self.str += args[-1] + ' \\\\\n'

	def save(self, outName='macros.tex'):
		''' Saves str to file outName ''' 

		out = open(outName, 'w')

		out.write(self.str)

		out.close()

		print('created table ', outName)

		self.str = '' # reset str 

if __name__ == '__main__':
	t = table()

	t.addLine('lineName', 'hello', 'world')
	t.addLine('test2', 'test2')

	t.save()