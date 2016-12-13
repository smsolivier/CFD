def writeNumber(x, fmt='{:.3}'):
	''' Returns the string in Latex format: 
			\num{fmt.format(x)} 
	''' 
	return '\\num{' + fmt.format(x) + '}' 
