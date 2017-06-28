# implement a simple floor function because WDL does not have one

import fileinput
import math

for line in fileinput.input():
	value = float(line)
	intvalue = int(math.floor(value))
	print(intvalue)
