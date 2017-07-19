import sys
import os

key = sys.argv[1]
filename = sys.argv[2]

contamination = float('nan')
lower = float('nan')
upper = float('nan')

try:
	f = open(filename,'r')
	line = f.readline()
	fields = line.split()
	
	contamination = float(fields[0])
	lower = float(fields[1])
	upper = float(fields[2])
	f.close()
except:
	pass
finally:
	print('%s\t%s\t%f\t%s\t%f\t%s\t%f' % (key, 'contamination_schmutzi', contamination,
					   'contamination_schmutzi_lower', lower, 
					   'contamination_schmutzi_upper', upper)
	   )
