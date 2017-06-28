import sys
import os

key = sys.argv[1]
filename = sys.argv[2]

contamination = float('nan')
lower = float('nan')
upper = float('nan')

try:
	f = open(filename,'r')
	line = f.readline() # ignore header line
	line = f.readline()
	fields = line.split()
	
	contamination = float(fields[3])
	lower = float(fields[4])
	upper = float(fields[5])
	f.close()
except:
	pass
finally:
	print('%s\t%s\t%f\t%s\t%f\t%s\t%f' % (key, 'contamination_rare_variant', contamination,
					   'contamination_rare_variant_lower', lower, 
					   'contamination_rare_variant_upper', upper)
	   )

