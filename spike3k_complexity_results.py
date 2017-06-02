import sys
import os

# read in 
filename = sys.argv[1]
with open(filename, "r") as f:
	for line in f:
		if line[0] != '#':
			fields = line.strip().split()
			key = fields[0]
			complexity_estimate = fields[3]
			print('%s\tspike3k_complexity\t%s' % (key, complexity_estimate) )
