import os
import sys
# read in pmdtools input --first damage statistics from stdin
# 1-2-Q3-Q4.bam
# C>T at first position and SE: 0.1   0.001

# first argument is label to use
label = sys.argv[1]

# iterate over passed files, which start at second argument
filenames = sys.argv[2: len(sys.argv)]

for filename in filenames:
	with open(filename) as f:
		# first line is key (not from pmdtools)
		filename = f.readline()
		key = os.path.splitext(os.path.basename(filename))[0] # filename without extension
		# second line is damage from pmdtools
		line = f.readline()
		# right fields are damage and error
		s = line.split(":")
		right = s[1]
		values = right.split()
		damage = float(values[0])
		stderr = values[1]

		print("%s\t%s\t%.3f" % (key, label, damage) )
