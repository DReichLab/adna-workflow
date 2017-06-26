# 1. Read the results of target coverage length
# 2. Normalize for reference length to produce customary coverage statistic
# 3. Output in a manner suitable for reading back in WDL

# The input file is in the aDNA counting statistics format

import sys

filename = sys.argv[1]
reference_length = int(sys.argv[2])
id = sys.argv[3]
field = sys.argv[4]

with open(filename) as f:
	# skip the first line
	f.readline()
	for line in f:
		#key label1 value1 label2 value2 ...
		tokens = str.strip(line).split('\t')
		sampleID = tokens[0]
		labels = tokens[1:len(tokens):2]
		values = tokens[2:len(tokens):2]
		
		if sampleID == id:
			for n in range(0, len(labels)):
				if field == labels[n]:
					coverage = float(values[n]) / reference_length
					print(coverage)