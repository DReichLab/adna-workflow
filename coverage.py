# 1. Read the results of target coverage length
# 2. Normalize for reference length to produce customary coverage statistic
# Output is a series of lines with tab-separated fields: sample_id, coverage

# The input file is in the aDNA counting statistics format

import sys

filename = sys.argv[1]
reference_length = int(sys.argv[2])
sample_id_arg = sys.argv[3]
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
		
		if sampleID == sample_id_arg:
			for n in range(0, len(labels)):
				if field == labels[n]:
					coverage = float(values[n]) / reference_length
					print("{}\t{:f}".format(sample_id_arg, coverage))
					sys.exit(0)

# if the value does not appear, there is no coverage so print 0
print("{}\t{:f}".format(sample_id_arg, 0.0))
