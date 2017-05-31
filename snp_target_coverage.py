import os
import sys
# read in a VCF file
# this assumes there is at least one read per position reported

# first argument is label to use
label = sys.argv[1]
# second argument is number of targets, to normalize coverage
number_of_targets = float(sys.argv[2])
# third argument is VCF file
filename = sys.argv[3]

with open(filename) as f:
	key = os.path.basename(filename)
	while key.count('.') > 0:
		key = os.path.splitext(key)[0] # filename without extension
	
	targets = {}
	for line in f:
		#ignore lines that start with '#'
		if line[0] == '#':
			continue
		fields = line.split('\t')
		chromosome = fields[0]
		position = fields[1]
		target_key = chromosome + '-' + position
		#print(target_key)
		if targets.has_key(target_key):
			targets[target_key] += 1
		else:
			targets[target_key] = 1

	coverage = len(targets) / number_of_targets
	print("%s\t%s\t%.3f" % (key, label, coverage) )
