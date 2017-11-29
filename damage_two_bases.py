import os
import sys
# read in pmdtools input --d damage statistics from stdin

# 1-2-Q3-Q4.bam
#z	CT	CA	CG	CC	GA	GT	GC	GG
#0	0.01663	0.00039	0.00063	0.98236	0.02126	0.00096	0.00022	0.97755	
#1	0.00329	0.00062	0.00066	0.99544	0.00271	0.00079	0.0009	0.9956	

# first argument is label to use
label = sys.argv[1]

# iterate over passed files, which start at second argument
filename = sys.argv[2]

damage_values = {}
with open(filename) as f:
	sample_key = os.path.splitext(os.path.basename(filename))[0] # filename without extension
	# first line is header, check it for expected columns
	line = f.readline()
	fields = line.split('\t')
	if fields[1] != 'CT':
		raise ValueError('Unexpected header field. Found %s, expected CT' % (fields[1]))
	if fields[5] != 'GA':
		raise ValueError('Unexpected header field. Found %s, expected GA' % (fields[5]))
	
	# second line is damage at first base pmdtools
	line = f.readline()
	fields = line.split('\t')
	damage_values[label + '_ct1'] = float(fields[1])
	damage_values[label + '_ga1'] = float(fields[5])
	
	# third line is damage at second base
	line = f.readline()
	fields = line.split('\t')
	damage_values[label + '_ct2'] = float(fields[1])
	damage_values[label + '_ga2'] = float(fields[5])

	damage_line = sample_key
	for damage_key in damage_values:
		damage_line = '%s\t%s\t%.3f' % (damage_line, damage_key, damage_values[damage_key])
	print(damage_line)
 
