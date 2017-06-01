import os
import sys
# read in a VCF file and compute:
# number of reads for autosomal, X, Y chromosomes
# sex determination using # reads for chromosomes
# coverage as a fraction of targets

# first argument is label to use
label = sys.argv[1]
# second argument is VCF file
filename = sys.argv[2]

with open(filename) as f:
	key = os.path.basename(filename)
	while key.count('.') > 0:
		key = os.path.splitext(key)[0] # filename without extension
	
	targets = {}
	autosome_read_count = 0
	x_read_count = 0
	y_read_count = 0
	
	for line in f:
		#ignore lines that start with '#'
		if line[0] == '#':
			continue
		fields = line.split('\t')
		chromosome = fields[0]
		position = fields[1]
		
		# read depth (number of reads) from info field
		info = fields[7]
		info_sections = info.split(';')
		depth = 0
		for info_field in info_sections:
			if info_field.find('DP=') >= 0:
				depth = int(info_field.replace('DP=',''))
				break;
		
		target_key = chromosome + '-' + position
		#print(target_key, depth)
		# keep track of distinct targets
		if targets.has_key(target_key):
			targets[target_key] += depth
		else:
			targets[target_key] = depth
			
		nchromosome = int(chromosome)
		if 1 <= nchromosome and nchromosome <= 22:
			autosome_read_count += depth
		elif nchromosome == 23:
			x_read_count += depth
		elif nchromosome == 24:
			y_read_count += depth
		
		minimum_autosomal_read_hits_for_sex_determination = 50
		female_threshold = 0.01
		male_threshold = 0.05
		sex = 'U'
		if autosome_read_count >= minimum_autosomal_read_hits_for_sex_determination:
			if float(y_read_count) / autosome_read_count < female_threshold:
				sex = 'F'
			elif float(y_read_count) / autosome_read_count >= male_threshold:
				sex = 'M'

	print("%s\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%s" % (key,
													  label + '_autosome', autosome_read_count,
													  label + '_x', x_read_count, 
													  label + '_y', y_read_count,
													  label + 'sex', sex) 
	)
