# Prepare an input file for Nick's spike3k complexity program
# Input file is a tab delimited file with three columns
# 1. sample ID
# 2. pre-dup-removal autosomal targets
# 3. post-dup-removal-autosomal targets. 
# Output file adds a final column which is predicted complexity

# This takes an input of the spike3k pre deduplication and post deduplication
# 47-1-Q19-Q40    spike3k_pre_autosome    849     spike3k_pre_x   0       spike3k_pre_y   82      spike3k_pre_sex      M
# 47-1-Q19-Q41    spike3k_pre_autosome    850     spike3k_pre_x   0       spike3k_pre_y   81      spike3k_pre_sex      M

# 47-1-Q19-Q40    spike3k_post_autosome   400     spike3k_post_x  0       spike3k_post_y  40      spike3k_post_sex     M
# 47-1-Q19-Q41    spike3k_post_autosome   200     spike3k_post_x  0       spike3k_post_y  20      spike3k_post_sex     M

import sys
import os

autosome_prededpuplication = {}
autosome_postdeduplication = {}

# read in 
preduplication_filename = sys.argv[1]
with open(preduplication_filename, "r") as preduplication:
	for line in preduplication:
		fields = line.split('\t')
		key = fields[0]
		labels = fields[1:len(fields):2]
		values = fields[2:len(fields):2]
		
		for n in range(0, len(labels)):
			if labels[n] == 'spike3k_pre_autosome':
				autosome_prededpuplication[key] = int(values[n])
				
postduplication_filename = sys.argv[2]
with open(postduplication_filename, "r") as postduplication:
	for line in postduplication:
		fields = line.split('\t')
		key = fields[0]
		labels = fields[1:len(fields):2]
		values = fields[2:len(fields):2]
		
		for n in range(0, len(labels)):
			if labels[n] == 'spike3k_post_autosome':
				autosome_postdeduplication[key] = int(values[n])

for key in autosome_prededpuplication.keys():
	if autosome_postdeduplication.has_key(key):
		print("%s\t%d\t%d" % (key, autosome_prededpuplication[key], autosome_postdeduplication[key]) )
