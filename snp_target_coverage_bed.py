import os
import sys

# The spike3k targets are now analyzed using not only the single SNP position, but also
# the surrounding region. This region is passed to samtools as a bed file. This takes the 
# counts for autosome, X, Y targets and estimates the sex and outputs this information
# keyed using the autosome filename with labels built from the first passed argument. 

# first argument is label to use
label = sys.argv[1]
# second through fourth arguments are counts of reads that align with targets
autosome_filename = sys.argv[2]
x_filename = sys.argv[3]
y_filename = sys.argv[4]

autosome_read_count = 0
x_read_count = 0
y_read_count = 0

key = os.path.basename(autosome_filename)
while key.count('.') > 0:
	key = os.path.splitext(key)[0] # filename without extension

with open(autosome_filename) as f:
	line = next(f)
	autosome_read_count = int(line)
	
with open(x_filename) as f:
	line = next(f)
	x_read_count = int(line)

with open(y_filename) as f:
	line = next(f)
	y_read_count = int(line)
		
			
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
													label + '_sex', sex) 
)
