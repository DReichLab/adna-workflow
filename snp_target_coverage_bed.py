import os
import sys

# Input: label and files containing autosome, X, and Y target counts
# Output: print to stdout keyed statistics for reporting with sex determination

# sex determination for 1240k is based on fraction of Y chromosome
def sex_determination(x_read_count, y_read_count):
	sex_chromosome_count = x_read_count + y_read_count
	minimum_sex_chromosome_count_for_sex_determination = 100
	female_threshold = 0.1
	male_threshold = 0.3
	sex = 'U'
	if sex_chromosome_count >= minimum_sex_chromosome_count_for_sex_determination:
		if float(y_read_count) / sex_chromosome_count <= female_threshold:
			sex = 'F'
		elif float(y_read_count) / sex_chromosome_count >= male_threshold:
			sex = 'M'
	return sex

if __name__ == '__main__':
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
		
	sex = sex_determination(x_read_count, y_read_count)
	
	print("%s\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%s" % (key,
														label + '_autosome', autosome_read_count,
														label + '_x', x_read_count, 
														label + '_y', y_read_count,
														label + '_sex', sex) 
)
