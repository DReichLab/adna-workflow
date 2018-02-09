# Produces a depth histogram for evaluating 1240k complexity 
# Input
# rows containing 3 columns, either file based or from stdin:
# chromosome position(1-based) depth
# Output
# histogram of number of positions with each read depth

import fileinput

depth_counts = {}
max_depth = 0
for line in fileinput.input():
	fields = line.split()
	chromosome = fields[0]
	position = fields[1]
	depth = int(fields[2])
	
	max_depth = max(max_depth, depth)
	current_depth_count = depth_counts.get(depth, 0) + 1
	depth_counts[depth] = current_depth_count
		
# for preseq, ignore targets with 0 depth
for depth in range(1, max_depth+1):
	count = depth_counts.get(depth, 0)
	print('{:d}\t{:d}'.format(depth, count))
