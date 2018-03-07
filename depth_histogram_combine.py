# combine multiple depth histograms into a single depth histogram
# example usage:
# python depth_histogram_combine.py autosome autosome.histogram X X.histogram Y Y.histogram 

import sys

def read_histogram(filename):
	max_depth = 0
	depth_counts = {}
	with open(filename) as f:
		for line in f:
			depth, count = [int(x) for x in line.split()]
			if count > 0:
				max_depth = max(max_depth, depth)
			depth_counts[depth] = count
	return max_depth, depth_counts

if __name__ == '__main__':
	labels = sys.argv[1::2]
	depth_filenames = sys.argv[2::2]

	if len(labels) != len(depth_filenames):
		raise ValueError('There should be an equal number of labels ({:d}) and depth filenames ({:d})'.format(len(labels), len(depth_filenames)))
	
	max_depth = 0
	depth_histograms = []
	for filename in depth_filenames:
		depth, histogram = read_histogram(filename)
		max_depth = max(max_depth, depth)
		depth_histograms.append(histogram)
		
	# print header
	print('depth\t{}'.format('\t'.join(labels)))
	for depth in range(0, max_depth+1):
		counts = []
		for histogram in depth_histograms:
			counts.append(str(histogram.get(depth, 0)))
		print('{:d}\t{}'.format(depth, '\t'.join(counts) ))
