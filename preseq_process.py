# Process the output of preseq to estimate the number of (additional) reads required for a library
# 

from __future__ import print_function
import argparse

# raw_reads 
# number of distinct targets hit
def preseq_analysis(reads_hitting_any_target, unique_reads, number_raw_reads, total_reads_hitting_any_target_actual, unique_targets_hit, minimum_raw_reads, empiricalTargetEstimator):
	length = len(reads_hitting_any_target)
	if length != len(unique_reads):
		raise ValueError('length mismatch between entries in reads_hitting_any_target and unique_reads')
	
	inverse_e = 0.3678794
	tenth = 0.1
	
	reads_inverse_e = float('inf')
	reads_tenth = float('inf')
	
	# not all demultiplexing reads merge, align, or map to targets
	# use a simple ratio to translate between reads for preseq and demultiplexing reads
	read_ratio = 0
	if total_reads_hitting_any_target_actual > 0:
		read_ratio = float(number_raw_reads) / total_reads_hitting_any_target_actual
	#print(read_ratio)
	
	raw_reads = [x * read_ratio for x in reads_hitting_any_target]
	estimated_targets = [empiricalTargetEstimator.unique_reads_to_empirical_targets(x) for x in unique_reads]
	
	raw_reads_inverse_e, ignored = find_xy_for_slope(raw_reads, estimated_targets, inverse_e)
	raw_reads_tenth, ignored = find_xy_for_slope(raw_reads, estimated_targets, tenth)
	
	# ratio of unique targets to raw reads
	thresholds = ['0.01', '0.075', '0.005']
	total_reads_required = {}
	expected_unique_targets_at_threshold = {}
	additional_reads_required = {}
	
	for threshold in thresholds:
		raw_reads_threshold, expected_unique_targets_at_threshold[threshold] = find_xy_for_slope(raw_reads, estimated_targets, float(threshold))
		total_reads_required[threshold] = max(raw_reads_threshold, minimum_raw_reads)
		additional_reads_required[threshold] = max(total_reads_required[threshold] - number_raw_reads, 0)
	
	values = {
		'number_raw_reads': number_raw_reads,
		'preseq_unique_targets_hit': unique_targets_hit,
		'preseq_raw_reads_inverse_e': raw_reads_inverse_e,
		'preseq_raw_reads_tenth': raw_reads_tenth,
	}
	for threshold in thresholds:
		values['preseq_total_reads_required_' + threshold] = total_reads_required[threshold]
		values['preseq_additional_reads_required_' + threshold] = additional_reads_required[threshold]
		values['preseq_expected_unique_targets_at_threshold_' + threshold] = expected_unique_targets_at_threshold[threshold]
	return values

# assuming decreasing slope, find x, y such that the slope is the slope_threshold
def find_xy_for_slope(x, y, slope_threshold):
	if len(x) != len(y):
		raise ValueError('length mismatch')
	length = len(y)
	slopes = [float('inf')]
	for i in range(1,length):
		slope = (y[i] - y[i-1]) / (x[i] - x[i-1])
		slopes.append(slope)
	
	i = len(slopes)-1
	while (i > 0) and (slope_threshold > slopes[i]):
		i -= 1
	# edge cases
	if i == 0:
		return 0, 0
	elif i == len(slopes)-1:
		return x[i], y[i]
	
	slope_change = slopes[i+1] - slopes[i]
	x_change = x[i+1] - x[i]
	x_at_threshold = x[i] + (slope_threshold - slopes[i]) / slope_change * x_change
	# this is not exact: we are mixing a piecewise linear model with acceleration
	# just use the threshold value
	y_at_threshold = y[i] + slope_threshold * (x_at_threshold - x[i])
	
	return x_at_threshold, y_at_threshold

def interpolate(X, Y, x):
	if len(X) != len(Y):
		raise ValueError('length mismatch')
	# find position where x[i] <= x < x[i+1]
	i = len(X)-1
	while (i > 0) and (x < X[i]):
		i -= 1
	
	if i == len(X)-1:
		if x > X[i]:
			return '>{:.1f}'.format(Y[i])
		else:
			return Y[i]
	
	slope = (Y[i+1] - Y[i]) / (X[i+1] - X[i])
	x_change = x - X[i]
	y = Y[i] + slope * x_change
	return y

# count number of reads overlapping any target and number of unique targets hit
# This histogram is not the input to preseq. It counts the number of reads
# for each target, which allows counting the actual number of targets hit
def total_and_unique_target_hits(target_histogram_filename):
	total_hits = 0
	unique_targets = 0
	try:
		with open(target_histogram_filename) as histogram:
			for line in histogram:
				target_count, occurrences = map(lambda x: int(x), line.split())
				total_hits += target_count * occurrences
				unique_targets += occurrences
	except:
		pass
	return total_hits, unique_targets
	
# defines the empirical relationship between the unique reads output from preseq and the number of unique 1240k targets hit
# empirical targets = 1 / (a + b * (unique reads) ^ power)
class EmpiricalTargetEstimator:
	total_autosome_targets = 1150639 # autosome targets in the 1240k target set
	a = 1
	b = 1
	power = -1;
	
	def __init__(self, a, b, power):
		self.a = a
		self.b = b
		self.power = power
	
	# historical empirical relationship between unique reads hitting targets and total targets covered
	def unique_reads_to_empirical_targets(self, unique_reads):
		estimated_hit_fraction = max(float(unique_reads) / self.total_autosome_targets, 1e-9)
		corrected_hit_faction = 1 / (self.a + self.b * (estimated_hit_fraction ** self.power))
		corrected_hit_estimate = corrected_hit_faction * self.total_autosome_targets
		return corrected_hit_estimate
	
def read_preseq_file(filename):
	reads_hitting_any_target = []
	unique_reads = []
	try:
		with open(filename) as preseq_file:
			# skip preseq header
			preseq_file.readline()
			
			# read in table of reads
			for line in preseq_file:
				reads, expected_distinct_reads, lower, upper = line.split()
				reads_hitting_any_target.append(float(reads))
				unique_reads.append(float(expected_distinct_reads))
	except:
		pass
	return reads_hitting_any_target, unique_reads

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	# number of raw demultiplexing reads that produced the histogram processed by preseq
	parser.add_argument("-n", "--number_raw_reads", help="number of raw reads demultiplexing", type=int, required=True)
	# minimum raw reads required for a library
	parser.add_argument("-m", "--minimum_raw_reads", help="minimum required number of raw reads demultiplexing", type=int, default=3e6)
	parser.add_argument("-r", "--total_reads_hitting_any_target_actual", help="Reads hitting any autosome target, used for computing read ratio", type=int, required=True)
	#parser.add_argument("-t", "--expected_targets_per_raw_read_threshold", help="threshold for ratio of expected unique targets hit to raw demultiplexing reads", type=float, default=0.01)
	
	# empirical model parameters
	parser.add_argument("-a", "--modelParameterA", help="parameter for empirical model of unique reads to total targets covered", type=float, default=1)
	parser.add_argument("-b", "--modelParameterB", help="parameter for empirical model of unique reads to total targets covered", type=float, default=1)
	parser.add_argument("-p", "--modelParameterPower", help="exponent parameter for empirical model of unique reads to total targets covered", type=float, default=-1)

	parser.add_argument("target_histogram_filename", help="Target histogram to count number of actual target hits. This is not the input to preseq.")
	parser.add_argument("preseq_filename", help="preseq file containing mapping of reads hitting targets to unique reads hitting targets. Name this using the index-barcode key.")
	parser.add_argument("-k", "--key", help="Key for aggregating statistics. Normally this should be index-barcode key.", required=True)
	args = parser.parse_args()
	
	empiricalTargetEstimator = EmpiricalTargetEstimator(args.modelParameterA, args.modelParameterB, args.modelParameterPower)
	
	# empirically measured parameters to extrapolate from
	target_hits_ignored, unique_targets_hit = total_and_unique_target_hits(args.target_histogram_filename)

	# read preseq table from file
	reads_hitting_any_target, unique_reads = read_preseq_file(args.preseq_filename)
	
	values = preseq_analysis(reads_hitting_any_target, unique_reads, args.number_raw_reads, args.total_reads_hitting_any_target_actual, unique_targets_hit, args.minimum_raw_reads, empiricalTargetEstimator)
	
	# output
	print(args.key, end='\t')
	for key in values:
		print('{}\t{:.0f}'.format(key, values[key]), end='\t')

