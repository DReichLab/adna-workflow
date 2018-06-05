# Process the output of preseq to estimate the number of (additional) reads required for a library
# 

from __future__ import print_function
import argparse

# preseq projects the number of unique reads given a number of reads, given a histogram of reads at a particular point
# We analyze four different quantities
# 1. raw_reads (reads demultiplexing for index-barcode key)
# 2. reads hitting a 1240k autosome target
# 3. unique reads hitting a 1240k autosome target
# 4. number of distinct 1240k targets hit

# 1 is assumed to be related to 2 by a constant ratio
# preseq computes the function between 2 and 3
# 4 is computed from 3 using an empirical model
def preseq_analysis(reads_hitting_any_target, unique_reads, number_raw_reads, total_reads_hitting_any_target_actual, unique_targets_hit, minimum_raw_reads, empiricalTargetEstimator):
	length = len(reads_hitting_any_target)
	if length != len(unique_reads):
		raise ValueError('length mismatch between entries in reads_hitting_any_target and unique_reads')
	
	inverse_e = 0.3678794
	tenth = 0.1
	
	raw_reads_inverse_e = float('inf')
	raw_reads_tenth = float('inf')
	
	# not all demultiplexing reads merge, align, or map to targets
	# use a simple ratio to translate between reads for preseq and demultiplexing reads
	read_ratio = 0
	if total_reads_hitting_any_target_actual > 0:
		read_ratio = float(number_raw_reads) / total_reads_hitting_any_target_actual
	#print(read_ratio)
	
	raw_reads = [x * read_ratio for x in reads_hitting_any_target]
	estimated_targets = [empiricalTargetEstimator.unique_reads_to_empirical_targets(x) for x in unique_reads]
	
	try:
		raw_reads_inverse_e, ignored = find_xy_for_slope(raw_reads, estimated_targets, inverse_e)
		raw_reads_tenth, ignored = find_xy_for_slope(raw_reads, estimated_targets, tenth)
	except:
		pass
	
	# ratio of unique targets to raw reads
	unique_target_thresholds = ['0.01', '0.0075', '0.005']
	
	values = {
		'number_raw_reads': number_raw_reads,
		'preseq_unique_targets_hit': unique_targets_hit,
		'preseq_raw_reads_inverse_e': raw_reads_inverse_e,
		'preseq_raw_reads_tenth': raw_reads_tenth,
	}
	# fill out values so something is always returned
	for threshold in unique_target_thresholds:
		raw_reads_threshold = 0
		expected_unique_targets_at_threshold = -1
		try:
			raw_reads_threshold, expected_unique_targets_at_threshold = find_xy_for_slope(raw_reads, estimated_targets, float(threshold))
		except:
			pass
		total_reads_required = max(raw_reads_threshold, minimum_raw_reads)
		additional_reads_required = max(total_reads_required - number_raw_reads, 0)
		
		# if we are already at 70% of raw reads goal, then stop
		if (additional_reads_required / total_reads_required) < 0.3:
			additional_reads_required = 0
		
		unique_reads_hitting_any_target_at_threshold = 0
		try:
			unique_reads_hitting_any_target_at_threshold = interpolate(raw_reads, unique_reads, raw_reads_threshold)
		except:
			pass
	
		values['preseq_total_reads_required_' + threshold] = total_reads_required
		values['preseq_additional_reads_required_' + threshold] = additional_reads_required
		values['preseq_expected_unique_targets_at_threshold_' + threshold] = expected_unique_targets_at_threshold
		values['preseq_target_coverage_at_threshold_' + threshold] = unique_reads_hitting_any_target_at_threshold / empiricalTargetEstimator.total_autosome_targets
		
		if len(raw_reads) > 0 and raw_reads_threshold == raw_reads[-1]: # outside of preseq projections
			values['preseq_total_reads_required_' + threshold] = '>{:.0f}'.format(values['preseq_total_reads_required_' + threshold])
			values['preseq_additional_reads_required_' + threshold] = '>{:.0f}'.format(values['preseq_additional_reads_required_' + threshold])
	
	# marginal uniqueness at current number of reads
	try:
		expected_distinct_reads, marginal_uniqueness = interpolate_base(reads_hitting_any_target, unique_reads, total_reads_hitting_any_target_actual)
		values['preseq_marginal_uniqueness'] = marginal_uniqueness
	except:
		pass
	# marginal uniqueness threshold
	marginal_uniqueness_thresholds = ['0.368', '0.10']
	for threshold in marginal_uniqueness_thresholds:
		reads_hitting_any_target_at_threshold = 0
		unique_reads_hitting_any_target_at_threshold = 0
		try:
			reads_hitting_any_target_at_threshold, unique_reads_hitting_any_target_at_threshold = find_xy_for_slope(reads_hitting_any_target, unique_reads, float(threshold))
		except:
			pass
		
		coverage = unique_reads_hitting_any_target_at_threshold / empiricalTargetEstimator.total_autosome_targets
		values['preseq_coverage_at_marginal_uniqueness_' + threshold] = coverage if (len(reads_hitting_any_target) > 0) and (reads_hitting_any_target_at_threshold < reads_hitting_any_target[-1]) else ('>{:.3f}'.format(coverage))
		
	return values

# assuming decreasing slope, find x, y such that the slope is the slope_threshold
# if slope is already less than threshold, report first value
def find_xy_for_slope(X, Y, slope_threshold):
	if len(X) != len(Y):
		raise ValueError('length mismatch')
	if len(X) < 2:
		raise ValueError('no slope can be computed')
	length = len(Y)
	slopes = [float('inf')]
	for i in range(1,length):
		slope = (Y[i] - Y[i-1]) / (X[i] - X[i-1])
		slopes.append(slope)
	
	i = len(slopes)-1
	while (i > 0) and (slope_threshold > slopes[i]):
		i -= 1
	# edge cases
	if i == 0:
		return X[0], Y[0]
	elif i == len(slopes)-1:
		return X[i], Y[i]
	
	slope_change = slopes[i+1] - slopes[i]
	x_change = X[i+1] - X[i]
	x_at_threshold = X[i] + (slope_threshold - slopes[i]) / slope_change * x_change
	# this is not exact: we are mixing a piecewise linear model with acceleration
	# just use the threshold value
	y_at_threshold = Y[i] + slope_threshold * (x_at_threshold - X[i])
	
	return x_at_threshold, y_at_threshold

def interpolate(X, Y, x):
	y, slope = interpolate_base(X, Y, x)
	return y

# return interpolated value of y at x and slope
# slope is defined on interval [ x_i, x_{i+1} )
def interpolate_base(X, Y, x):
	if len(X) != len(Y):
		raise ValueError('length mismatch')
	if len(X) == 0:
		raise ValueError('no data to interpolate')
	# find position where x[i] <= x < x[i+1]
	i = len(X)-1
	while (i > 0) and (x < X[i]):
		i -= 1
	
	if i == len(X)-1:
		if x > X[i]:
			return '>{:.1f}'.format(Y[i]), None
		else:
			return Y[i], None
	
	slope = (Y[i+1] - Y[i]) / (X[i+1] - X[i])
	x_change = x - X[i]
	y = Y[i] + slope * x_change
	return y, slope

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
		print(key, end='\t')
		if str(values[key]).startswith('>'):
			print('{}'.format(values[key]), end='\t')
		elif 'coverage' in key:
			print('{:.3f}'.format(values[key]), end='\t')
		else:
			print('{:.0f}'.format(values[key]), end='\t')

