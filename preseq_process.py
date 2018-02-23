# Process the output of preseq to estimate the number of (additional) reads required for a library
# 

import sys
import argparse

# raw_reads 
# number of distinct targets hit
def preseq_analysis(reads_hitting_any_target, unique_reads, number_raw_reads, total_reads_hitting_any_target_actual, unique_targets_hit, minimum_raw_reads, expected_targets_per_raw_read_threshold):
	length = len(reads_hitting_any_target)
	if length != len(unique_reads):
		raise ValueError('length mismatch between entries in reads_hitting_any_target and unique_reads')
	
	inverse_e = 0.3678794
	tenth = 0.1
	
	reads_inverse_e = float('inf')
	reads_tenth = float('inf')
	total_reads_required = minimum_raw_reads
	expected_unique_targets_at_threshold = -1
	
	# not all demultiplexing reads merge, align, or map to targets
	# use a simple ratio to translate between reads for preseq and demultiplexing reads
	read_ratio = float(number_raw_reads) / total_reads_hitting_any_target_actual
	#print(read_ratio)
	
	for i in range(1,length):
		#slope_uncorrected_targets = float(unique_reads[i] - unique_reads[i-1]) / (reads_hitting_any_target[i] - reads_hitting_any_target[i-1])
			
		slope_corrected_targets = (target_hits_estimated_to_empirical(unique_reads[i]) - target_hits_estimated_to_empirical(unique_reads[i-1])) / ((reads_hitting_any_target[i] - reads_hitting_any_target[i-1]) * read_ratio)
		#print("{:f}\t{:f}\t{:f}\t{:f}".format(reads_hitting_any_target[i], reads_hitting_any_target[i] * read_ratio, unique_reads[i], slope_corrected_targets))
		if slope_corrected_targets <= inverse_e:
			reads_inverse_e = min(reads_inverse_e, reads_hitting_any_target[i])
		if slope_corrected_targets <= tenth:
			reads_tenth = min(reads_tenth, reads_hitting_any_target[i])
		
		if slope_corrected_targets > expected_targets_per_raw_read_threshold:
			total_reads_required = max(total_reads_required, reads_hitting_any_target[i] * read_ratio)
			expected_unique_targets_at_threshold = target_hits_estimated_to_empirical(unique_reads[i])
		
			
	raw_reads_inverse_e = reads_inverse_e * read_ratio
	raw_reads_tenth = reads_tenth * read_ratio
	additional_reads_required = total_reads_required - number_raw_reads
	
	values = {
		'number_raw_reads': number_raw_reads,
		'unique_targets_hit': unique_targets_hit,
		'raw_reads_inverse_e': raw_reads_inverse_e,
		'raw_reads_tenth': raw_reads_tenth,
		'total_reads_required': total_reads_required,
		'additional_reads_required': additional_reads_required,
		'expected_unique_targets_at_threshold': expected_unique_targets_at_threshold
		}
	return values

# count number of reads overlapping any target and number of unique targets hit
# This histogram is not the input to preseq. It counts the number of reads
# for each target, which allows counting the actual number of targets hit
def total_and_unique_target_hits(target_histogram_filename):
	with open(target_histogram_filename) as histogram:
		total_hits = 0
		unique_targets = 0
		for line in histogram:
			target_count, occurrences = map(lambda x: int(x), line.split())
			total_hits += target_count * occurrences
			unique_targets += occurrences
		return total_hits, unique_targets
	
# historical empirical relationship between computed number of target hits and actual target hits
def target_hits_estimated_to_empirical(estimated_hits):
	total_autosome_targets = 1150639 # autosome targets in the 1240k target set
	estimated_hit_fraction = max(float(estimated_hits) / total_autosome_targets, 1e-9)
	corrected_hit_faction = 1 / (1 + 1 / estimated_hit_fraction)
	corrected_hit_estimate = corrected_hit_faction * total_autosome_targets
	return corrected_hit_estimate

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	# number of raw demultiplexing reads that produced the histogram processed by preseq
	parser.add_argument("-n", "--number_raw_reads", help="number of raw reads demultiplexing", type=int, required=True)
	# minimum raw reads required for a library
	parser.add_argument("-m", "--minimum_raw_reads", help="minimum required number of raw reads demultiplexing", type=int, default=3e6)
	parser.add_argument("-t", "--expected_targets_per_raw_read_threshold", help="threshold for ratio of expected unique targets hit to raw demultiplexing reads", type=float, default=0.01)
	
	parser.add_argument("target_histogram_filename", help="Target histogram to count number of actual target hits. This is not the input to preseq.")
	parser.add_argument("preseq_filename", help="preseq file containing mapping of reads hitting targets to unique reads hitting targets")
	args = parser.parse_args()
	
	# empirically measured parameters to extrapolate from
	total_reads_hitting_any_target_actual, unique_targets_hit = total_and_unique_target_hits(args.target_histogram_filename)

	reads_hitting_any_target = []
	unique_reads = []
	with open(args.preseq_filename) as preseq_file:
		# skip preseq header
		preseq_file.readline()
		
		# read in table of reads
		for line in preseq_file:
			reads, expected_distinct_reads, lower, upper = line.split()
			reads_hitting_any_target.append(float(reads))
			unique_reads.append(float(expected_distinct_reads))
			
	values = preseq_analysis(reads_hitting_any_target, unique_reads, args.number_raw_reads, total_reads_hitting_any_target_actual, unique_targets_hit, args.minimum_raw_reads, args.expected_targets_per_raw_read_threshold)
	for key in values:
		print('{}\t{:.1f}'.format(key, values[key]))

