# Process the output of preseq to estimate the number of (additional) reads required for a library
# 

import sys

# raw_reads 
# number of distinct targets hit
def preseq_analysis(target_reads, unique_targets_uncorrected, number_raw_reads, total_target_hits_actual, unique_targets_actual, minimum_raw_reads, corrected_targets_per_raw_read_threshold):
	length = len(target_reads)
	if length != len(unique_targets_uncorrected):
		raise ValueError('length mismatch between entries in target_reads and unique_targets_uncorrected')
	
	inverse_e = 0.3678794
	tenth = 0.1
	
	reads_inverse_e = float('inf')
	reads_tenth = float('inf')
	total_reads_required = minimum_raw_reads
	expected_unique_targets_at_saturation = -1
	
	# not all demultiplexing reads merge, align, or map to targets
	# use a simple ratio to translate between reads for preseq and demultiplexing reads
	read_ratio = float(number_raw_reads) / total_target_hits_actual
	
	for i in range(1,length):
		slope_uncorrected_targets = float(unique_targets_uncorrected[i] - unique_targets_uncorrected[i-1]) / (target_reads[i] - target_reads[i-1])
		if slope_uncorrected_targets <= inverse_e:
			reads_inverse_e = min(reads_inverse_e, target_reads[i])
		if slope_uncorrected_targets <= tenth:
			reads_tenth = min(reads_tenth, target_reads[i])
			
		slope_corrected_targets = (target_hits_estimated_to_empirical(unique_targets_uncorrected[i]) - target_hits_estimated_to_empirical(unique_targets_uncorrected[i-1])) / ((target_reads[i] - target_reads[i-1]) * read_ratio)
		#print("{:f}\t{:f}\t{:f}\t{:f}".format(target_reads[i], target_reads[i] * read_ratio, unique_targets_uncorrected[i], slope_corrected_targets))
		
		if slope_corrected_targets > corrected_targets_per_raw_read_threshold:
			total_reads_required = max(total_reads_required, target_reads[i] * read_ratio)
			expected_unique_targets_at_saturation = target_hits_estimated_to_empirical(unique_targets_uncorrected[i])
		
			
	raw_reads_inverse_e = reads_inverse_e * read_ratio
	raw_reads_tenth = reads_tenth * read_ratio
	additional_reads_required = total_reads_required - number_raw_reads
	
	values = {
		'number_raw_reads': number_raw_reads,
		'unique_targets_actual': unique_targets_actual,
		'raw_reads_inverse_e': raw_reads_inverse_e,
		'raw_reads_tenth': raw_reads_tenth,
		'total_reads_required': total_reads_required,
		'additional_reads_required': additional_reads_required,
		'expected_unique_targets_at_saturation': expected_unique_targets_at_saturation
		}
	return values

# count number of reads overlapping any target and number of unique targets hit
def total_and_unique_target_hits(histogram_filename):
	with open(histogram_filename) as histogram:
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
	corrected_hit_faction = 1 / (1.1 + 1 / estimated_hit_fraction)
	corrected_hit_estimate = corrected_hit_faction * total_autosome_targets
	return corrected_hit_estimate

if __name__ == '__main__':
	histogram_filename = sys.argv[1]
	# empirically measured parameters to extrapolate from
	total_target_hits_actual, unique_targets_actual = total_and_unique_target_hits(histogram_filename)
	
	preseq_filename = sys.argv[2]
	number_raw_reads = int(sys.argv[3]) # number of raw reads that produced the histogram processed by preseq
	minimum_raw_reads = int(sys.argv[4]) # minimum raw reads required for a library
	corrected_targets_per_raw_read_threshold = float(sys.argv[5])

	target_reads = []
	unique_targets_uncorrected = []
	with open(preseq_filename) as preseq_file:
		# skip preseq header
		preseq_file.readline()
		
		# read in table of reads
		for line in preseq_file:
			reads, expected_distinct_reads, lower, upper = line.split()
			target_reads.append(float(reads))
			unique_targets_uncorrected.append(float(expected_distinct_reads))
			
	values = preseq_analysis(target_reads, unique_targets_uncorrected, number_raw_reads, total_target_hits_actual, unique_targets_actual, minimum_raw_reads, corrected_targets_per_raw_read_threshold)
	for key in values:
		print('{}\t{:.1f}'.format(key, values[key]))

