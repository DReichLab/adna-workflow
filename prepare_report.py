# combine statistics from samples with damage reports into a tab-separated file readable by MS Excel
import sys
import math
import argparse

READS_THRESHOLD_TO_REPORT_KEY = 25000

headersToReport = [
	'sample_sheet_key',
	'library_id',
	'plate_id',
	'experiment',
	'i5',
	'i7',
	'p5_barcode',
	'p7_barcode',
	'sample_sheet_i5',
	'sample_sheet_i7',
	'sample_sheet_p5_barcode',
	'sample_sheet_p7_barcode',
				   'raw', 
				   'merged', 
				   'endogenous_pre',
				   'autosome_pre', 'autosome_pre-coverageLength', 
				   'autosome_post', 'autosome_post-coverageLength',
				   'X_pre', 'X_pre-coverageLength', 
				   'X_post', 'X_post-coverageLength',
				   'Y_pre', 'Y_pre-coverageLength',
				   'Y_post', 'Y_post-coverageLength',
				   'duplicates_nuclear',
				   'median_nuclear',
				   'mean_nuclear',
				   'damage_nuclear_ct1',
				   'damage_nuclear_ct2',
				   'damage_nuclear_ga1',
				   'damage_nuclear_ga2',
				   'MT_pre', 'MT_pre-coverageLength',
				   'MT_post', 'MT_post-coverageLength',
				   'duplicates_rsrs',
				   'median_rsrs',
				   'mean_rsrs',
				   'damage_rsrs_ct1',
				   'damage_rsrs_ct2',
				   'damage_rsrs_ga1',
				   'damage_rsrs_ga2',
				   'MT_Haplogroup',
				   'MT_Haplogroup_rank',
				   'MT_Haplogroup_NotFoundPolys',
				   'MT_Haplogroup_FoundPolys',
				   'MT_Haplogroup_RemainingPolys',
				   #spike3k is obsolete
				   #'spike3k_pre_autosome',
				   #'spike3k_pre_x',
				   #'spike3k_pre_y',
				   #'spike3k_post_autosome',
				   #'spike3k_post_x',
				   #'spike3k_post_y',
				   #'spike3k_post_sex',
				   #'spike3k_complexity',
				   #'recommendation_spike3k',
#				   'contamination_schmutzi',
#				   'contamination_schmutzi_lower',
#				   'contamination_schmutzi_upper',
#				   'contamination_rare_variant',
#				   'contamination_rare_variant_lower',
#				   'contamination_rare_variant_upper',
				   'contamination_contammix',
				   'contamination_contammix_lower',
				   'contamination_contammix_upper',
				   'contamination_contammix_gelman',
				   'contamination_contammix_inferred_error',
				   '1240k_pre_autosome',
				   '1240k_post_autosome',
				   '1240k_pre_x',
				   '1240k_post_x',
				   '1240k_pre_y',
				   '1240k_post_y',
				   '1240k_post_sex',
				   'angsd_nsites',
				   'angsd_MoM',
				   'angsd_SE(MoM)',
				   'angsd_ML',
				   'angsd_SE(ML)',
				   'preseq_unique_targets_hit',
				   'preseq_raw_reads_inverse_e',
				   'preseq_raw_reads_tenth',
				   'preseq_coverage_at_marginal_uniqueness_0.10',
				   ]
preseq_unique_target_per_raw_read_thresholds = ['0.01', '0.0075', '0.005']
for threshold in preseq_unique_target_per_raw_read_thresholds:
	headersToReport.append('preseq_total_reads_required_{}'.format(threshold))
	headersToReport.append('preseq_additional_reads_required_{}'.format(threshold))
	headersToReport.append('preseq_expected_unique_targets_at_threshold_{}'.format(threshold))
	headersToReport.append('preseq_target_coverage_at_threshold_{}'.format(threshold))

samples = dict()

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def addToSamples(f):
	for line in f:
		#key label1 value1 label2 value2 ...
		tokens = str.strip(line).split('\t')
		# to workaround Cromwell problem where task output filenames cannot contain ':'
		# we re-replace alternate character with ':' for analysis
		sampleID = tokens[0].replace('-', ':')
		labels = tokens[1:len(tokens):2]
		values = tokens[2:len(tokens):2]
		
		if sampleID not in samples:
			samples[sampleID] = {}
		
		for n in range(0, len(labels)):
			samples[sampleID][labels[n]] = values[n]
			
def isCoverageLabel(label):
	return "coverageLength" in label
				
# if a coverage label, normalize the coverage length sum to a coverage multiplier
# otherwise, do nothing
def coverageNormalization(originalLabel, value):
	label = originalLabel
	length = 1
	if isCoverageLabel(label):
		label = originalLabel.replace('Length', '')
		if label.startswith('MT'):
			length = 16569
		elif label.startswith('autosome'):
			length = 2881033286
		elif label.startswith('X'):
			length = 155270560
		elif label.startswith('Y'):
			length = 59373566
	return label, (float(value) / length)

def printSample(sampleID, thisSample):
	print(sampleID, end='\t')
	for label in headersToReport:
		if thisSample is not None and label in thisSample:
			# adjust coverage values, and pass through non-coverage values
			value = thisSample[label]
			if isCoverageLabel(label):
				ignoredLabel, value = coverageNormalization(label, value)
				value = '%.4g' % value
			print(value, end='')
		print('\t', end='')
	print('')

# match a barcode pair with the sample sheet barcodes
# each of the Q barcodes comprises 4 7-base-pair sequences
# for example: Q2 = {13, 14, 15, 16}
# the demultiplexing will use the Q sequences, so 13 will appear as Q2
# if the sample sheet uses 13, match a Q2 with 13
# returns:
# 	1: index-barcode 4 tuple
#	2: sample/extract/library id if any, or empty string
def findSampleSheetEntry(sampleID, keyMapping):
	libraryID, plateID, experiment = keyMapping.get(sampleID, ['','',''])
	sampleSheetID = ''
	if libraryID != '':
		sampleSheetID = sampleID
	else:
		i5, i7, p5_set, p7_set = sampleID.split('_')
		for p5 in p5_set.split(':'):
			for p7 in p7_set.split(':'):
				trialSampleID = '{}_{}_{}_{}'.format(i5, i7, p5, p7)
				trialLibraryID, trialPlateID, trialExperiment = keyMapping.get(trialSampleID, ['','',''])
				
				if libraryID == '':
					if trialLibraryID != '':
						libraryID = trialLibraryID
						sampleSheetID = trialSampleID
						plateID = trialPlateID
						experiment = trialExperiment
				# if there is more than one libraryID that matches, we have a nonprogramming problem
				elif trialLibraryID != '':
					MULTIPLE = 'MULTIPLE'
					libraryID = MULTIPLE
					sampleSheetID = MULTIPLE
					plateID = MULTIPLE
					experiment = MULTIPLE
					
	return sampleSheetID, libraryID, plateID, experiment

# Make an initial recommendation for whether this sample should continue in processing based on spike3k metrics, assuming it is UDG-half treated
def recommendation_spike3k(sample):
	# There are five possible outcomes
	PENDING = -2 # contamination estimate is not complete
	LOW_DATA = -1
	FAIL = 0
	WARNING = 1
	PASS = 2
	
	# if any criterion is poor, we downgrade recommendation
	recommendation_value = PASS
	
	if int(sample.get('merged', 0)) < 1000:
		return 'low data'
	
	# damage assuming UDG-half treatment
	DAMAGE_FAIL = 0.01
	DAMAGE_WARN = 0.03
	damage_first_base = (float(sample.get('damage_rsrs_ct1', -1.0)) + float(sample.get('damage_rsrs_ga1', -1.0)))/2
	if damage_first_base < 0:
		recommendation_value = PENDING
	elif damage_first_base < DAMAGE_FAIL:
		recommendation_value = min(recommendation_value, FAIL)
	elif damage_first_base < DAMAGE_WARN:
		recommendation_value = min(recommendation_value, WARNING)
	# contamination
	CONTAMINATION_FAIL = 0.6
	CONTAMINATION_WARN = 0.95
	contamination_fraction_matching_consensus = float(sample.get('contamination_contammix', 0))
	if math.isnan(contamination_fraction_matching_consensus) or contamination_fraction_matching_consensus < CONTAMINATION_FAIL:
		recommendation_value = min(recommendation_value, FAIL)
	elif contamination_fraction_matching_consensus < CONTAMINATION_WARN:
		recommendation_value = min(recommendation_value, WARNING)
	# MT coverage
	COVERAGE_FAIL = 0.5
	COVERAGE_WARN = 2.0
	label = 'MT_post-coverageLength'
	shortened_label, mt_coverage = coverageNormalization(label, sample.get(label, 0))
	if mt_coverage < COVERAGE_FAIL:
		recommendation_value = min(recommendation_value, FAIL)
	elif mt_coverage < COVERAGE_WARN:
		recommendation_value = min(recommendation_value, WARNING)
	# spike coverage estimate
	COMPLEXITY_FAIL = 0.005
	COMPLEXITY_WARN = 0.05
	spike_complexity = float(sample.get('spike3k_complexity', 0))
	if spike_complexity < COMPLEXITY_FAIL:
		recommendation_value = min(recommendation_value, FAIL)
	elif spike_complexity < COMPLEXITY_WARN:
		recommendation_value = min(recommendation_value, WARNING)
	# failed samples that have substantial good data are promoted to warning
	RESCUE_PRODUCT = 0.05
	if spike_complexity * damage_first_base > RESCUE_PRODUCT:
		recommendation_value = max(recommendation_value, WARNING)
	
	# return a string representing the recommendation
	if recommendation_value == FAIL:
		return 'fail'
	elif recommendation_value == WARNING:
		return 'warning'
	elif recommendation_value == PASS:
		return 'pass'
	elif recommendation_value == PENDING:
		return 'pending'
	
# parse an index-barcode key into labels
def parse_index_barcode_key_into_labels(key, i5_labels, i7_labels, barcode_labels):
	i5 = i7 = p5 = p7 = ''
	try:
		i5_sequence, i7_sequence, p5_sequence_set, p7_sequence_set = key.split('_')
		i5 = sequence_to_label(i5_sequence, i5_labels)
		i7 = sequence_to_label(i7_sequence, i7_labels)
		p5 = sequence_to_label(p5_sequence_set, barcode_labels)
		p7 = sequence_to_label(p7_sequence_set, barcode_labels)
	except:
		pass
	return i5, i7, p5, p7

# if any component sequences are not in the dictionary, then return the sequence
def sequence_to_label(sequence, labels):
	if sequence in labels:
		return labels[sequence]
	else:
		return sequence

# return a map of sequences to labels
# Labels are much easier to read than the sequences
def sequence_labels(filename):
	mapping = {}
	try:
		with open(filename) as f:
			for line in f:
				fields = line.split()
				sequence_string = fields[0]
				label = fields[1]
				mapping[sequence_string] = label # set (or singleton) is always labeled
				if ':' in sequence_string: # for Q barcodes, each individual sequence has a name
					sequences = sequence_string.split(':')
					for i in range(len(sequences)):
						mapping[sequences[i]] = fields[2+i]
	except:
		pass
	return mapping

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Combine ancient DNA analysis outputs into a single file keyed by index-barcode key")
	parser.add_argument("--barcode_labels_filename", help="File containing mapping of barcode sequences (P5 and P7) to labels")
	parser.add_argument("--i5_labels_filename", help="File containing mapping of i5 sequences to labels")
	parser.add_argument("--i7_labels_filename", help="File containing mapping of i7 sequences to labels")
	parser.add_argument("statistics_filename", help="File containing read count total and per sample")
	parser.add_argument("index_barcode_filename", help="Map from index-barcode key to library ID, plate ID, and experiment")
	parser.add_argument("keyed_values", help="Any number of files containing keyed values. First column is index-barcode key, followed by pairs of (field name, value) columns", nargs='*')
	args = parser.parse_args()
	
	barcode_labels = sequence_labels(args.barcode_labels_filename)
	i5_labels = sequence_labels(args.i5_labels_filename)
	i7_labels = sequence_labels(args.i7_labels_filename)
	
	# read from stats
	statsFilename = args.statistics_filename
	with open(statsFilename, "r") as f:
		line = f.readline()
		total_reads = int(line)
		print('Total reads: ', total_reads)
		addToSamples(f)
		
	# mapping from index and barcodes to sample/extract/library ID and plate ID
	keyMapping = dict()
	keyMappingFilename = args.index_barcode_filename
	with open(keyMappingFilename, "r") as f:
		for line in f:
			fields = [x.strip() for x in line.split('\t')]
			index_barcode_key = fields[0]
			sample_extract_library_id = fields[1]
			plate_id = fields[2]
			experiment = fields[3]
			keyMapping[index_barcode_key] = [sample_extract_library_id, plate_id, experiment]

	# read from damages, medians, haplogroups
	# these are all files with the index barcode keys and additional keyed fields
	filenames = args.keyed_values
	for filename in filenames:
		with open(filename, "r") as f:
			addToSamples(f)
			
	# populate additional sample fields
	for sampleID in samples:
		singleSample = samples[sampleID]
		# parse index-barcode fields into components
		singleSample['i5'], singleSample['i7'], singleSample['p5_barcode'], singleSample['p7_barcode'] = parse_index_barcode_key_into_labels(sampleID, i5_labels, i7_labels, barcode_labels)
		# add sample/extract/library ID, if available
		samples[sampleID]['sample_sheet_key'], samples[sampleID]['library_id'], samples[sampleID]['plate_id'], samples[sampleID]['experiment'] = findSampleSheetEntry(sampleID, keyMapping)
		singleSample['sample_sheet_i5'], singleSample['sample_sheet_i7'], singleSample['sample_sheet_p5_barcode'], singleSample['sample_sheet_p7_barcode'] = parse_index_barcode_key_into_labels(singleSample.get('sample_sheet_key', ''), i5_labels, i7_labels, barcode_labels)
		# add % endogenous
		
		if ('autosome_pre' in singleSample
		or 'X_pre' in singleSample
		or 'Y_pre' in singleSample
		or 'MT_pre' in singleSample):
			samples[sampleID]['endogenous_pre'] = float(
				int(samples[sampleID].get('autosome_pre', '0'))
				+ int(samples[sampleID].get('X_pre', '0'))
				+ int(samples[sampleID].get('Y_pre', '0'))
				+ int(samples[sampleID].get('MT_pre', '0')) ) / int(samples[sampleID]['merged'])
		# add recommendation concerning future processing
		#singleSample['recommendation_spike3k'] = recommendation_spike3k(singleSample)

	# print headers
	print ('Index-Barcode Key', end='\t')
	for header in headersToReport:
		headerToPrint = header
		# if a coverage label, shorten it
		if isCoverageLabel(header):
			headerToPrint, unusedValue = coverageNormalization(header, 0)
		print(headerToPrint, end='\t')
	print ('') # includes newline

	sorted_samples = sorted(samples, key=lambda x: int(samples[x]['raw']), reverse=True)
	# output each sample with data using preset header order
	for sampleID in sorted_samples:
		thisSample = samples[sampleID]
		try:
			if int(samples[sampleID]['raw']) >= READS_THRESHOLD_TO_REPORT_KEY or sampleID in keyMapping:
				printSample(sampleID, thisSample)
		except KeyError:
			eprint('KeyError', sampleID)
	# samples that are expected, but do not have results
	for sampleID in keyMapping:
		if sampleID not in samples:
			libraryID, plateID, experiment = keyMapping[sampleID]
			print('{0}\t{0}\t{1}\t{2}\t{3}'.format(sampleID, libraryID, plateID, experiment) )
