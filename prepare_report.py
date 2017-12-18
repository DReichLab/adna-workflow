from __future__ import print_function
# combine statistics from samples with damage reports into a tab-separated file readable by MS Excel
import sys

headersToReport = ['sample_sheet_key',
				   'library_id',
				   'raw', 
				   'merged', 
				   'endogenous_pre',
				   'autosome_pre', 'autosome_pre-coverageLength', 
				   'autosome_post', 'autosome_post-coverageLength',
				   'X_pre', 'X_pre-coverageLength', 
				   'X_post', 'X_post-coverageLength',
				   'Y_pre', 'Y_pre-coverageLength',
				   'Y_post', 'Y_post-coverageLength',
				   'duplicates_hs37d5',
				   'median_hs37d5',
				   'mean_hs37d5',
				   'damage_hs37d5_ct1',
				   'damage_hs37d5_ct2',
				   'damage_hs37d5_ga1',
				   'damage_hs37d5_ga2',
				   'MT_pre', 'MT_pre-coverageLength',
				   'MT_post', 'MT_post-coverageLength',
				   'duplicates_rsrs',
				   'median_rsrs',
				   'mean_rsrs',
				   'damage_rsrs_ct1',
				   'damage_rsrs_ct2',
				   'damage_rsrs_ga1',
				   'damage_rsrs_ga2',
				   'Haplogroup',
				   'Haplogroup_rank',
				   'spike3k_pre_autosome',
				   'spike3k_pre_x',
				   'spike3k_pre_y',
				   'spike3k_post_autosome',
				   'spike3k_post_x',
				   'spike3k_post_y',
				   'spike3k_post_sex',
				   'spike3k_complexity',
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
				   'contamination_contammix_inferred_error'
				   ]

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
	libraryID = keyMapping.get(sampleID, '')
	sampleSheetID = ''
	if libraryID != '':
		sampleSheetID = sampleID
	else:
		i5, i7, p5_set, p7_set = sampleID.split('_')
		for p5 in p5_set.split(':'):
			for p7 in p7_set.split(':'):
				trialSampleID = '{}_{}_{}_{}'.format(i5, i7, p5, p7)
				trialLibraryID = keyMapping.get(trialSampleID, '')
				
				if libraryID == '':
					if trialLibraryID != '':
						libraryID = trialLibraryID
						sampleSheetID = trialSampleID
				# if there is more than one libraryID that matches, we have a nonprogramming problem
				elif trialLibraryID != '':
					libraryID = 'MULTIPLE'
					sampleSheetID = 'MULTIPLE'
					
	return sampleSheetID, libraryID

if __name__ == '__main__':
	# read from stats
	statsFilename = sys.argv[1]
	with open(statsFilename, "r") as f:
		line = f.readline()
		total_reads = int(line)
		print('Total reads: ', total_reads)
		addToSamples(f)
		
	# mapping from index and barcodes to sample/extract/library ID
	keyMapping = dict()
	keyMappingFilename = sys.argv[2]
	with open(keyMappingFilename, "r") as f:
		for line in f:
			index_barcode_key, sample_extract_library_id = line.split('\t')
			keyMapping[index_barcode_key] = sample_extract_library_id.strip()

	# read from damages, medians, haplogroups
	# these are all files with the index barcode keys and additional keyed fields
	filenames = sys.argv[3:len(sys.argv)]
	for filename in filenames:
		with open(filename, "r") as f:
			addToSamples(f)
			
	# populate additional sample fields
	for sampleID in samples:
		# add sample/extract/library ID, if available
		samples[sampleID]['sample_sheet_key'], samples[sampleID]['library_id'] = findSampleSheetEntry(sampleID, keyMapping)
		# add % endogenous
		singleSample = samples[sampleID]
		if ('autosome_pre' in singleSample
		or 'X_pre' in singleSample
		or 'Y_pre' in singleSample
		or 'MT_pre' in singleSample):
			samples[sampleID]['endogenous_pre'] = float(
				int(samples[sampleID].get('autosome_pre', '0'))
				+ int(samples[sampleID].get('X_pre', '0'))
				+ int(samples[sampleID].get('Y_pre', '0'))
				+ int(samples[sampleID].get('MT_pre', '0')) ) / int(samples[sampleID]['merged'])

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
			if int(samples[sampleID]['raw']) >= 500 or sampleID in keyMapping:
				printSample(sampleID, thisSample)
		except KeyError:
			eprint('KeyError', sampleID)
	# samples that are expected, but do not have results
	for sampleID in keyMapping:
		if sampleID not in samples:
			print('%s\t%s' % (sampleID, keyMapping[sampleID]) )
