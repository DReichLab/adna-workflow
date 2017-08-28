from __future__ import print_function
# combine statistics from samples with damage reports into a tab-separated file readable by MS Excel
import sys

headersToReport = ['raw', 
				   'merged', 
				   'autosome_pre', 'autosome_pre-coverageLength', 
				   'autosome_post', 'autosome_post-coverageLength',
				   'X_pre', 'X_pre-coverageLength', 
				   'X_post', 'X_post-coverageLength',
				   'Y_pre', 'Y_pre-coverageLength',
				   'Y_post', 'Y_post-coverageLength',
				   'duplicates_hs37d5',
				   'median_hs37d5',
				   'mean_hs37d5',
				   'damage_hs37d5',
				   'MT_pre', 'MT_pre-coverageLength',
				   'MT_post', 'MT_post-coverageLength',
				   'duplicates_rsrs',
				   'median_rsrs',
				   'mean_rsrs',
				   'damage_rsrs',
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
				   'contamination_schmutzi',
				   'contamination_schmutzi_lower',
				   'contamination_schmutzi_upper',
				   'contamination_rare_variant',
				   'contamination_rare_variant_lower',
				   'contamination_rare_variant_upper',
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
		sampleID = tokens[0]
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

# read from stats
statsFilename = sys.argv[1]
with open(statsFilename, "r") as f:
	line = f.readline()
	total_reads = int(line)
	print('Total reads: ', total_reads)
	addToSamples(f)

# read from damages, medians, haplogroups
filenames = sys.argv[2:len(sys.argv)]
for filename in filenames:
	with open(filename, "r") as f:
		addToSamples(f)

# print headers
print ('Index-Barcode Key', end='\t')
for header in headersToReport:
	headerToPrint = header
	# if a coverage label, shorten it
	if isCoverageLabel(header):
		headerToPrint, unusedValue = coverageNormalization(header, 0)
	print(headerToPrint, end='\t')
print ('') # includes newline
# output in preset header order
for sampleID in samples:
	thisSample = samples[sampleID]
	try:
		if int(samples[sampleID]['raw']) >= 500:
			print(sampleID, end='\t')
			for label in headersToReport:
				if label in thisSample:
					# adjust coverage values, and pass through non-coverage values
					value = thisSample[label]
					if isCoverageLabel(label):
						ignoredLabel, value = coverageNormalization(label, value)
						value = '%.4g' % value
					print(value, end='')
				print('\t', end='')
			print('')
	except KeyError:
		eprint('KeyError', sampleID)
