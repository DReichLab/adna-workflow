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
				   'MT_pre', 'MT_pre-coverageLength',
				   'MT_post', 'MT_post-coverageLength']

# read from stats
statsFilename = sys.argv[1]
samples = dict()

with open(statsFilename, "r") as f:
	line = f.readline()
	total_reads = int(line)
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

# read from damages

# print headers
print ('Index-Barcode Key', end='\t')
for header in headersToReport:
	print(header, end='\t')
print ('') # includes newline
# output in preset header order
for sampleID in samples:
	print(sampleID, end='\t')
	thisSample = samples[sampleID]
	for label in headersToReport:
		if label in thisSample:
			print(thisSample[label], end='')
		print('\t', end='')
	print('')
