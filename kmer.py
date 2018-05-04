import argparse
import sys
from collections import Counter

from prepare_report import parse_index_barcode_key_into_labels, sequence_to_label, sequence_labels

# Analyze the distribution of barcode pairs for each index pair
# Each index pair is output on a row, with the 10 most common barcode pairs in descending order of counts

DELIMITER = '_'

def parse_index_barcode_key(index_barcode_key):
	i5, i7, p5, p7 = index_barcode_key.split(DELIMITER)
	index_pair = i5 + DELIMITER + i7
	barcode_pair = p5 + DELIMITER + p7
	return index_pair, barcode_pair

# compute a fraction, where the denominator is superset of the numerator
# if denominator is zero, the fraction is zero
# if the numerator is greater than the denominator, this is an error
def fraction_part(numerator, denominator):
	if denominator > 0:
		return numerator / denominator
	elif numerator == 0:
		return 0
	if numerator > denominator:
		raise ValueError('invalid fraction')

def output_kmer_line(index_pair, barcode_pair=None, sample_sheet_entry=None):
	library_id = sample_sheet_entry.library_id if sample_sheet_entry != None else ''
	experiment = sample_sheet_entry.experiment if sample_sheet_entry != None else ''
	print(library_id, end='\t')
	print(experiment, end='\t')
	
	if index_pair in index_pairs and barcode_pair in index_pairs[index_pair]:
		demultiplexed_reads = index_pairs[index_pair][barcode_pair]
	else:
		demultiplexed_reads = 0
	print("{:d}".format(demultiplexed_reads), end='\t')
	
	index_pair_count = index_pair_counts.get(index_pair, 0)
	print("{:.3f}".format(fraction_part(demultiplexed_reads, index_pair_count)), end='\t')
	
	barcode_pair_counts = index_pairs.get(index_pair, {})
	sorted_barcode_pairs = sorted(barcode_pair_counts, key=lambda x: barcode_pair_counts[x], reverse=True)
	i5, i7 = index_pair.split(DELIMITER)
	# statistics for this index pair
	print('{}\t{}\t{}\t{}\t'.format(i5, sequence_to_label(i5, i5_labels), i7, sequence_to_label(i7, i7_labels)), end='')
	print('{:d}\t{:.4f}'.format(index_pair_count, float(index_pair_count) / total_reads), end='')
	# statistics for top 10 barcodes for this index pair
	for i in range(0, number_top_barcode_pairs):
		barcode_pair = sorted_barcode_pairs[i] if i < len(sorted_barcode_pairs) else DELIMITER
		barcode_pair_count = barcode_pair_counts.get(barcode_pair, 0)
		index_barcode_key = '{}_{}'.format(index_pair, barcode_pair)
		
		sample_sheet_entry = sample_sheet.get(index_barcode_key, SampleSheetEntry('', ''))
		libraryID = sample_sheet_entry.library_id
		experiment = sample_sheet_entry.experiment
		p5, p7 = barcode_pair.split(DELIMITER)
		fraction = fraction_part(barcode_pair_count, index_pair_count)
		print ('\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:d}'.format(libraryID, experiment, p5, barcodes.get(p5, ''), p7, barcodes.get(p7, ''), fraction, barcode_pair_count), end='')
	# unknown barcodes
	unknown_barcode_pair_counts = unknown_barcodes.get(index_pair, Counter()).most_common(number_top_barcode_pairs)
	for i in range(0, number_top_barcode_pairs):
		barcode_pair, barcode_pair_count = unknown_barcode_pair_counts[i] if i < len(unknown_barcode_pair_counts) else (DELIMITER, 0)
		p5, p7 = barcode_pair.split(DELIMITER)
		fraction = fraction_part(barcode_pair_count, index_pair_count)
		print('\t{}\t{}\t{:.3f}\t{:d}'.format(p5, p7, fraction, barcode_pair_count), end='')
	print('')
	
class SampleSheetEntry:
	def __init__(self, library_id, experiment):
		self.library_id = library_id
		self.experiment = experiment

index_pairs = dict()
index_pair_counts = Counter()
sample_sheet = dict()
unknown_barcodes = dict()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Produce kmer analysis to check demultiplexing results")
	parser.add_argument("barcode_labels_filename", help="File containing mapping of barcode sequences (P5 and P7) to labels")
	parser.add_argument("i5_labels_filename", help="File containing mapping of i5 sequences to labels")
	parser.add_argument("i7_labels_filename", help="File containing mapping of i7 sequences to labels")
	parser.add_argument("counts_filename", help="Statistics file containing count of reads for each key")
	parser.add_argument("sample_sheet", help="Sample sheet filename mapping index-barcode keys to library IDs")
	parser.add_argument("unknown_barcodes", help="File containing unknown barcodes with counts by index pair")
	args = parser.parse_args()
	parser.add_argument('-t', "--threshold_to_report", help="threshold number of reads to include an index pair in the kmer analysis report", type=int, default=25000, nargs=1)
	args = parser.parse_args()
	
	barcodes = sequence_labels(args.barcode_labels_filename)
	i5_labels = sequence_labels(args.i5_labels_filename)
	i7_labels = sequence_labels(args.i7_labels_filename)
	counts_filename = args.counts_filename
	sample_sheet_filename = args.sample_sheet
	unknown_barcodes_filename = args.unknown_barcodes

	# make one pass through all data to count the number of reads for each index pair
	with open(counts_filename) as f:
		line = f.readline() # total reads
		total_reads = int(line)
		for line in f:
			fields = line.split('\t')
			index_barcode_key = fields[0]
			raw_label = fields[1]
			raw_count = int(fields[2])
			if raw_label != 'raw':
				raise ValueError("expected 'raw', but found %s", (raw_label) )
			
			index_pair, barcode_pair = parse_index_barcode_key(index_barcode_key)
			# create an entry for this index pair
			if index_pair not in index_pairs:
				index_pairs[index_pair] = Counter()
				
			# accounting for barcode pair within index pair
			index_pair_counts[index_pair] += int(raw_count)
			index_pairs[index_pair][barcode_pair] = int(raw_count)
			
	# sample sheet
	with open(sample_sheet_filename) as f:
		for line in f:
			fields = [x.strip() for x in line.split('\t')]
			sampleID = fields[0]
			libraryID = fields[1]
			batchID = fields[2]
			experiment = fields[3]
			sample_sheet[sampleID] = SampleSheetEntry(libraryID, experiment)
			
	# unknown barcodes: these are the barcodes that appear most frequently among reads that do not have known barcodes
	with open(unknown_barcodes_filename) as f:
		for line in f:
			fields = line.strip().split('\t')
			index_barcode_key = fields[0].replace('-', ':') # 
			index_pair, barcode_pair = parse_index_barcode_key(index_barcode_key)
			if index_pair not in unknown_barcodes:
				unknown_barcodes[index_pair] = Counter()
				
			unknown_barcode_pairs = fields[1::2]
			unknown_barcode_counts = [int(x) for x in fields[2::2]]
			for i in range(len(unknown_barcode_pairs)):
				unknown_barcodes[index_pair][unknown_barcode_pairs[i]] = unknown_barcode_counts[i]
			
	number_top_barcode_pairs = 10
	# print field headers
	print ('library id\texperiment\tdemultiplexed reads\tfrac index pair reads\ti5 index\ti5 label\ti7 index\ti7 label\tnum reads for index pair\tindex pair frac total run\t', end='')
	for i in range(1, number_top_barcode_pairs + 1):
		print('top library ID {0}\texperiment {0}\ttop barcode pair {0} p5\tp5 name\ttop barcode pair {0} p7\tp7 name\tfrac index pair reads {0}\tnum reads {0}'.format(i), end='\t')
	for i in range(1, number_top_barcode_pairs + 1):
		print('top unknown barcode pair {0} p5\ttop unknown barcode pair {0} p7\tfrac index pair reads {0}\tnum reads {0}'.format(i), end='\t')
	print('')
	
	# first print all sample sheet entries, with the corresponding index pair kmer analysis
	sample_sheet_index_pairs = [parse_index_barcode_key(key)[0] for key in sample_sheet]
	# sort sample sheet by library id
	sorted_sample_sheet = sorted(sample_sheet, key=lambda x: sample_sheet[x].library_id)
	for index_barcode_key in sorted_sample_sheet:
		index_pair, barcode_pair = parse_index_barcode_key(index_barcode_key)
		output_kmer_line(index_pair, barcode_pair, sample_sheet[index_barcode_key])

	# now print any commonly occurring index pairs that have not been printed yet
	# counts and fractions of barcode pairs
	sorted_index_pair_counts = sorted(index_pair_counts, key=lambda x: index_pair_counts[x], reverse=True)
	for index_pair in sorted_index_pair_counts:
		if index_pair not in sample_sheet_index_pairs and index_pair_counts[index_pair] >= args.threshold_to_report:
			output_kmer_line(index_pair)
