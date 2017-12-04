from __future__ import print_function
import sys

# Analyze the distribution of barcode pairs for each index pair
# Each index pair is output on a row, with the 10 most common barcode pairs in descending order of counts

barcodes_filename = sys.argv[1]
counts_filename = sys.argv[2]

def parse_index_barcode_key(index_barcode_key):
	i5, i7, p5, p7 = index_barcode_key.split('_')
	index_pair = i5 + '_' + i7
	barcode_pair = p5 + '_' + p7
	return index_pair, barcode_pair

index_pairs = dict()
index_pair_counts = dict()
barcodes = dict()

with open(barcodes_filename) as f:
	for line in f:
		fields = line.strip().split('\t')
		barcode_sequence = fields[0]
		name = fields[1]
		barcodes[barcode_sequence] = name

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
			index_pairs[index_pair] = dict()
			
		# accounting for barcode pair within index pair
		index_pair_counts[index_pair] = index_pair_counts.get(index_pair, 0) + int(raw_count)
		index_pairs[index_pair][barcode_pair] = int(raw_count)
		
number_top_barcode_pairs = 10
# print field headers
print ('i5 index\ti7 index\tnum reads\tfrac total\t', end='')
for i in range(1, number_top_barcode_pairs + 1):
	print('top barcode pair {0} p5\tp5 name\ttop barcode pair {0} p7\tp7 name\tfrac index pair reads {0}\tnum reads {0}'.format(i) ,end='\t')
print('')

# counts and fractions of barcode pairs
sorted_index_pair_counts = sorted(index_pair_counts, key=lambda x: index_pair_counts[x], reverse=True)
for index_pair in sorted_index_pair_counts:
	index_pair_count = index_pair_counts[index_pair]
	barcode_pair_counts = index_pairs[index_pair]
	sorted_barcode_pairs = sorted(barcode_pair_counts, key=lambda x: barcode_pair_counts[x], reverse=True)
	i5, i7 = index_pair.split('_')
	print('{}\t{}\t{:d}\t{:.4f}'.format(i5, i7, index_pair_count, float(index_pair_count) / total_reads), end='')
	for barcode_pair in sorted_barcode_pairs[:number_top_barcode_pairs]:
		barcode_pair_count = barcode_pair_counts[barcode_pair]
		p5, p7 = barcode_pair.split('_')
		print ('\t{}\t{}\t{}\t{}\t{:.3f}\t{:d}'.format(p5, barcodes.get(p5, ''), p7, barcodes.get(p7, ''), float(barcode_pair_count) / index_pair_count, barcode_pair_count), end='')
	print('')
	
