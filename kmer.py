from __future__ import print_function
import sys

filename = sys.argv[1]

def parse_index_barcode_key(index_barcode_key):
	i5, i7, p5, p7 = index_barcode_key.split('_')
	index_pair = i5 + '_' + i7
	barcode_pair = p5 + '_' + p7
	return index_pair, barcode_pair

index_pairs = dict()
index_pair_counts = dict()

with open(filename) as f:
	line = f.readline() # total reads
	totalReads = int(line)
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

sorted_index_pair_counts = sorted(index_pair_counts, key=lambda x: index_pair_counts[x], reverse=True)
for index_pair in sorted_index_pair_counts:
	index_pair_count = index_pair_counts[index_pair]
	barcode_pair_counts = index_pairs[index_pair]
	sorted_barcode_pairs = sorted(barcode_pair_counts, key=lambda x: barcode_pair_counts[x], reverse=True)
	print('%s\t%d' % (index_pair, index_pair_count), end='')
	for barcode_pair in sorted_barcode_pairs[:10]:
		barcode_pair_count = barcode_pair_counts[barcode_pair]
		print ('\t%s\t%.3f\t%d' % (barcode_pair, float(barcode_pair_count) / index_pair_count, barcode_pair_count), end='')
	print('')
	
