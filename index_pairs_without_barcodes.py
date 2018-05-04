import argparse

class KeyedStatistics:
	def __init__(self):
		self.total_reads = 0
		self.samples = {}
	
	def read_file(self, filename):
		with open(filename) as f:
			# first line has read count
			line = f.readline()
			self.total_reads += int(line)
			
			for line in f:
				fields = line.split('\t')
				sample_key = fields[0]
				keys = fields[1::2]
				values = fields[2::2]
				
				keyed_values = {}
				for i in range(len(keys)):
					keyed_values[keys[i]] = int(values[i])
				self.samples[sample_key] = keyed_values
				

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Find index pairs that do not appear to have known barcodes. From barcode count statistics, generate an index-barcode file to use with demultiplexing to perform kmer analysis of reads that do not match known barcodes.")
	parser.add_argument("-t", "--threshold", help="threshold to demultiplex and examine unidentified barcodes", type=int, default=10000)
	parser.add_argument("statistics_file", help="keyed statistics file from barcode counting")
	args = parser.parse_args()
	
	statistics = KeyedStatistics()
	statistics.read_file(args.statistics_file)
	
	# filter to only index pairs that have significant reads
	filtered_statistics = [sample_key for sample_key in statistics.samples if statistics.samples[sample_key].get('without_barcodes', 0) >= args.threshold ]
	for sample_key in filtered_statistics:
		print("{}\t{}\t{}\t{}".format(sample_key, 'for', 'kmer analysis', 'only'))
