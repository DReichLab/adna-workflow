import argparse
import pysam
from Bio.Seq import Seq
from multiprocessing import Pool
from collections import Counter
import os

# If there are no barcodes, this may require LOTS of memory. For 6 base pair barcodes, there may be 4^12 ~= 16.8 million integers and labels, which is basically linear in number of reads
def most_common_barcode_pairs(filename, barcode_length, num_top_barcode_pairs):
	counts = Counter()
	bam = pysam.AlignmentFile(filename, "rb")
	for read in bam.fetch(until_eof=True):
		#print(read)
		barcode_p5 = read.query_sequence[0:barcode_length]
		barcode_p7 = Seq(read.query_sequence[-barcode_length:]).reverse_complement()
		#print(barcode_p5)
		#print(barcode_p7)
		key = "{}_{}".format(barcode_p5, barcode_p7)
		counts[key] += 1
	bam.close()
	top_barcode_pairs = counts.most_common(num_top_barcode_pairs)
	return filename, top_barcode_pairs

def output_result(args):
	fullpath = args[0]
	top_barcode_pairs = args[1]
	
	filename = os.path.splitext(os.path.basename(fullpath))[0]
	
	print(filename, end='\t')
	for (barcode_pair, count) in top_barcode_pairs:
		print("{}\t{:d}".format(barcode_pair, count), end='\t')
	print('')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Find the most common barcode pairs for a sam/bam file.")
	parser.add_argument("-b", "--barcode_length", help="keyed statistics file from barcode counting", type=int, default=6)
	parser.add_argument("-t", "--top_barcode_pairs", help="number of top barcode pairs to return", type=int, default=10)
	parser.add_argument("-p", "--pool_size", help="size of thread pool to process", type=int, default=4)
	parser.add_argument("bams", help="bam file(s) with potential barcodes inline", nargs='+')
	args = parser.parse_args()

	pool = Pool(processes=args.pool_size)
	results = [pool.apply_async(most_common_barcode_pairs, args=(bam, args.barcode_length, args.top_barcode_pairs), callback=output_result) for bam in args.bams]
	pool.close()
	pool.join()
	[x.get() for x in results] # check for failures
