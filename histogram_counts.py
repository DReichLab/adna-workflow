import argparse
import sys

from collections import Counter
from math import sqrt
from random import random, seed

def load_histogram(target_histogram_filename):
	counter = Counter()
	try:
		with open(target_histogram_filename) as histogram:
			for line in histogram:
				target_count, occurrences = map(lambda x: int(x), line.split())
				counter[target_count] = occurrences
	except:
		pass
	return counter

# given a read with n copies, and probability of retaining each copy, how many copies get retained
def retained_copies(copies, retain_probability):
	if retain_probability < 0 or retain_probability > 1.0:
		raise ValueError('bad retain_probability {:f}'.format(retain_probability))
	
	count = 0
	for i in range(copies):
		v = random()
		#print(v, file=sys.stderr)
		if v <= retain_probability:
			count += 1
	if count > copies:
		raise ValueError('{:d} is outside expected range of [0, {:d}]'.format(count, copies))
	return count

def total_and_unique_target_hits(counter, retain_probability=1.0):
	max_key = max(counter.keys())
	total_hits = 0
	unique_targets = 0
	for num_copies in counter:
		for i in range(counter[num_copies]):
			hits = retained_copies(num_copies, retain_probability)
			total_hits += hits
			if hits > 0:
				unique_targets += 1
	return total_hits, unique_targets

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="For analysis of histograms associated with preseq. From the histogram of the number of appearances of each read, report total hits and unique reads. Options relate to downsampling for computing preseq-like marginal uniquness computation.")
	parser.add_argument("-f", "--fraction", type=float, default=1.0, help='For downsampling, the fraction of reads that will be retained.')
	parser.add_argument("histogram_filename", help="preseq input histogram", nargs='+')
	parser.add_argument("-s", "--seed", help="random number seed for downsampling randomization", default=2020)
	parser.add_argument("-n", "--iterations", type=int, default=1, help='Number of iterations to average counts of reads and unique reads when downsampling')
	args = parser.parse_args()
	
	for histogram_filename in args.histogram_filename:
		seed(args.seed)
		counter = load_histogram(histogram_filename)
		total_hits = 0
		unique_targets = 0
		unique_targets_squared = 0
		for i in range(args.iterations):
			interation_hits, iteration_unique_targets = total_and_unique_target_hits(counter, args.fraction)
			total_hits += interation_hits
			unique_targets += iteration_unique_targets
			unique_targets_squared += iteration_unique_targets ** 2
		average_total_hits = total_hits / args.iterations
		average_unique_targets = unique_targets / args.iterations
		variance = unique_targets_squared / args.iterations - average_unique_targets ** 2
		standard_deviation = sqrt(variance)
		lower_interval = average_unique_targets - 1.96 * standard_deviation
		upper_interval = average_unique_targets + 1.96 * standard_deviation
		
		if args.iterations > 1:
			interval = '{:.3f}\t{:.3f}'.format(lower_interval, upper_interval)
			print(f'{average_total_hits}\t{average_unique_targets}\t{interval}')
		else:
			print(f'{average_total_hits}\t{average_unique_targets}')
