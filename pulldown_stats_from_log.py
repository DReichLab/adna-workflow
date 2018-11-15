import argparse

# read a pulldown log file and return a map of ids to SNPs
def pulldown_snp_stats(filename):
	library_targets = {}
	with open(filename) as f:
		for line in f:
			if 'coverage' in line:
				fields = line.split()
				library_id = fields[0]
				targets = int(fields[-1])
				library_targets[library_id] = targets
	return library_targets
				

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Find number of 1240k targets hit by pulldown for samples")
	parser.add_argument('-l', "--log", help="pulldown log file to parse", required=True)
	args = parser.parse_args()
	
	library_targets = pulldown_snp_stats(args.log)
	print(library_targets)
