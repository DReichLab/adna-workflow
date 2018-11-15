import argparse

# read a pulldown log file and return a map of ids to SNPs
def pulldown_snp_stats(filename):
	targets_by_instance = {}
	with open(filename) as f:
		for line in f:
			if 'coverage' in line:
				fields = line.split()
				instance_id = fields[0]
				targets = int(fields[-1])
				targets_by_instance[instance_id] = targets
	return targets_by_instance
				

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Find number of 1240k targets hit by pulldown for samples")
	parser.add_argument('-l', "--log", help="pulldown log file to parse", required=True)
	args = parser.parse_args()
	
	targets_by_instance = pulldown_snp_stats(args.log)
	for instance, targets in library_targets.items():
		print("{}\tpulldown_coverage\t{:d}".format(instance, targets))
