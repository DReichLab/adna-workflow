# Process the results of angsd and return them in the keyed statistics format
# We use "new method 1" results

"Method1: new_llh Version: MoM:0.006946 SE(MoM):2.138220e-03 ML:0.006724 SE(ML):7.099542e-15"

import sys
import re

def parse_angsd_results(filename):
	nSNP_sites_pattern = re.compile('[\s]*We have nSNP sites:[\s]+([\d]+),')
	stats_search = "Method1: new_llh Version: "
	angsd = {
		"nsites": 0,
	}

	with open(filename) as f:
		for line in f:
			site_result = nSNP_sites_pattern.match(line)
			if site_result:
				angsd["nsites"] = int(site_result.group(1))
			if line.startswith(stats_search):
				key_value_pairs = line[len(stats_search):].split()
				for pair in key_value_pairs:
					key, value = pair.split(':')
					angsd[key] = float(value)
	return angsd

if __name__ == '__main__':

	angsd_filename = sys.argv[1]
	results = {}
	try:
		results = parse_angsd_results(angsd_filename)
	finally:
		print("angsd_{}\t{:d}".format("nsites", results.get("nsites", 0) ), end='\t')
		print("angsd_{}\t{:.3g}".format("MoM", results.get("MoM", -1.0) ), end='\t')
		print("angsd_{}\t{:.2g}".format("SE(MoM)", results.get("SE(MoM)", -1.0) ), end='\t')
		print("angsd_{}\t{:.3g}".format("ML", results.get("ML", -1.0) ), end='\t')
		print("angsd_{}\t{:.2g}".format("SE(ML)", results.get("SE(ML)", -1.0) ), end='\t')
		print('')
