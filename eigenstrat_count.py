import argparse
import sys
from pathlib import Path

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Simple count of autosome SNPs in eigenstrat files")
	
	parser.add_argument("input", help="geno,ind,snp file root")
	args = parser.parse_args()

	geno_file = args.input + '.geno'
	ind_file = args.input + '.ind'
	snp_file = args.input + '.snp'

	with open(snp_file) as f:
		# assume chromosome sorted snp file
		autosome_snps = 0
		for line in f:
			fields = line.strip().split()
			try:
				chromosome = int(fields[1])
				if 1 <= chromosome and chromosome <= 22:
					autosome_snps += 1
			except:
				break
	print('{:d} autosome snps'.format(autosome_snps), file=sys.stderr)
	
	with open(ind_file) as f:
		ids = [line.split()[0] for line in f]
	print(len(ids), file=sys.stderr)
	
	counts = {}
	for instance_id in ids:
		counts[instance_id] = 0
	VALID_GENO = frozenset('029')
	DATA_SNP = frozenset('02')
	with open(geno_file) as f:
		line_number = 0
		for line in f:
			line_number += 1
			if line_number > autosome_snps:
				break
			line = line.strip()
			if len(ids) != len(line):
				raise ValueError('length mismatch on line {:d}: {:d} {:d}'.format(line_number, len(ids), len(line)))

			for index in range(len(ids)):
				if line[index] not in VALID_GENO:
					raise ValueError('bad genotype {} on line {:d}'.format(line[index], line_number))
				else:
					if line[index] in DATA_SNP:
						counts[ids[index]] += 1
	
	for instance_id in ids:
		print('{}\t{:d}'.format(instance_id, counts[instance_id]))
