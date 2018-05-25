import argparse
import filecmp
import fileinput
import shutil

# return true if snp files are the same
# return false otherwise
def compare_snp_files(snp_file_array):
	for other_file in snp_file_array[1:]:
		if filecmp.cmp(snp_file_array[0], other_file, shallow=False) == False:
			return False
	return True

def check_snp_file(snp_filename):
	BASES = set('ACGT')
	count = 0
	with open(snp_filename, 'r') as snps:
		for snp_line in snps:
			fields = snp_line.split()
			if len(fields) != 6:
				raise ValueError("does not have expected 6 fields: '{}'".format(snp_line))
			chromosome = fields[1]
			if not chromosome.isdigit() or (int(chromosome) < 1 or int(chromosome) > 24):
				raise ValueError("bad chromosome '{}', {}".format(snp_line, chromosome))
			
			try:
				cm = float(fields[2])
			except:
				raise ValueError("nonfloat '{}': {}".format(snp_line, fields[2]))
			
			position = fields[3]
			if not position.isdigit():
				raise ValueError("bad position '{}': {}".format(snp_line, position))
			
			allele1 = fields[4]
			allele2 = fields[5]
			if (allele1 not in BASES) or (allele2 not in BASES) or (allele1 == allele2):
				raise ValueError("bad alleles '{}': {} {}".format(snp_line, allele1, allele2))
			count += 1
	return count

# combine the contents of multiple individual files into one
def write_ind_file(output_filename, input_ind_filenames):
	ind_counts = []
	individuals = set()
	VALID_SEX = frozenset('MFU')
	with open(output_filename, 'w') as out:
		for ind_filename in input_ind_filenames:
			count = 0
			with open(ind_filename, 'r') as ind_file:
				for line in ind_file:
					fields = line.split()
					if len(fields) > 0:
						out.write(line)
						count += 1
						# make sure there are no duplicate entries for an individual
						individual = fields[0]
						if individual in individuals:
							raise ValueError('individual {} appears more than once'.format(individual))
						# sex is one of M,F,U
						sex = fields[1]
						if sex not in VALID_SEX:
							raise ValueError('invalid sex {} in {}'.format(sex, line))
						individuals.add(individual)
			ind_counts.append(count)
	return ind_counts

def combine_geno_files(output_filename, input_geno_filenames, ind_counts, num_snps):
	geno_files = []
	VALID_GENO = frozenset('029')
	try:
		# open files
		for geno_filename in input_geno_filenames:
			f = open(geno_filename, 'r')
			geno_files.append(f)
		with open(output_filename, 'w') as out:
			for i in range(num_snps):
				for file_index in range(len(geno_files)):
					geno = geno_files[file_index].readline().strip()
					# check that there is one value for each individual
					if len(geno) != ind_counts[file_index]:
						raise ValueError('{}: expected {} genotypes, but found {}'.format(input_geno_filenames[file_index], ind_counts[file_index], len(geno)))
					# check geno values a
					values_appearing = set(geno)
					if not (values_appearing <= VALID_GENO):
						raise ValueError('invalid genotype value(s) {}'.format(values_appearing - VALID_GENO))
					# write geno values in order
					out.write(geno)
				out.write('\n')
			# check that there is no leftover values
			for file_index in range(len(geno_files)):
				geno = geno_files[file_index].readline().strip()
				if len(geno) > 0:
					raise ValueError('extra values in {}'.format(input_geno_filenames[file_index]))
	finally:
		for f in geno_files:
			f.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Merge pulldown eigenstrat results. SNP sets must be the same.")
	
	parser.add_argument('-i', "--input_stems", help="List of stems, with each stem having three files, for example: a0.{geno,snp,ind} or b0.b1.{geno,snp,ind}.", nargs='+', required=True)
	parser.add_argument('-o', "--output_stem", help="Stem for output files: output.{geno,snp,ind}.", default='out')
	args = parser.parse_args()

	geno_filenames = ["{}.geno".format(stem) for stem in args.input_stems]
	snp_filenames = ["{}.snp".format(stem) for stem in args.input_stems]
	ind_filenames = ["{}.ind".format(stem) for stem in args.input_stems]
	
	num_snps = check_snp_file(snp_filenames[0])
	if not compare_snp_files(snp_filenames):
		raise ValueError('SNP file mismatch')
	
	shutil.copyfile(snp_filenames[0], "{}.snp".format(args.output_stem) )
	ind_counts = write_ind_file("{}.ind".format(args.output_stem), ind_filenames)
	combine_geno_files("{}.geno".format(args.output_stem), geno_filenames, ind_counts, num_snps)
