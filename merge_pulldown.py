import argparse
import filecmp
import fileinput
import shutil
from collections import OrderedDict

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
def merge_ind_files(output_filename, input_ind_filenames, max_overlap=1):
	ind_counts = [] # individual counts for each individual file
	individuals_geno_offsets_dict = OrderedDict() # mapping of geno sets to individuals by geno offset, treating geno data as one array
	VALID_SEX = frozenset('MFU')
	with open(output_filename, 'w') as out:
		overall_line_count = 0
		for ind_filename in input_ind_filenames:
			prior_files_total = 0 if len(ind_counts) == 0 else overall_line_count
			with open(ind_filename, 'r') as ind_file:
				for line in ind_file:
					fields = line.split()
					if len(fields) > 0:
						individual = fields[0]
						if individual not in individuals_geno_offsets_dict:
							out.write(line)
							individuals_geno_offsets_dict[individual] = [overall_line_count]
						else:
							# make sure there are no duplicate entries for an individual
							if len(individuals_geno_offsets_dict[individual]) >= max_overlap:
								raise ValueError('individual {} appears too many times'.format(individual))
							else:
								individuals_geno_offsets_dict[individual].append(overall_line_count)
						overall_line_count += 1
						# sex is one of M,F,U
						sex = fields[1]
						if sex not in VALID_SEX:
							raise ValueError('invalid sex {} in {}'.format(sex, line))
			
			ind_counts.append(overall_line_count - prior_files_total)
	individuals_geno_offsets = [individuals_geno_offsets_dict[x] for x in individuals_geno_offsets_dict]
	return ind_counts, individuals_geno_offsets

# First values are preferred, array must be non-empty
def merge_genotypes(array):
	geno = '9'
	VALID_GENO = frozenset('029')
	for value in array:
		if geno not in VALID_GENO:
			raise ValueError('invalid genotype value {}'.format(geno))
		geno = value
		if geno != '9':
			break
	return geno

def merge_geno_files(output_filename, input_geno_filenames, ind_counts, individuals_geno_offsets, num_snps):
	geno_files = []
	VALID_GENO = frozenset('029')
	try:
		# open files
		for geno_filename in input_geno_filenames:
			f = open(geno_filename, 'r')
			geno_files.append(f)
		with open(output_filename, 'w') as out:
			for i in range(num_snps):
				# read genotypes from all files into memory
				raw_genotypes_for_snp = ''
				for file_index in range(len(geno_files)):
					geno = geno_files[file_index].readline().strip()
					# check that there is one value for each individual
					if len(geno) != ind_counts[file_index]:
						raise ValueError('{}: expected {} genotypes, but found {}'.format(input_geno_filenames[file_index], ind_counts[file_index], len(geno)))
					# check geno values a
					values_appearing = set(geno)
					if not (values_appearing <= VALID_GENO):
						raise ValueError('invalid genotype value(s) {}'.format(values_appearing - VALID_GENO))
					raw_genotypes_for_snp += geno
				# merge genotypes
				merged_genotypes = []
				for offsets_for_individual in individuals_geno_offsets:
					individual_genos = [raw_genotypes_for_snp[offset] for offset in offsets_for_individual]
					individual_geno = merge_genotypes(individual_genos)
					merged_genotypes.append(individual_geno)
				out.write(''.join(str(g) for g in merged_genotypes))
				out.write('\n')
			# check that there is no leftover values
			for file_index in range(len(geno_files)):
				geno = geno_files[file_index].readline().strip()
				if len(geno) > 0:
					raise ValueError('extra values in {}'.format(input_geno_filenames[file_index]))
	finally:
		for f in geno_files:
			f.close()

# Perform snp file comparison (currently SNP files must be identical)
# Merge ind and genotype files
def merge_geno_snp_ind(geno_filenames, snp_filenames, ind_filenames, output_stem, max_overlap):
	num_snps = check_snp_file(snp_filenames[0])
	if not compare_snp_files(snp_filenames):
		raise ValueError('SNP file mismatch')
	
	shutil.copyfile(snp_filenames[0], "{}.snp".format(output_stem) )
	ind_counts, individuals_geno_offsets = merge_ind_files("{}.ind".format(output_stem), ind_filenames, max_overlap)
	merge_geno_files("{}.geno".format(output_stem), geno_filenames, ind_counts, individuals_geno_offsets, num_snps)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Merge pulldown eigenstrat results. SNP sets must be the same.")
	
	parser.add_argument('-i', "--input_stems", help="List of stems, with each stem having three files, for example: a0.{geno,snp,ind} or b0.b1.{geno,snp,ind}.", nargs='+', required=True)
	parser.add_argument('-o', "--output_stem", help="Stem for output files: output.{geno,snp,ind}.", default='out')
	parser.add_argument('-m', "--max_overlap", help="Number of results files an individual can appear in. For UDG half+minus, this should be 2. If an individual should only appear in one result set, this should be 1.", type=int, default=1)
	args = parser.parse_args()

	geno_filenames = ["{}.geno".format(stem) for stem in args.input_stems]
	snp_filenames = ["{}.snp".format(stem) for stem in args.input_stems]
	ind_filenames = ["{}.ind".format(stem) for stem in args.input_stems]
	output_stem = args.output_stem

	merge_geno_snp_ind(geno_filenames, snp_filenames, ind_filenames, output_stem, args.max_overlap)
