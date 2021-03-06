import subprocess
import pathlib
import sys
import os
import argparse
from multiprocessing import Pool
import filecmp

from read_groups_from_bam import read_groups_from_bam
from release_libraries import LibraryParameters
from merge_pulldown import merge_geno_snp_ind

SNP_SETS = {
	'1240k' : '1240kSNP.snp',
	'BigYoruba' : 'bigYRI_v3__excluding__1240k_v3.snp',
	'BigYoruba+1240k' : 'bigyoruba_1240k.snp'
}

# Nick's pulldown program has an input file 
pulldown_parameters_base = '''BASE: /n/groups/reich/matt/pipeline/static
indivname:      PULLDOWN_STRING.ind
snpname:        BASE/SNP_SET
indivoutname:     PULLDOWN_STRING.out.prelim.ind   
snpoutname:       PULLDOWN_STRING.out.snp
genotypeoutname:  PULLDOWN_STRING.out.geno
outputformat:     eigenstrat
threshtable:      BASE/pulldown_thresholds 
defstring:        capture_UDG
dbbam:            LABEL.dblist
readbam:          BASE/readbam
majmode:          NO
printcount:       NO
udgmode: UDG
minlen: 30
maxlen: 123
oldpullmode:	YES
'''

non_damage_restricted_options = 'pmdscore:         NO'
# For damage pulldown
damage_restricted_options = '''pmdscore: YES
pmdlo: 2 
'''

class PulldownResults:
	def __init__(self, ind, snp, geno):
		self.ind = ind
		self.snp = snp
		self.geno = geno

# build map for sex from report
def sex_from_report(report_filename):
	sex_by_index_barcode_key = {}
	# consider doing this from the index-barcode key instead of library ID and experiment
	with open(report_filename) as report:
		report.readline() # ignore read count line
		header_line = report.readline()
		headers = header_line.split('\t')
		index_barcode_key_index = headers.index('Index-Barcode Key')
		sex_index = headers.index('1240k_post_sex')
		sex = 'U'
		for line in report:
			fields = line.split('\t')
			try:
				index_barcode_key_field = fields[index_barcode_key_index]
				if len(fields[sex_index]) > 0:
					sex = fields[sex_index]
				sex_by_index_barcode_key[index_barcode_key_field] = sex
			except:
				pass
	return sex_by_index_barcode_key

def create_individual_file(filename, library_parameters, udg_filter, sex_by_index_barcode_key, exclude_library_list):
	count = 0
	with open(filename, 'w') as individual_file:
		for parameters in library_parameters:
			if parameters.udg == udg_filter and parameters.library_id not in exclude_library_list:
				count += 1
				library_id = parameters.library_id
				sex = sex_by_index_barcode_key[parameters.index_barcode_key]
				individual_file.write("{0}\t{1}\t{0}\n".format(library_id, sex))
	return count

# To reuse dblist for damage-restricted pulldown, we modify the individual file paired with eigenstrat output 
# after pulldown to create unique identifiers
def rewrite_individual_file(input_filename, output_filename, to_append):
	with open(input_filename, 'r') as original_file:
		with open(output_filename, 'w') as modified_file:
			for line in original_file:
				fields = line.rsplit(None, 2) # preserve whitespace at left
				fields[0] = fields[0] + to_append # add modifier
				modified_file.write(" ".join(fields) + "\n")

NORMAL = 'normal'
DAMAGE_RESTRICTED = 'damage_restricted'
pulldown_damage_options = {NORMAL : non_damage_restricted_options, DAMAGE_RESTRICTED : damage_restricted_options}
def create_pulldown_parameter_file(pulldown_label, pulldown_base_filename, udg_type, damage_type, snp_set_filename):
	pulldown_parameters = pulldown_parameters_base.replace('LABEL', pulldown_label).replace('SNP_SET', snp_set_filename).replace('UDG', udg_type).replace('PULLDOWN_STRING', pulldown_base_filename) + pulldown_damage_options[damage_type]
	pulldown_parameter_filename_nopath = "{}.parameters".format(pulldown_base_filename)
	with open("{}".format(pulldown_parameter_filename_nopath), 'w') as pulldown_parameter_file:
		pulldown_parameter_file.write(pulldown_parameters)
	return pulldown_parameter_filename_nopath

def prepare_pulldown(library_parameters, args):
	# input files for combined pulldown
	sex_by_index_barcode_key = sex_from_report(args.report)
	pulldown_label = args.pulldown_label
	
	# Different UDG treatments determine pulldown parameters
	# Separate pulldown by udg type
	MINUS = 'minus'
	HALF = 'half'
	PLUS = 'plus'
	ALLOWED_UDG_VALUES = [MINUS, HALF, PLUS]
	
	# check that all UDG values are allowed
	for parameters in library_parameters:
		if parameters.udg not in ALLOWED_UDG_VALUES:
			raise ValueError('Unhandled udg {}'.format(parameters.udg))
		
	# build Nick-style database file
	exclude_library_list = []
	current_directory = os.getcwd()
	with open("{}.dblist".format(pulldown_label), 'w') as db_file:
		for parameters in library_parameters:
			# filename must match the deduplicated bam
			release_library_name = parameters.get_release_library_name()
			
			library_bam_filename  = "{}/{}/{}".format(args.release_directory, parameters.library_id, release_library_name)
			read_groups = read_groups_from_bam(library_bam_filename)
			if len(read_groups) > 0:
				db_line = "{0}\t{0}\t{1}\t{2}\n".format(parameters.library_id, library_bam_filename, ":".join(read_groups))
				db_file.write(db_line)
			else:
				exclude_library_list.append(parameters.library_id)
	
	parameter_indices = [
		(HALF, NORMAL),
		(HALF, DAMAGE_RESTRICTED),
		(MINUS, NORMAL),
		(MINUS, DAMAGE_RESTRICTED),
		(PLUS, NORMAL) # UDG plus cannot be damage restricted
		]
	parameter_file_outputs = []
	for (udg_type, damage_type) in parameter_indices:
		pulldown_base_filename = "{}.{}.{}".format(pulldown_label, udg_type, damage_type)
		# Generate individual files
		count = create_individual_file(pulldown_base_filename + ".ind", library_parameters, udg_type, sex_by_index_barcode_key, exclude_library_list)
		# build parameter files
		snp_set_filename = SNP_SETS[args.snp_set]
		pulldown_parameter_filename_nopath = create_pulldown_parameter_file(pulldown_label, pulldown_base_filename, udg_type, damage_type, snp_set_filename)
		if count > 0:
			parameter_file_outputs.append(pulldown_parameter_filename_nopath)
	
	return parameter_file_outputs

def read_from_parameter_file(parameter_filename):
	individual_filename = None
	snpout_filename = None
	geno_filename = None
	damage_restricted = False
	with open(parameter_filename) as f:
		for line in f:
			key, value = line.split()
			if key == "indivoutname:":
				individual_filename = value
			elif key == "pmdscore:" and value == 'YES':
				damage_restricted = True
			elif key == "snpoutname:":
				snpout_filename = value
			elif key == "genotypeoutname:":
				geno_filename = value
	return individual_filename, snpout_filename, geno_filename, damage_restricted

def pulldown(pulldown_parameter_filename, pulldown_executable):
	with open('{}.stdout'.format(pulldown_parameter_filename), 'w') as stdout_pulldown, \
		open('{}.stderr'.format(pulldown_parameter_filename), 'w') as stderr_pulldown:
			subprocess.run([pulldown_executable, "-p", pulldown_parameter_filename], check=True, stdout=stdout_pulldown, stderr=stderr_pulldown)
	# modify individual file, if damage restricted. Otherwise simply copy
	individual_filename, snpout_filename, geno_filename, damage_restricted = read_from_parameter_file(pulldown_parameter_filename)
	modified_individual_filename = individual_filename.replace('.prelim', '')
	modifier = '_d' if damage_restricted else ''
	rewrite_individual_file(individual_filename, modified_individual_filename, modifier)
	
	return PulldownResults(modified_individual_filename, snpout_filename, geno_filename)

def merge_pulldowns(pulldown_label, pulldown_file_sets, max_overlap):
	geno_files = [pulldown_results.geno for pulldown_results in pulldown_file_sets]
	snp_files = [pulldown_results.snp for pulldown_results in pulldown_file_sets]
	ind_files = [pulldown_results.ind for pulldown_results in pulldown_file_sets]
	
	output_stem = "{}.combined".format(pulldown_label)
	merge_geno_snp_ind(geno_files, snp_files, ind_files, output_stem, max_overlap)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Prepare pulldown input files for a list of libraries and run pulldown.")
	
	# pulldown is optional
	# 1240k(+) libraries require pulldown. MT libraries do not. 
	parser.add_argument('-p', "--pulldown_executable", help="executable to run pulldown", required=True)
	parser.add_argument('-l', "--pulldown_label", help="label for pulldown filenames", required=True)
	parser.add_argument('-r', "--release_directory", help="parent directory to read released libraries", required=True)
	
	parser.add_argument('--snp_set', choices=['1240k', 'BigYoruba', 'BigYoruba+1240k'], help="SNP set to use for pulldown", default='1240k')
	# TODO hook up these arguments correctly
	parser.add_argument("--minimum_length", help="minimum length of read to include in pulldown", type=int, default=30)
	parser.add_argument("--maximum_length", help="maximum length of read to include in pulldown", type=int, default=123)
	
	parser.add_argument("bam_list", help="Each line contains the parameters to build a library bam for release. This includes the library ID, the individual ID, experiment, read group description (sequencing run name with experiment type and udg treatment), experiment, and (bam, sequencing run date) pairs ")
	parser.add_argument("report", help="report is used for looking up sex")
	args = parser.parse_args()
	
	bam_list_filename = args.bam_list
	# read list of files
	with open(bam_list_filename) as f:
		unfiltered_library_parameters = [LibraryParameters(line) for line in f]
		library_parameters = [par for par in unfiltered_library_parameters if (args.snp_set in par.experiment)]
	
	# prepare pulldown input files
	if args.pulldown_label is not None:
		parameter_files = prepare_pulldown(library_parameters, args)
	
		# run pulldown
		pool = Pool(processes=2)
		results = [pool.apply_async(pulldown, args=(parameter_file, args.pulldown_executable)) for parameter_file in parameter_files]
		pool.close()
		pulldown_file_sets = [result.get() for result in results] # check for exceptions
		pool.join()
		
		# merge pulldown results
		merge_pulldowns(args.pulldown_label, pulldown_file_sets, 1) # libraries have only 1 UDG treatment
