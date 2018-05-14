import subprocess
import pathlib
import sys
import os
import argparse
from multiprocessing import Pool
import filecmp

from read_groups_from_bam import read_groups_from_bam
from release_libraries import LibraryParameters

# Nick's pulldown program has an input file 
pulldown_parameters_base = '''BASE: /home/np29
TT:  BASE/tables
BB:  BASE/o2bin
DIR_MACRO:           .
indivname:      DIR_MACRO/PULLDOWN_STRING.ind
snpname:        /n/groups/reich/matt/pipeline/static/1240kSNP.snp
indivoutname:     DIR_MACRO/PULLDOWN_STRING.out.ind   
snpoutname:       DIR_MACRO/PULLDOWN_STRING.out.snp
genotypeoutname:  DIR_MACRO/PULLDOWN_STRING.out.geno
outputformat:     eigenstrat
threshtable:      TT/defaultthresh 
defstring:        capture
dbbam:            DIR_MACRO/LABEL.dblist
readbam:          BB/readbam
majmode:          NO
printcount:       NO
udgmode: UDG
'''

non_damage_restricted_options = 'pmdscore:         NO'
# For damage pulldown
damage_restricted_options = '''pmdscore: YES
pmdlo: 2 
'''

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
				sex = fields[sex_index]
				sex_by_index_barcode_key[index_barcode_key_field] = sex
			except:
				pass
	return sex_by_index_barcode_key

def identity(string):
	return string

def damage_restricted_d(string):
	return "{}_d".format(string)

def create_individual_file(filename, library_parameters, udg_filter, id_transform_function, sex_by_index_barcode_key):
	count = 0
	with open(filename, 'w') as individual_file:
		for parameters in library_parameters:
			if parameters.udg == udg_filter:
				count += 1
				library_id = id_transform_function(parameters.library_id)
				sex = sex_by_index_barcode_key[parameters.index_barcode_key]
				individual_file.write("{0}\t{1}\t{0}\n".format(library_id, sex))
	return count

NORMAL = 'normal'
DAMAGE_RESTRICTED = 'damage_restricted'
pulldown_damage_options = {NORMAL : non_damage_restricted_options, DAMAGE_RESTRICTED : damage_restricted_options}
def create_pulldown_parameter_file(pulldown_label, pulldown_base_filename, udg_type, damage_type):
	pulldown_parameters = pulldown_parameters_base.replace('LABEL', pulldown_label).replace('UDG', udg_type).replace('PULLDOWN_STRING', pulldown_base_filename) + pulldown_damage_options[damage_type]
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
		name_manipulation_function = identity if damage_type == NORMAL else damage_restricted_d
		count = create_individual_file(pulldown_base_filename + ".ind", library_parameters, udg_type, name_manipulation_function, sex_by_index_barcode_key)
		# build parameter files
		pulldown_parameter_filename_nopath = create_pulldown_parameter_file(pulldown_label, pulldown_base_filename, udg_type, damage_type)
		if count > 0:
			parameter_file_outputs.append(pulldown_parameter_filename_nopath)
		
	# build Nick-style database file
	current_directory = os.getcwd()
	with open("{}.dblist".format(pulldown_label), 'w') as db_file:
		for parameters in library_parameters:
			# filename must match the deduplicated bam
			release_library_name = parameters.get_release_library_name()
			
			library_bam_filename  = "{}/{}/{}".format(args.release_directory, parameters.library_id, release_library_name)
			read_groups = read_groups_from_bam(library_bam_filename)
			db_line = "{0}\t{0}\t{1}\t{2}\n".format(parameters.library_id, library_bam_filename, ":".join(read_groups))
			db_file.write(db_line)
	
	return parameter_file_outputs

def pulldown(pulldown_parameter_filename, pulldown_executable):
	with open('{}.stdout'.format(pulldown_parameter_filename), 'w') as stdout_pulldown, \
		open('{}.stderr'.format(pulldown_parameter_filename), 'w') as stderr_pulldown:
			subprocess.run([pulldown_executable, "-p", pulldown_parameter_filename], check=True, stdout=stdout_pulldown, stderr=stderr_pulldown)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Prepare pulldown input files for a batch.")
	
	# pulldown is optional
	# 1240k(+) libraries require pulldown. MT libraries do not. 
	parser.add_argument('-p', "--pulldown_executable", help="executable to run pulldown")
	parser.add_argument('-l', "--pulldown_label", help="label for pulldown filenames")
	#parser.add_argument('-n', "--num_threads", help="size of thread pool", type=int, default =10)
	parser.add_argument('-r', "--release_directory", help="parent directory to put released libraries")
	
	parser.add_argument("bam_list", help="Each line contains the parameters to build a library bam for release. This includes the library ID, the individual ID, experiment, read group description (sequencing run name with experiment type and udg treatment), experiment, and (bam, sequencing run date) pairs ")
	parser.add_argument("report", help="report is used for looking up sex")
	args = parser.parse_args()
	
	bam_list_filename = args.bam_list
	# read list of files
	with open(bam_list_filename) as f:
		library_parameters = [LibraryParameters(line) for line in f]
	
	# prepare pulldown input files
	if args.pulldown_label is not None:
		parameter_files = prepare_pulldown(library_parameters, args)
	
	pool = Pool(processes=2)
	results = [pool.apply_async(pulldown, args=(parameter_file, args.pulldown_executable)) for parameter_file in parameter_files]
	pool.close()
	[result.get() for result in results] # check for exceptions
	pool.join()
