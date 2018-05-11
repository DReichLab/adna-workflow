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
	with open(filename, 'w') as individual_file:
		for parameters in library_parameters:
			if parameters.udg == udg_filter:
				library_id = id_transform_function(parameters.library_id)
				sex = sex_by_index_barcode_key[parameters.index_barcode_key]
				individual_file.write("{0}\t{1}\t{0}".format(library_id, sex))
	return filename

pulldown_options = {'normal' : non_damage_restricted_options, 'damage_restricted' : damage_restricted_options}
def create_pulldown_parameter_file(pulldown_label, pulldown_type):
	pulldown_string = "{}.{}".format(pulldown_label, pulldown_type)
	pulldown_parameters = pulldown_parameters_base.replace('LABEL', pulldown_label).replace('PULLDOWN_STRING', pulldown_string) + pulldown_options[pulldown_type]
	pulldown_parameter_filename_nopath = "{}.parameters".format(pulldown_string)
	with open("{}".format(pulldown_parameter_filename_nopath), 'w') as pulldown_parameter_file:
		pulldown_parameter_file.write(pulldown_parameters)

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
	
	# Generate individual files
	# normal pulldown for udg half individual file
	create_individual_file("{}.half.normal.ind".format(pulldown_label), library_parameters, HALF, identity, sex_by_index_barcode_key)
	# damage restricted pulldown for udg half individual file
	create_individual_file("{}.half.damage_restricted.ind".format(pulldown_label), library_parameters, HALF, damage_restricted_d, sex_by_index_barcode_key)
	
	# normal pulldown for udg minus individual file
	create_individual_file("{}.minus.normal.ind".format(pulldown_label), library_parameters, MINUS, identity, sex_by_index_barcode_key)
	# damage restricted pulldown for udg minus individual file
	create_individual_file("{}.minus.damage_restricted.ind".format(pulldown_label), library_parameters, MINUS, damage_restricted_d, sex_by_index_barcode_key)
	
	# UDG plus cannot be damage restricted
	create_individual_file("{}.plus.normal.ind".format(pulldown_label), library_parameters, PLUS, identity, sex_by_index_barcode_key)
	
	for parameters in library_parameters:
		if parameters.udg not in ALLOWED_UDG_VALUES:
			raise ValueError('Unhandled udg {}'.format(parameters.udg))
			
	# build Nick-style database file
	current_directory = os.getcwd()
	with open("{}.dblist".format(pulldown_label), 'w') as db_file:
		for parameters in library_parameters:
			# filename must match the deduplicated bam
			release_library_name = library_parameters.get_release_library_name()
			
			library_bam_filename  = "{}/{}/{}".format(args.release_directory, parameters.library_id, release_library_name)
			db_line = "{0}\t{0}\t{1}\t{2}".format(library_id, library_bam_filename, ":".join(read_groups))
			db_file.write(db_line)
			
	# build parameter files
	pulldown_options = {'normal' : non_damage_restricted_options, 'damage_restricted' : damage_restricted_options}
	for pulldown_type in pulldown_options:
		pulldown_string = "{}.{}".format(pulldown_label, pulldown_type)
		pulldown_parameters = pulldown_parameters_base.replace('LABEL', pulldown_label).replace('PULLDOWN_STRING', pulldown_string) + pulldown_options[pulldown_type]
		pulldown_parameter_filename_nopath = "{}.parameters".format(pulldown_string)
		with open("{}".format(pulldown_parameter_filename_nopath), 'w') as pulldown_parameter_file:
			pulldown_parameter_file.write(pulldown_parameters)
		
	# damage restricted
	damage_restricted_parameters = pulldown_parameters_library.replace('SAMPLE_MACRO', library_id + '_d') + damage_restricted_options
	pulldown_parameter_damage_restricted_filename_nopath = "{}.damage_restricted.parameters".format(library_id)
	with open("{}/{}".format(working_directory, pulldown_parameter_damage_restricted_filename_nopath), 'w') as pulldown_parameter_file:
		pulldown_parameter_file.write(damage_restricted_parameters)

def pulldown(pulldown_parameter_filename):
	with open('{}.stdout'.format(pulldown_parameter_filename), 'w') as stdout_pulldown, \
		open('{}.stderr'.format(pulldown_parameter_filename), 'w') as stderr_pulldown:
			subprocess.run([pulldown_executable, "-p", pulldown_parameter_filename], check=True, cwd=working_directory, stdout=stdout_pulldown, stderr=stderr_pulldown)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Prepare pulldown input files for a batch.")
	
	# pulldown is optional
	# 1240k(+) libraries require pulldown. MT libraries do not. 
	parser.add_argument('-p', "--pulldown_binary", help="binary to run pulldown")
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
		prepare_pulldown(library_parameters, args)
