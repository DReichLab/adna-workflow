import subprocess
import pathlib
import sys
import os
import argparse
from multiprocessing import Pool
import filecmp

from read_groups_from_bam import read_groups_from_bam

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
pmdlo: 3 
'''

def add_read_groups(adna_jar_filename, demultiplexed_bam_filename, output_bam_filename, bam_date_string, label, library_id, individual, working_directory):
	subprocess.run("java -Xmx5000m -jar {} AssignReadGroups -i {} -o {} -s {} -x {} -d {} -l {}".format(
		adna_jar_filename, 
		demultiplexed_bam_filename, 
		output_bam_filename,
		individual,
		label,
		bam_date_string,
		library_id
		)
	, shell=True, check=True, cwd=working_directory)

def build_release_library(adna_jar_filename, demultiplexed_bam_filenames, bam_date_strings, label, library_id, experiment, individual, picard_jar, working_directory):
	if len(demultiplexed_bam_filenames) != len(bam_date_strings):
		raise ValueError('Each bam needs a date string')
	# make a working directory for this library
	pathlib.Path(working_directory).mkdir(exist_ok=True)
	
	# add read groups for each library component
	count = 0
	library_component_bams = []
	for i in range(len(demultiplexed_bam_filenames)):
		count += 1
		demultiplexed_bam_filename = demultiplexed_bam_filenames[i]
		bam_date_string = bam_date_strings[i]
		output_bam_filename = "{0}_{1:d}.bam".format(library_id, count)
		add_read_groups(adna_jar_filename, demultiplexed_bam_filename, output_bam_filename, bam_date_string, label, library_id, individual, working_directory)
		library_component_bams.append(output_bam_filename)
	
	# use stderr by library
	with open('{}/stdout_build_release_library'.format(working_directory), 'w') as stdout_build, \
		open('{}/stderr_build_release_library'.format(working_directory), 'w') as stderr_build:
		# merge 
		library_with_duplicates_filename = "{0}.duplicates.bam".format(library_id)
		subprocess.run("java -Xmx5500m -jar {} MergeSamFiles I={} O={} SORT_ORDER=coordinate".format(picard_jar, ' I='.join(library_component_bams), library_with_duplicates_filename), shell=True, check=True, cwd=working_directory, stdout=stdout_build, stderr=stderr_build)
		
		# deduplicate
		library_filename = "{0}.{}.bam".format(library_id, experiment)
		subprocess.run("java -Xmx5500m -jar {0} MarkDuplicates I={1} O={2} M={2}.dedup_stats REMOVE_DUPLICATES=true BARCODE_TAG=XD ADD_PG_TAG_TO_READS=false MAX_FILE_HANDLES=1000".format(picard_jar, library_with_duplicates_filename, library_filename), shell=True, check=True, cwd=working_directory, stdout=stdout_build, stderr=stderr_build)
		
	# discover read groups
	read_groups = read_groups_from_bam("{}/{}".format(working_directory, library_filename))
		
	return library_filename, read_groups

def copy_release_library(library_parameters, destination_parent_directory):
	# create library directory if it does not exist
	destination_library_path = library_parameters.get_release_library_path(destination_parent_directory)
	pathlib.Path(destination_library_path).mkdir(exist_ok=True)
	shutil.copy(library_parameters.release_library_name, destination_library_path)

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
				individual_file_damage_restricted.write("{0}\t{1}\t{0}".format(library_id, sex))
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
			library_bam_filename  = "{}/{}/{}".format(args.release_directory, parameters.library_id, parameters.release_library_name)
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
'''
	# verify that snp files are identical
	snp_file_check = filecmp.cmp("{}/{}.snp".format(working_directory, library_id), 
							  "{}/{}_d.snp".format(working_directory, library_id), shallow=False)
	if not snp_file_check:
		raise ValueError('snp files are different for normal and damage-restricted pulldown')
'''

def index_library(library_parameters, release_directory):
	release_library_path = library_parameters.get_release_library_path(destination_parent_directory)
	library_bam_filename = library_parameters.release_library_name
	subprocess.run(['samtools', 'index', library_bam_filename], check=True, cwd=release_library_path)
	
class LibraryParameters:
	def __init__(self, line):
		parameters_for_library = line.split('\t')
		self.index_barcode_key = parameters_for_library[0]
		self.library_id = parameters_for_library[1]
		self.individual_id = parameters_for_library[2]
		self.read_group_description = parameters_for_library[3]
		self.experiment = parameters_for_library[4]
		self.udg = parameters_for_library[5]
		self.bam_filenames = parameters_for_library[6::2]
		self.bam_date_strings = [date_string.strip() for date_string in parameters_for_library[7::2]]
	
	def set_release_library_name(self, release_library_name):
		self.release_library_name = release_library_name
		
	def get_release_library_path(self, release_directory):
		release_library_path = "{}/{}".format(release_directory, self.library_id)
		return release_library_path
			
# each line contains the parameters to build one library for pulldown
def build_library_and_pulldown(line, args):
	adna_jar_filename = args.adna_jar
	picard_jar = args.picard_jar
	report_filename = args.report
	
	if args.pulldown is not None:
		pulldown_executable = args.pulldown[0]
		
	
	
	working_directory = library_id
	
	release_library_filename = build_release_library(adna_jar_filename, bam_filenames, bam_date_strings, read_group_description, library_id, individual_id, experiment, picard_jar, working_directory)
	
	if args.pulldown is not None:
		pulldown(release_library_filename, library_id, report_filename, experiment, pulldown_executable, working_directory, index_barcode_key)
	# TODO copy pulldown results
	# TODO check SNP file to see whether it needs to be copied
	return release_library_filename

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Combine ancient DNA analysis outputs into a single file keyed by index-barcode key. Merges component bams and deduplicates library. Runs pulldown.")
	
	# pulldown is optional
	# 1240k(+) libraries require pulldown. MT libraries do not. 
	parser.add_argument('-p', "--pulldown_binary", help="binary to run pulldown", nargs=1)
	parser.add_argument('-l', "--pulldown_label", help="label for pulldown filenames", nargs=1)
	parser.add_argument('-n', "--num_threads", help="size of thread pool", type=int, default =10, nargs=1)
	parser.add_argument('-r', "--release_directory", help="parent directory to put released libraries", nargs=1)
	
	parser.add_argument("bam_list", help="Each line contains the parameters to build a library bam for release. This includes the library ID, the individual ID, experiment, read group description (sequencing run name with experiment type and udg treatment), experiment, and (bam, sequencing run date) pairs ")
	parser.add_argument("adna_jar", help="jar file for assigning read groups to each read in a bam")
	parser.add_argument("picard_jar", help="jar file the Broad Institute Picard toolset, used to merge bams and deduplicate")
	parser.add_argument("report", help="report is used for looking up sex")
	args = parser.parse_args()
	
	bam_list_filename = args.bam_list
	# read list of files
	with open(bam_list_filename) as f:
		library_parameters = [LibraryParameters(line) for line in f]
	# build libraries
	adna_jar_filename = args.adna_jar
	picard_jar = args.picard_jar
	
	pool = Pool(processes=args.num_threads)
	for parameters in library_parameters:
		working_directory = library_id
		parameters.build_result = pool.apply_async(build_release_library, args=(adna_jar_filename, parameters.bam_filenames, parameters.bam_date_strings, parameters.read_group_description, parameters.library_id, parameters.individual_id, parameters.experiment, picard_jar, working_directory))
	pool.close()
	pool.join()
	
	for parameters in library_parameters:
		library_filename, read_groups = parameters.build_result.get()
		parameters.set_release_library_name(library_filename)
		parameters.read_groups = read_groups
	# copy libraries
	pool = Pool(processes=args.num_threads)
	if args.release_directory is not None:
		[pool.apply_async(copy_release_library, args=(parameters, args.release_directory)) for parameters in library_parameters]
	pool.close()
	pool.join()
	
	# index released libraries
	pool = Pool(processes=args.num_threads)
	for parameters in library_parameters:
		pool.apply_async(index_library, args=(parameters, args.release_directory))
	pool.close()
	pool.join()
	
	# prepare pulldown input files
	if args.pulldown_label is not None:
		prepare_pulldown(library_parameters, args)
