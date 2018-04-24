import subprocess
import pathlib
import sys
import os
import argparse
from multiprocessing import Pool

from read_groups_from_bam import read_groups_from_bam

# Nick's pulldown program has an input file 
pulldown_parameters_base = '''BASE: /home/np29
TT:  BASE/tables
BB:  BASE/o2bin
DIR_MACRO:           .
indivname:      DIR_MACRO/SAMPLE_EXTRACT_LIBRARY_ID.ind
snpname:        /n/groups/reich/matt/pipeline/static/1240kSNP.snp
indivoutname:     DIR_MACRO/SAMPLE_MACRO.ind   
snpoutname:       DIR_MACRO/SAMPLE_MACRO.snp
genotypeoutname:  DIR_MACRO/SAMPLE_MACRO.geno
outputformat:     eigenstrat
threshtable:      TT/defaultthresh 
defstring:        capture
dbbam:            DIR_MACRO/SAMPLE_EXTRACT_LIBRARY_ID.dblist
readbam:          BB/readbam
majmode:          NO
printcount:       NO
udgmode: half
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
	
	# merge 
	library_with_duplicates_filename = "{0}.duplicates.bam".format(library_id)
	subprocess.run("java -Xmx5500m -jar {} MergeSamFiles I={} O={} SORT_ORDER=coordinate".format(picard_jar, ' I='.join(library_component_bams), library_with_duplicates_filename), shell=True, check=True, cwd=working_directory)
	
	# deduplicate
	library_filename = "{0}.bam".format(library_id)
	subprocess.run("java -Xmx5500m -jar {0} MarkDuplicates I={1} O={2} M={2}.dedup_stats REMOVE_DUPLICATES=true BARCODE_TAG=XD ADD_PG_TAG_TO_READS=false MAX_FILE_HANDLES=1000".format(picard_jar, library_with_duplicates_filename, library_filename), shell=True, check=True, cwd=working_directory)
	return library_filename

def pulldown(library_bam_filename, library_id, report_filename, experiment, pulldown_executable, working_directory, index_barcode_key):
	# index bam
	subprocess.run(['samtools', 'index', library_bam_filename], check=True, cwd=working_directory)
	
	# discover read groups
	read_groups = read_groups_from_bam("{}/{}".format(working_directory, library_bam_filename))
	
	# build Nick-style database file
	current_directory = os.getcwd()
	db_line = "{}\t{}\t{}\t{}".format(working_directory, library_id, library_bam_filename, ":".join(read_groups))
	with open("{}/{}.dblist".format(working_directory, library_id), 'w') as db_file:
		db_file.write(db_line)
		
	# build individual file
	# lookup sex from report
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
				if index_barcode_key == index_barcode_key_field:
					sex = fields[sex_index]
					break
			except:
				pass
	with open("{}/{}.ind".format(working_directory, library_id), 'w') as individual_file:
		individual_file.write("{0}\t{1}\t{0}".format(library_id, sex))
			
	# build parameter files
	pulldown_parameters_library = pulldown_parameters_base.replace('SAMPLE_EXTRACT_LIBRARY_ID', library_id)
	# normal (non-damage restricted)
	pulldown_parameter_filename_nopath = "{}.normal.parameters".format(library_id)
	normal_parameters = pulldown_parameters_library.replace('SAMPLE_MACRO', library_id) + non_damage_restricted_options
	with open("{}/{}".format(working_directory, pulldown_parameter_filename_nopath), 'w') as pulldown_parameter_file:
		pulldown_parameter_file.write(normal_parameters)
		
	# damage restricted
	damage_restricted_parameters = pulldown_parameters_library.replace('SAMPLE_MACRO', library_id + '_d') + damage_restricted_options
	pulldown_parameter_damage_restricted_filename_nopath = "{}.damage_restricted.parameters".format(library_id)
	with open("{}/{}".format(working_directory, pulldown_parameter_damage_restricted_filename_nopath), 'w') as pulldown_parameter_file:
		pulldown_parameter_file.write(damage_restricted_parameters)
	
	# pulldown
	with open('{}/stdout_normal'.format(working_directory), 'w') as stdout_pulldown, \
		open('{}/stderr_normal'.format(working_directory), 'w') as stderr_pulldown:
			subprocess.run([pulldown_executable, "-p", pulldown_parameter_filename_nopath], check=True, cwd=working_directory, stdout=stdout_pulldown, stderr=stderr_pulldown)
	# damage restricted pulldown
	with open('{}/stdout_damage_restricted'.format(working_directory), 'w') as stdout_pulldown, \
		open('{}/stderr_damage_restricted'.format(working_directory), 'w') as stderr_pulldown:
			subprocess.run([pulldown_executable, "-p", pulldown_parameter_damage_restricted_filename_nopath], check=True, cwd=working_directory, stdout=stdout_pulldown, stderr=stderr_pulldown)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Combine ancient DNA analysis outputs into a single file keyed by index-barcode key. Merges component bams and deduplicates library. Runs pulldown.")
	
	# pulldown is optional
	# 1240k(+) libraries require pulldown. MT libraries do not. 
	parser.add_argument('-p', "--pulldown", help="executable for running pulldown, computing randomized genotypes from bam", nargs=1)
	
	parser.add_argument("bam_list", help="Each line contains the parameters to build a library bam for release. This includes the library ID, the individual ID, experiment, read group description (sequencing run name with experiment type and udg treatment), experiment, and (bam, sequencing run date) pairs ")
	parser.add_argument("adna_jar", help="jar file for assigning read groups to each read in a bam")
	parser.add_argument("picard_jar", help="jar file the Broad Institute Picard toolset, used to merge bams and deduplicate")
	parser.add_argument("report", help="report is used for looking up sex")
	args = parser.parse_args()
	
	bam_list_filename = args.bam_list
	adna_jar_filename = args.adna_jar
	picard_jar = args.picard_jar
	report_filename = args.report
	
	if args.pulldown is not None:
		pulldown_executable = args.pulldown[0]
	# read list of files

	with open(bam_list_filename) as f:
		for line in f:
			parameters_for_library = line.split('\t')
			index_barcode_key = parameters_for_library[0]
			library_id = parameters_for_library[1]
			individual_id = parameters_for_library[2]
			read_group_description = parameters_for_library[3]
			experiment = parameters_for_library[4]
			bam_filenames = parameters_for_library[5::2]
			bam_date_strings = [date_string.strip() for date_string in parameters_for_library[6::2]]
			
			working_directory = library_id
			
			release_library_filename = build_release_library(adna_jar_filename, bam_filenames, bam_date_strings, read_group_description, library_id, individual_id, experiment, picard_jar, working_directory)
			# TODO copy library
			if args.pulldown is not None:
				pulldown(release_library_filename, library_id, report_filename, experiment, pulldown_executable, working_directory, index_barcode_key)
			# TODO copy pulldown results
			# TODO check SNP file to see whether it needs to be copied
