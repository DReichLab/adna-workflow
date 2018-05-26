import subprocess
import pathlib
import sys
import os
import argparse
import shutil
from pathlib import Path
from multiprocessing import Pool

def add_read_groups(adna_jar_filename, demultiplexed_bam_filename, output_bam_filename, bam_date_string, label, library_id, individual, working_directory):
	subprocess.run(["java", "-Xmx5000m", "-jar", adna_jar_filename,
				 "AssignReadGroups", 
				 "-i", demultiplexed_bam_filename,
				 "-o", output_bam_filename,
				 "-s", individual,
				 "-x", label,
				 "-d", bam_date_string,
				 "-l", library_id], check=True, cwd=working_directory)

def build_release_library(adna_jar_filename, picard_jar, working_directory, library_parameters):
	# make a working directory for this library
	pathlib.Path(working_directory).mkdir(exist_ok=True)
	
	library_id = library_parameters.library_id
	experiment = library_parameters.experiment
	reference = library_parameters.reference
	# add read groups for each library component
	count = 0
	library_component_bams = []
	for i in range(len(library_parameters.bam_filenames)):
		count += 1
		demultiplexed_bam_filename = library_parameters.bam_filenames[i]
		bam_date_string = library_parameters.bam_date_strings[i]
		output_bam_filename = "{0}_{1:d}.{2}.{3}.bam".format(library_id, count, experiment, reference)
		label = "{}_{}".format(library_parameters.read_group_description, library_id)
		add_read_groups(adna_jar_filename, demultiplexed_bam_filename, output_bam_filename, bam_date_string, label, library_id, library_parameters.individual_id, working_directory)
		library_component_bams.append(output_bam_filename)
	
	# use stderr by library
	with open('{}/stdout_build_release_library'.format(working_directory), 'w') as stdout_build, \
		open('{}/stderr_build_release_library'.format(working_directory), 'w') as stderr_build:
		# merge 
		library_with_duplicates_filename = "{0}.{1}.{2}.duplicates.bam".format(library_id, experiment, reference)
		subprocess.run("java -Xmx5500m -jar {} MergeSamFiles I={} O={} SORT_ORDER=coordinate".format(picard_jar, ' I='.join(library_component_bams), library_with_duplicates_filename), shell=True, check=True, cwd=working_directory, stdout=stdout_build, stderr=stderr_build)
		
		# deduplicate
		library_filename = library_parameters.get_release_library_name()
		subprocess.run("java -Xmx5500m -jar {0} MarkDuplicates I={1} O={2} M={2}.dedup_stats REMOVE_DUPLICATES=true BARCODE_TAG=XD ADD_PG_TAG_TO_READS=false MAX_FILE_HANDLES=1000".format(picard_jar, library_with_duplicates_filename, library_filename), shell=True, check=True, cwd=working_directory, stdout=stdout_build, stderr=stderr_build)
		
	return library_filename

def copy_release_library(library_parameters, destination_parent_directory, release_library_name, working_directory):
	if destination_parent_directory != None:
		# create library directory if it does not exist
		destination_library_path = library_parameters.get_release_library_path(destination_parent_directory)
		pathlib.Path(destination_library_path).mkdir(exist_ok=True)
		source_file = Path(working_directory) / release_library_name
		# handle case where this library already exists
		if not os.path.exists(source_file):
			shutil.copy(source_file, destination_library_path)
		else:
			sys.stderr.write('{} already exists'.format(source_file))

def index_library(library_parameters, release_parent_directory, release_library_name):
	if release_parent_directory != None:
		release_library_path = library_parameters.get_release_library_path(release_parent_directory)
		subprocess.run(['samtools', 'index', release_library_name], check=True, cwd=release_library_path)
	
class LibraryParameters:
	def __init__(self, line):
		parameters_for_library = line.split('\t')
		self.index_barcode_key = parameters_for_library[0]
		self.library_id = parameters_for_library[1]
		self.individual_id = parameters_for_library[2]
		self.read_group_description = parameters_for_library[3]
		self.experiment = parameters_for_library[4]
		self.udg = parameters_for_library[5]
		self.reference = parameters_for_library[6]
		self.version = int(parameters_for_library[7])
		self.wetlab_notes = parameters_for_library[8]
		self.bam_filenames = parameters_for_library[9::2]
		self.bam_date_strings = [date_string.strip() for date_string in parameters_for_library[10::2]]
		
		if len(self.bam_filenames) != len(self.bam_date_strings):
			raise ValueError('Each bam needs a date string')
	
	def get_release_library_name(self):
		release_library_name = "{0}.{1}.{2}.v{3:d}.bam".format(self.library_id, self.experiment, self.reference, self.version)
		return release_library_name
		
	def get_release_library_path(self, release_directory):
		release_library_path = "{}/{}".format(release_directory, self.library_id)
		return release_library_path
	
def process_library(adna_jar_filename, picard_jar, release_directory, library_parameters):
	working_directory = library_parameters.library_id
	
	release_library_filename = build_release_library(adna_jar_filename, picard_jar, working_directory, library_parameters)
	copy_release_library(library_parameters, release_directory, release_library_filename, working_directory)
	index_library(library_parameters, release_directory, release_library_filename)
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Combine ancient DNA analysis outputs into a single file keyed by index-barcode key. Merges component bams and deduplicates library. Copy library to release directory and index it in preparation for pulldown.")
	
	parser.add_argument('-n', "--num_threads", help="size of thread pool", type=int, default =10)
	parser.add_argument('-d', "--release_directory", help="parent directory to put released libraries")
	
	parser.add_argument("bam_list", help="Each line contains the parameters to build a library bam for release. This includes the library ID, the individual ID, experiment, read group description (sequencing run name with experiment type and udg treatment), experiment, and (bam, sequencing run date) pairs ")
	parser.add_argument("adna_jar", help="jar file for assigning read groups to each read in a bam")
	parser.add_argument("picard_jar", help="jar file the Broad Institute Picard toolset, used to merge bams and deduplicate")
	args = parser.parse_args()
	
	bam_list_filename = args.bam_list
	# read list of files
	with open(bam_list_filename) as f:
		library_parameters = [LibraryParameters(line) for line in f]
	# build libraries
	adna_jar_filename = args.adna_jar
	picard_jar = args.picard_jar
	
	pool = Pool(processes=args.num_threads)
	results = [pool.apply_async(process_library, args=(adna_jar_filename, picard_jar, args.release_directory, parameters)) for parameters in library_parameters]
	pool.close()
	[result.get() for result in results] # check for exceptions
	pool.join()
	
