import subprocess
import pathlib
import sys
import os
import argparse
import shutil
from pathlib import Path
from multiprocessing import Pool
from duplicates_tags import bam_has_XD_tag
from has_read_groups import read_group_checks
from bam_finder import ShopVersion

# Check to see if a bam has any reads
def bam_has_aligned_reads(bam_filename):
	result = subprocess.run(['samtools', 'view', '-F', '4', '-c', bam_filename], check=True, universal_newlines=True, stdout=subprocess.PIPE).stdout.strip()
	return int(result) > 0

def aligned_reads_only(input_bam, output_bam):
	subprocess.run(['samtools', 'view', '-h', '-b', '-F', '4', '-o', output_bam, input_bam], check=True)

def add_read_groups(adna_jar_filename, demultiplexed_bam_filename, output_bam_filename, bam_date_string, label, library_id, individual, working_directory, jvm_mem_string, leniency):
	command_read_groups = ["java", jvm_mem_string, "-jar", adna_jar_filename,
				 "AssignReadGroups", 
				 "-i", demultiplexed_bam_filename,
				 "-o", output_bam_filename,
				 "-s", individual,
				 "-x", label,
				 "-d", bam_date_string,
				 "-l", library_id]
	if leniency:
		command_read_groups += ['--lenient']
	subprocess.run(command_read_groups, check=True, cwd=working_directory)

def build_release_library(adna_jar_filename, picard_jar, working_directory, library_parameters, jvm_mem_string, leniency):
	# make a working directory for this library
	pathlib.Path(working_directory).mkdir(exist_ok=True)
	
	library_id = library_parameters.library_id
	experiment = library_parameters.experiment
	reference = library_parameters.reference
	# add read groups for each library component
	count = 0
	library_component_bams = []
	component_bam_missing_duplicate_tag = False
	for input_bam, bam_date_string in zip(library_parameters.bam_filenames, library_parameters.bam_date_strings):
		# only bams with reads need to be merged
		if bam_has_aligned_reads(input_bam):
			if not bam_has_XD_tag(input_bam):
				component_bam_missing_duplicate_tag = True
			count += 1
			# filter aligned reads only
			component_bam_filename = '{}_{:d}.bam'.format(Path(input_bam).stem, count)
			component_bam_path = '{}/{}'.format(working_directory, component_bam_filename)
			aligned_reads_only(input_bam, component_bam_path)
			
			output_bam_filename = "{0}_{1:d}.{2}.{3}.bam".format(library_id, count, experiment, reference)
			# Demultiplexed, but unreleased bams need read groups added
			has_read_groups, has_real_library_name, date_string = read_group_checks(component_bam_path)
			if not has_read_groups:
				label = "{}_{}".format(library_parameters.read_group_description, library_id)
				add_read_groups(adna_jar_filename, component_bam_filename, output_bam_filename, bam_date_string, label, library_id, library_parameters.individual_id, working_directory, jvm_mem_string, leniency)
			# Shop's bams need read groups rewritten
			elif not has_real_library_name:
				shop_version = ShopVersion(component_bam_filename)
				bam_date_string = shop_version.date_string
				label = "{}_{}".format(library_parameters.read_group_description, library_id)
				add_read_groups(adna_jar_filename, component_bam_filename, output_bam_filename, bam_date_string, label, library_id, library_parameters.individual_id, working_directory, jvm_mem_string)
			# Previously released libraries already have read groups, and do not need them added
			# Simply include these in the merge
			else:
				os.symlink(component_bam_filename, '{}/{}'.format(working_directory, output_bam_filename))
			library_component_bams.append(output_bam_filename)
	
	# use stderr by library
	with open('{}/stdout_build_release_library'.format(working_directory), 'w') as stdout_build, \
		open('{}/stderr_build_release_library'.format(working_directory), 'w') as stderr_build:
			
		leniency_string = "VALIDATION_STRINGENCY=LENIENT" if leniency else ""
		
		library_filename = library_parameters.get_release_library_name()
		if len(library_component_bams) > 0: # merge any bams with reads and mark duplicates
			# merge 
			library_with_duplicates_filename = "{0}.{1}.{2}.duplicates.bam".format(library_id, experiment, reference)
			subprocess.run("java {} -jar {} MergeSamFiles I={} O={} SORT_ORDER=coordinate {}".format(jvm_mem_string, picard_jar, ' I='.join(library_component_bams), library_with_duplicates_filename, leniency_string), shell=True, check=True, cwd=working_directory, stdout=stdout_build, stderr=stderr_build)
			
			# if any component bam is missing the XD tag, we cannot use the tag to deduplicate properly with barcodes
			# We treat this conservatively by removing all barcode information and deduplicating based on position and length
			if component_bam_missing_duplicate_tag:
				library_with_duplicates_tag_rewritten_filename = "{0}.{1}.{2}.duplicates.tagxd.bam".format(library_id, experiment, reference)
				subprocess.run(['java', jvm_mem_string, '-jar', adna_jar_filename, 'DuplicatesTagRewrite', 
					'-i', library_with_duplicates_filename,
					'-o', library_with_duplicates_tag_rewritten_filename], check=True, cwd=working_directory)
				to_deduplicate_filename = library_with_duplicates_tag_rewritten_filename
			else:
				to_deduplicate_filename = library_with_duplicates_filename
			
			# deduplicate
			subprocess.run("java {0} -jar {1} MarkDuplicates I={2} O={3} M={3}.dedup_stats REMOVE_DUPLICATES=true BARCODE_TAG=XD ADD_PG_TAG_TO_READS=false MAX_FILE_HANDLES=1000 COMPRESSION_LEVEL=9 {4}".format(jvm_mem_string, picard_jar, to_deduplicate_filename, library_filename, leniency_string), shell=True, check=True, cwd=working_directory, stdout=stdout_build, stderr=stderr_build)
		else: # There are no reads, so use an empty bam. First bam should exist and be empty, so return a copy of that, ensuring aligned reads only
			aligned_reads_only(library_parameters.bam_filenames[0], working_directory + '/' + library_filename)
		
	return library_filename

def copy_release_library(library_parameters, destination_parent_directory, release_library_name, working_directory):
	if destination_parent_directory != None:
		# create library directory if it does not exist
		destination_library_path = library_parameters.get_release_library_path(destination_parent_directory)
		pathlib.Path(destination_library_path).mkdir(mode=0o775, exist_ok=True)
		source_file = Path(working_directory) / release_library_name
		# handle case where this library already exists
		if not os.path.exists(destination_library_path + "/" + release_library_name):
			created = shutil.copy(source_file, destination_library_path)
			os.chmod(created, 0o664)
		else:
			sys.stderr.write('{} already exists\n'.format(source_file))
			raise FileExistsError(source_file)

def index_library(library_parameters, release_parent_directory, release_library_name):
	if release_parent_directory != None:
		release_library_path = library_parameters.get_release_library_path(release_parent_directory)
		# indexing empty files fails
		if bam_has_aligned_reads('{}/{}'.format(release_library_path, release_library_name)):
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
		
	def __repr__(self):
		parameters = [self.index_barcode_key,
				self.library_id,
				self.individual_id,
				self.read_group_description,
				self.experiment,
				self.udg,
				self.reference,
				str(self.version),
				self.wetlab_notes]
		interleaved = self.bam_filenames + self.bam_date_strings
		interleaved[::2] = self.bam_filenames
		interleaved[1::2] = self.bam_date_strings
		return '\t'.join(parameters + interleaved)
	
	def get_release_library_name(self):
		release_library_name = "{0}.{1}.{2}.v{3:d}.bam".format(self.library_id, self.experiment, self.reference, self.version)
		return release_library_name
		
	def get_release_library_path(self, release_directory):
		release_library_path = "{}/{}".format(release_directory, self.library_id)
		return release_library_path
	
def process_library(adna_jar_filename, picard_jar, release_directory, library_parameters, jvm_mem_string, leniency):
	working_directory = library_parameters.library_id
	
	release_library_filename = build_release_library(adna_jar_filename, picard_jar, working_directory, library_parameters, jvm_mem_string, leniency)
	copy_release_library(library_parameters, release_directory, release_library_filename, working_directory)
	index_library(library_parameters, release_directory, release_library_filename)
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Combine ancient DNA analysis outputs into a single file keyed by index-barcode key. Merges component bams and deduplicates library. Copy library to release directory and index it in preparation for pulldown.")
	
	parser.add_argument('-n', "--num_threads", help="size of thread pool", type=int, default =10)
	parser.add_argument('-d', "--release_directory", help="parent directory to put released libraries")
	parser.add_argument('-m', "--mem", help='memory string for java virtual machine, for example "5000m" or "10G"', default='5500m')
	parser.add_argument('-l', "--lenient", action='store_true', help='Picard leniency flag for unmapped reads with > 0 mapping quality')
	
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
	jvm_mem_string = '-Xmx{}'.format(args.mem)
	leniency = args.lenient
	
	pool = Pool(processes=args.num_threads)
	results = [pool.apply_async(process_library, args=(adna_jar_filename, picard_jar, args.release_directory, parameters, jvm_mem_string, leniency)) for parameters in library_parameters]
	pool.close()
	[result.get() for result in results] # check for exceptions
	pool.join()
	
