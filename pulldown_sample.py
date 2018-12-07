import argparse
import sys
import os
from multiprocessing import Pool
from pulldown import create_pulldown_parameter_file, pulldown, merge_pulldowns
from read_groups_from_bam import read_groups_and_libraries_from_bam

# return a set of the library ids in the given file
def read_library_ids(filename):
	library_ids = set()
	if filename != None:
		with open(filename) as f:
			for line in f:
				library_id = line.strip()
				library_ids.add(library_id)
	return library_ids

def create_individual_file(filename, instances, sex_by_instance_id, udg_filter):
	count = 0
	with open(filename, 'w') as individual_file:
		for instance_id, instance in instances.items():
			if udg_filter == instance.udg:
				if instance_id != instance.instance_id:
					raise ValueError('instance id is not consistent')
				individual_file.write("{0}\t{1}\t{0}\n".format(instance_id, sex_by_instance_id[instance_id]))
				count += 1
	return count

class Instance:
	def __init__(self, instance_id):
		self.instance_id = instance_id
		self.bam = None
			
# File format is [individual_id] [instance id] [library list]
# libraries are ignored because we are reusing the file that specifies merges
# return a dictionary of instance_ids -> Instance
def instances_from_file(filename):
	instances = {}
	with open(filename) as f:
		for line in f:
			fields = line.split('\t')
			instance_id = fields[0]
			individual_id = fields[1]
			instance = Instance(instance_id)
			instances[instance_id] = instance
	return instances

# Read an adna keyed value aggregated statistics file by line to find sex
def sex_from_file(filename):
	sex_by_instance_id = {}
	with open(filename) as f:
		for line in f:
			fields = line.split('\t')
			instance_id = fields[0]
			keys = fields[1::2]
			values = fields[2::2]
			sex_index = keys.index('1240k_post_sex')
			sex = values[sex_index]
			
			sex_by_instance_id[instance_id] = sex.strip()
	return sex_by_instance_id

# based on lists of minus and plus libraries, determine UDG treatment
# a library has exactly one UDG treatment
# if a library is not minus or plus, then it is half
def udg_for_library(library_id, minus_libraries, plus_libraries):
	is_minus = library_id in minus_libraries
	is_plus = library_id in plus_libraries
	
	if is_minus and is_plus:
		raise ValueError('library cannot have multiple UDG treatments')
	elif is_minus:
		return MINUS
	elif is_plus:
		return PLUS
	else:
		return HALF

def prepare_pulldown(pulldown_label, instances, sex_by_instance_id, bam_paths, minus_libraries, plus_libraries):
	# generate one dblist file per UDG treatment
	# dblist bam paths should be permanent
	for udg in ALLOWED_UDG_VALUES:
		with open("{}.{}.dblist".format(pulldown_label, udg), 'w') as db_file:
			for bam_path in bam_paths:
				read_groups_to_libraries = read_groups_and_libraries_from_bam(bam_path)
				bam_filename = os.path.basename(bam_path)
				# filename contains reference, for example I8861.hg19, remove this
				instance_id_with_reference = os.path.splitext(bam_filename)[0]
				instance_id = instance_id_with_reference.rsplit('.', 1)[0]

				instance = instances[instance_id]
				if instance_id != instance.instance_id:
					raise ValueError('instance id is not consistent')
				
				# only include read groups for libraries which match the UDG treatment
				read_groups_matching_udg = []
				for read_group in sorted(read_groups_to_libraries.keys()):
					library_id = read_groups_to_libraries[read_group]
					if udg_for_library(library_id) == udg:
						read_groups_matching_udg.append(read_group)

				if len(read_groups_matching_udg) > 0:
					db_line = "{0}\t{0}\t{1}\t{2}\n".format(instance_id, os.path.realpath(bam_path), ":".join(read_groups_matching_udg))
					db_file.write(db_line)
				else:
					print("{} has no read groups".format(instance_id), file=sys.stderr)
	

	# order matters for preference of results when merging mixed UDG reads
	parameter_indices = [
		(PLUS, NORMAL) # UDG plus cannot be damage restricted
		(HALF, NORMAL),
		(HALF, DAMAGE_RESTRICTED),
		(MINUS, NORMAL),
		(MINUS, DAMAGE_RESTRICTED),
		]
	
	parameter_file_outputs = []
	for (udg_type, damage_type) in parameter_indices:
		pulldown_base_filename = "{}.{}.{}".format(pulldown_label, udg_type, damage_type)
		count = create_individual_file("{}.ind".format(pulldown_base_filename), instances, sex_by_instance_id, udg_type)
		pulldown_parameter_filename_nopath = create_pulldown_parameter_file('{}.{}'.format(pulldown_label, udg_type), pulldown_base_filename, udg_type, damage_type)
		if count > 0:
			parameter_file_outputs.append(pulldown_parameter_filename_nopath)
	return parameter_file_outputs

# Different UDG treatments determine pulldown parameters
# Separate pulldown by udg type
MINUS = 'minus'
HALF = 'half'
PLUS = 'plus'
ALLOWED_UDG_VALUES = [MINUS, HALF, PLUS]
NORMAL = 'normal'
DAMAGE_RESTRICTED = 'damage_restricted'
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Prepare pulldown input files for by-sample pulldown batch.")
	
	parser.add_argument('-p', "--pulldown_executable", help="executable to run pulldown")
	parser.add_argument('-l', "--pulldown_label", help="label for pulldown filenames")
	parser.add_argument('-r', "--release_directory", help="parent directory to read released libraries")
	parser.add_argument("--minus_libraries", help="file with list of UDG minus libraries")
	parser.add_argument("--plus_libraries", help="file with list of UDG minus libraries")
	
	parser.add_argument("sample_bam_list", help="Each line contains the instance id and its list of component bams")
	parser.add_argument("sex", help="sex by instance id")
	parser.add_argument("bams", help="bam files for pulldown, labeled beginning with instance id", nargs='+')
	
	args = parser.parse_args()
	
	instances = instances_from_file(args.sample_bam_list)
	sex_by_instance_id = sex_from_file(args.sex)
	minus_libraries = read_library_ids(args.minus_libraries)
	plus_libraries = read_library_ids(args.plus_libraries)
	parameter_files = prepare_pulldown(args.pulldown_label, instances, sex_by_instance_id, args.bams, minus_libraries, plus_libraries)
	
	# run pulldown
	pool = Pool(processes=2)
	results = [pool.apply_async(pulldown, args=(parameter_file, args.pulldown_executable)) for parameter_file in parameter_files]
	pool.close()
	pool.join()
	pulldown_file_sets = [result.get() for result in results] # check for exceptions
		
	
	# merge pulldown results
	merge_pulldowns(args.pulldown_label, pulldown_file_sets)
