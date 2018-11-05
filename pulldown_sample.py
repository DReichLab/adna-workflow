import argparse
import sys
import os
from multiprocessing import Pool
from pulldown import create_pulldown_parameter_file, pulldown, merge_pulldowns
from read_groups_from_bam import read_groups_from_bam

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
	def __init__(self, instance_id, udg):
		self.instance_id = instance_id
		self.udg = udg
		self.bam = None
			
# File format is [individual_id] [instance id] [udg] [library list]
# libraries are ignored because we are reusing the file that specifies merges
# return a dictionary of instance_ids -> Instance
def instances_from_file(filename):
	instances = {}
	with open(filename) as f:
		for line in f:
			fields = line.split('\t')
			individual_id = fields[0]
			instance_id = fields[1]
			udg = fields[2]
			if udg not in ALLOWED_UDG_VALUES:
				raise ValueError('Invalid UDG value: {}'.format(udg))
			instance = Instance(instance_id, udg)
			instances[instance_id] = instance
	return instances
''' This is a file containing only keys and sex
def sex_from_file(filename):
	sex_by_instance_id = {}
	with open(filename) as f:
		for line in f:
			fields = line.split('\t')
			instance_id = fields[0]
			sex = fields[1]
			sex_by_instance_id[instance_id] = sex
	return sex_by_instance_id
'''
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

def prepare_pulldown(pulldown_label, instances, sex_by_instance_id, bam_paths):
	# dblist bam paths should be permanent
	with open("{}.dblist".format(pulldown_label), 'w') as db_file:
		for bam_path in bam_paths:
			read_groups = read_groups_from_bam(bam_path)
			bam_filename = os.path.basename(bam_path)
			instance_id = os.path.splitext(bam_filename)[0]

			instance = instances[instance_id]
			if instance_id != instance.instance_id:
				raise ValueError('instance id is not consistent')

			if len(read_groups) > 0:
				db_line = "{0}\t{0}\t{1}\t{2}\n".format(instance_id, bam_path, ":".join(read_groups))
				db_file.write(db_line)
			else:
				print("{} has no read groups".format(instance_id), file=sys.stderr)
	


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
		count = create_individual_file("{}.ind".format(pulldown_base_filename), instances, sex_by_instance_id, udg_type)
		pulldown_parameter_filename_nopath = create_pulldown_parameter_file(pulldown_label, pulldown_base_filename, udg_type, damage_type)
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
	parser.add_argument('-m', "--minus_libraries", help="file with list of UDG minus libraries")
	
	parser.add_argument("sample_bam_list", help="Each line contains the instance id and its list of component bams")
	parser.add_argument("sex", help="sex by instance id")
	parser.add_argument("bams", help="bam files for pulldown, labeled beginning with instance id", nargs='+')
	
	args = parser.parse_args()
	
	instances = instances_from_file(args.sample_bam_list)
	sex_by_instance_id = sex_from_file(args.sex)
	parameter_files = prepare_pulldown(args.pulldown_label, instances, sex_by_instance_id, args.bams)
	
	# run pulldown
	pool = Pool(processes=2)
	results = [pool.apply_async(pulldown, args=(parameter_file, args.pulldown_executable)) for parameter_file in parameter_files]
	pool.close()
	pulldown_file_sets = [result.get() for result in results] # check for exceptions
	pool.join()
		
	
	# merge pulldown results
	merge_pulldowns(args.pulldown_label, pulldown_file_sets)
