import argparse
import sys

# File format is [individual_id] [instance id] [library list]
# libraries are ignored because we are reusing the file that specifies merges
# return a dictionary of instance_ids -> library id list
def instances_to_libraries_from_file(filename):
	instances = {}
	instance_to_individual = {}
	with open(filename) as f:
		for line in f:
			fields = line.split()
			instance_id = fields[0]
			individual_id = fields[1]
			library_ids = fields[2:]
			instances[instance_id] = library_ids
			instance_to_individual[instance_id] = individual_id
	return instances, instance_to_individual

def split_pulldowns(instances_to_libraries):
	# construct a list of instances that each library appears in
	library_id_to_instance = {}
	for instance_id in instances_to_libraries:
		for library_id in instances_to_libraries[instance_id]:
			# make a blank list if first appearance of library_id
			if library_id not in library_id_to_instance:
				library_id_to_instance[library_id] = []
			# add instance to list
			library_id_to_instance[library_id].append(instance_id)
	
	# To pass read group checks, we need each read group to appear at most once in each pulldown
	num_pulldowns = max([len(library_id_to_instance[library_id]) for library_id in library_id_to_instance])
	# setup enough blank lists to accomodate the most frequent library
	# we could also increase this for parallelization
	pulldown_instances = []
	for x in range(num_pulldowns):
		pulldown_instances.append([])
	count = 0
	# divide instances among num_pulldowns
	# no library will appear in a pulldown twice
	used_instances = set()
	for library_id, instance_list in library_id_to_instance.items():
		for instance in instance_list:
			if instance not in used_instances:
				pulldown_instances[count % num_pulldowns].append(instance)
				count += 1
				used_instances.add(instance)
	return pulldown_instances

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Prepare pulldown input files for by-sample pulldown batch.")
	
	parser.add_argument("sample_bam_list", help="Each line contains the instance id and its list of component libraries")
	#parser.add_argument("bams", help="bam files for pulldown, labeled beginning with instance id", nargs='+')
	
	args = parser.parse_args()
	
	instances_to_libraries, instances_to_individual = instances_to_libraries_from_file(args.sample_bam_list)
	pulldowns_instances = split_pulldowns(instances_to_libraries)
	pulldown_index = 1
	for pulldown_instances in pulldowns_instances:
		filename = "pulldown_instances{:02d}".format(pulldown_index)
		with open(filename, "w") as f:
			for instance in pulldown_instances:
				print("{}\t{}\t{}".format(instance, instances_to_individual[instance], '\t'.join(instances_to_libraries[instance]).strip()), file=f)
		pulldown_index += 1
	
