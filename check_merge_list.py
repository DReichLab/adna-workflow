import sys
import argparse


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Check a merge list for pulldown for duplicate instance ids and duplicate libraries within instances.")
	parser.add_argument("merge_list", help="[instance_id] [individual_id] [libraries]")
	
	args = parser.parse_args()
	
	filename = args.merge_list

	instance_ids = set()
	libraries_to_individual = {}
	with open(filename) as f:
		for line in f:
			fields = line.strip().split('\t')
			instance_id = fields[0]
			individual_id = fields[1]
			libraries = fields[2:]
			
			# no instance id should appear more than once
			if instance_id in instance_ids:
				raise ValueError('Duplicate instance id: {}'.format(instance_id))
			else:
				instance_ids.add(instance_id)
				
			library_set = set()
			for library in libraries:
				# libraries should only appear once for any instance id
				if library in library_set:
					print('Library {} appears twice for instance id: {}'.format(library, instance_id), file=sys.stderr)
				else:
					library_set.add(library)
				# libraries should only appear for one individual
				if library in libraries_to_individual and libraries_to_individual[library] != individual_id:
					print('Library {} appears for individuals {} and {}'.format(library, individual_id, libraries_to_individual[library]), file=sys.stderr)
				else:
					libraries_to_individual[library] = individual_id
