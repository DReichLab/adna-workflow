import argparse

def read_keyed_values(filename):
	samples = {}
	longest_keys = []
	with open(filename) as f:
		for line in f:
			fields = line.split('\t')
			if len(fields) > 0:
				# to workaround Cromwell problem where task output filenames cannot contain ':'
				# we re-replace alternate character with ':' for analysis
				main_id = fields[0].replace('-', ':')
				keys = fields[1::2]
				values = fields[2::2]
				if len(keys) > len(longest_keys):
					longest_keys = keys
				# 
				if main_id not in samples:
					samples[main_id] = {}
				for key, value in zip(keys, values):
					samples[main_id][key] = value
	return samples, longest_keys
				
def combine(result_sets):
	master = {}
	for result_set in result_sets:
		for main_id, keyed_values in result_set.items():
			if main_id not in master:
				master[main_id] = {}
			current = master[main_id]
			for key, value in keyed_values.items():
				current[key] = value.strip()
					#raise ValueError('{} appears in multiple files'.format(key))
	return master

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument("results_files", help="Each file has one entry per line. Each line begins with an identifier, followed by key-value pairs.", nargs='+')
	args = parser.parse_args()
	
	keyed_value_filenames = args.results_files
	results_set_keys_set_pairs = [read_keyed_values(filename) for filename in keyed_value_filenames]
	results_sets = [results_set for results_set, keys_set in results_set_keys_set_pairs]
	keys_sets = [keys_set for results_set, keys_set in results_set_keys_set_pairs]
	results = combine(results_sets)
	
	# print headers
	headers = [key for keys in keys_sets for key in keys]
	print('\t'.join(['ID'] + headers))
	# print data
	for main_id in results:
		to_print = [main_id]
		instance = results[main_id]
		values = [instance.get(header, '') for header in headers]
		print('\t'.join(to_print + values))
			
