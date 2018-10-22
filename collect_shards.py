import sys
import os

# Retrieve and cat the requested filename from each of the shard directories of a cromwell run
# This is useful for collecting return code values (rc files), for example
# This is not used directly for the production workflow

filename = sys.argv[1]

with os.scandir('.') as top_directory:
	shard_directories = [x for x in top_directory if x.is_dir() and x.name.startswith('shard-')]
	for shard in shard_directories:
		path = shard.path + "/execution/" + filename
		print(shard.name, end='\t')
		try:
			with open(path) as f:
				for line in f:
					print(line.strip(), end='\t')
		except IOError:
			print('NO FILE', end='\t')
		finally:
			print('')
