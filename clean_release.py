import sys
import argparse
import shutil
import os
from pathlib import Path
import datetime

from release_libraries import LibraryParameters

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Move released bam files corresponding to a release list. This is to ensure they get replaced.")
	parser.add_argument("-p", "--parent", help="parent directory for released files", default="/n/groups/reich/matt/pipeline/released_libraries")
	parser.add_argument("-d", "--directory", help="destination directory to move files", required=True)
	parser.add_argument("-a", "--age", help="Maximum age of files to move in days. ", type=int)
	parser.add_argument("release_list", help="bamlist used for assembling bams", nargs='+')
	
	args = parser.parse_args()
	
	bamlists = args.release_list
	if args.age:
		allowed_age = datetime.timedelta(days=args.age)
	now = datetime.datetime.now()

	instance_ids = set()
	parent_path = Path(args.parent)
	
	for bamlist in bamlists:
		with open(bamlist) as f:
			library_parameters = [LibraryParameters(line) for line in f]
	
		for library in library_parameters:
			source = '{}/{}'.format(library.get_release_library_path(args.parent), library.get_release_library_name())
			print(source)
			if os.path.exists(source):
				if args.age is None:
					shutil.move(source, args.directory)
				else:
					modified = datetime.datetime.fromtimestamp(os.stat(source).st_mtime)
					if (now - modified) < allowed_age: 
						shutil.move(source, args.directory)
					else:
						print('{} is older than requested age'.format(source), file=sys.stderr)
