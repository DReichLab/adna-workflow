import sys
import argparse
import shutil
import os
from pathlib import Path

from release_libraries import LibraryParameters

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Move released bam files corresponding to a release list. This is to ensure they get replaced.")
	parser.add_argument("-p", "--parent", help="parent directory for released files", default="/n/groups/reich/matt/pipeline/released_libraries")
	parser.add_argument("-d", "--directory", help="destination directory to move files", required=True)
	parser.add_argument("release_list", help="bamlist used for assembling bams")
	
	args = parser.parse_args()
	
	filename = args.release_list

	instance_ids = set()
	parent_path = Path(args.parent)
	
	with open(filename) as f:
		library_parameters = [LibraryParameters(line) for line in f]
	
	for library in library_parameters:
		source = library.get_release_library_path(args.parent)
		if os.path.exists(source):
			shutil.move(source, args.directory)
