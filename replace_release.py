import sys
import argparse
import shutil
import os
import re
import filecmp
import subprocess
from pathlib import Path

from release_libraries import LibraryParameters

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Move released bam file back to release directory")
	parser.add_argument("-p", "--parent", help="parent directory for released files", default="/n/groups/reich/matt/pipeline/released_libraries")
	parser.add_argument("bams", help="bamlist used for assembling bams", nargs='+')
	
	args = parser.parse_args()
	
	parent_path = Path(args.parent)
	
	for bam in args.bams:
		path = Path(bam)
		filename_only = path.name
		# read off library ID at beginning
		match = re.search('S[0-9]+\.(Y[0-9]+\.)?E[0-9]+\.L[0-9]+[a-z]?', bam)
		library_id = match.group(0)
		print(library_id)
		
		destination = '{}/{}/{}'.format(args.parent, library_id, filename_only)
		if os.path.exists(destination):
			if filecmp.cmp(bam, destination):
				print('{} already present'.format(bam))
			else:
				raise ValueError('{} exists'.format(destination))
		else:
			print('move {} to {}'.format(bam, destination))
			# this gives failures for metadata for files without write access
			#shutil.move(bam, destination)
			subprocess.run(['mv', bam, destination], check=True)
