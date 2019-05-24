import sys
import argparse
import shutil
import os
from pathlib import Path

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Move released bam files corresponding to a merge list. This is to ensure they get replaced.")
	parser.add_argument("-p", "--parent", help="parent directory for merged files", default="/n/groups/reich/matt/pipeline/sample_merge")
	parser.add_argument("-d", "--directory", help="directory to move files", required=True)
	parser.add_argument("merge_list", help="[instance_id] [individual_id] [libraries]")
	
	references = ['hg19', 'rsrs']
	
	args = parser.parse_args()
	
	filename = args.merge_list

	instance_ids = set()
	parent_path = Path(args.parent)
	
	with open(filename) as f:
		for line in f:
			fields = line.split('\t')
			instance_id = fields[0]
			individual_id = fields[1]
			libraries = fields[2:]
			
			for reference in references:
				source_path = parent_path / individual_id / "{}.{}.bam".format(instance_id, reference)
				source = str(source_path)
				if os.path.exists(source):
					shutil.move(source, args.directory)
