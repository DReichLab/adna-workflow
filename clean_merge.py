import sys
import argparse
import shutil
import os
from pathlib import Path
import datetime

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Move released bam files corresponding to a merge list. This is to ensure they get replaced.")
	parser.add_argument("-p", "--parent", help="parent directory for merged files", default="/n/groups/reich/matt/pipeline/sample_merge")
	parser.add_argument("-d", "--directory", help="directory to move files", required=True)
	parser.add_argument("-a", "--age", help="Maximum age of files to move in days. ", type=int)
	parser.add_argument("merge_list", help="[instance_id] [individual_id] [libraries]")
	
	references = ['hg19', 'rsrs']
	
	args = parser.parse_args()
	
	filename = args.merge_list
	
	if args.age:
		allowed_age = datetime.timedelta(days=args.age)
	now = datetime.datetime.now()

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
				print(source)
				if os.path.exists(source):
					if args.age is None:
						shutil.move(source, args.directory)
					else:
						modified = datetime.datetime.fromtimestamp(os.stat(source).st_mtime)
						if (now - modified) < allowed_age: 
							shutil.move(source, args.directory)
						else:
							print('{} is older than requested age', file=sys.stderr)
