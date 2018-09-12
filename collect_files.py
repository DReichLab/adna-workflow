from pathlib import Path
import argparse
import shutil
import sys

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Collect bam files associated with provided sample list")
	
	parser.add_argument('-f', "--filename", help="Filename containing list of identifiers (for example: library identifiers)")
	parser.add_argument('-p', "--parent_path", help="Parent path containing subdirectories and bam files", default='/n/groups/reich/matt/pipeline/released_libraries')
	parser.add_argument("id_list", nargs='*')
	
	args = parser.parse_args()
	
	id_list = []

	filename = args.filename
	# read list of ids from file
	if filename is not None:
		with open(filename) as f:
			for line in f:
				id_list.append(line.strip())
	
	id_list.extend(args.id_list)
	
	parent_path = Path(args.parent_path)
	for single_id in id_list:
		# copy whole genome and MT bams to current directory
		genome_filename = "{}.1240k_plus.hg19.v1.bam".format(single_id)
		mt_filename = "{}.1240k_plus.rsrs.v1.bam".format(single_id)
		
		genome_src_path = parent_path / single_id / genome_filename
		try:
			shutil.copyfile(genome_src_path, genome_filename)
		except Exception as e:
			print(e, file=sys.stderr)
		
		mt_src_path = parent_path / single_id / mt_filename
		try:
			shutil.copyfile(mt_src_path, mt_filename)
		except Exception as e:
			print(e, file=sys.stderr)
