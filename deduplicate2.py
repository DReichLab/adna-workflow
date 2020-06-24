import argparse
from multiprocessing import Pool
from pathlib import Path
import subprocess
import sys

from release_libraries import add_read_groups
from library_id import LibraryID
		
def deduplicate_bam(source, destination, picard_jar):
	subprocess.run(['java', '-Xmx7g', '-jar', picard_jar, 'MarkDuplicates', 'I={}'.format(source), 'O={}'.format(destination), 'M={}.dedup_stats'.format(destination), 'REMOVE_DUPLICATES=true', 'BARCODE_TAG=XD', 'ADD_PG_TAG_TO_READS=false', 'MAX_FILE_HANDLES=1000', 'VALIDATION_STRINGENCY=LENIENT'], check=True)

class LibraryInfo:
	pass
# white space delimited fields in input file are:
# input bam file (demultiplexed, non-deduplicated)
# deduplicated bam file (intermediate)
# final bam file
# label (experiment, udg treatment)
# date
# library id

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="For command line deduplication of presorted bams")
	parser.add_argument("-n", "--num_threads", help="", type=int, default=1)
	parser.add_argument("-f", "--file", help="source bam and destination bam on each line", required=True)
	parser.add_argument("-p", "--prefix", help="source prefix")
	parser.add_argument("-w", "--working_directory", help="", default='.')
	
	parser.add_argument("--picard", help="Broad Picard jar file", default='/n/groups/reich/matt/pipeline/static/picard-v2.17.10.jar')
	parser.add_argument("--adna_jar", help="jar file for assigning read groups to each read in a bam", default='/n/groups/reich/matt/pipeline/static/adnascreen-1.10.0-SNAPSHOT.jar')
	args = parser.parse_args()
	
	leniency = "VALIDATION_STRINGENCY=LENIENT"
	jvm_mem_string = '-Xmx7G'

	if args.prefix:
		prefix_path = Path(args.prefix)
	else:
		prefix_path = Path('.')
	
	libraries = []
	with open(args.file) as f:
		for line in f:
			fields = line.split()
			library_info = LibraryInfo()
			library_info.input_bam = fields[0]
			library_info.deduplicated_bam = fields[1]
			library_info.final_bam = fields[2]
			library_info.label = fields[3]
			library_info.date = fields[4]
			library_info.library_id = fields[5]
			libraries.append(library_info)

	# deduplicate
	pool = Pool(processes=args.num_threads)
	results = []
	for library in libraries:
		source = str(prefix_path / library.input_bam)
		destination = library.deduplicated_bam
		print('{}\t{}'.format(source, destination), file=sys.stderr)
		result = pool.apply_async(deduplicate_bam, args=(source, destination, args.picard))
		results.append(result)
	pool.close()
	pool.join()
	for result in results:
		result.get()
		
	# read groups
	pool = Pool(processes=args.num_threads)
	results = []
	for library in libraries:
		destination = fields[1]
		lib_obj = LibraryID(library.library_id)
		individual = 'I{:04d}'.format(lib_obj.sample)
		print('{}'.format(library.library_id), file=sys.stderr)
		result = pool.apply_async(add_read_groups, args=(args.adna_jar, library.deduplicated_bam, library.final_bam, library.date, library.label, library.library_id, individual, args.working_directory, jvm_mem_string, leniency))
		results.append(result)
	pool.close()
	pool.join()
	for result in results:
		result.get()
