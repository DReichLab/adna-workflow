import argparse
from multiprocessing import Pool
from os.path import basename, splitext
import subprocess
		
def deduplicate_bam(source, destination, picard_jar):
	subprocess.run(['java', '-Xmx7g', '-jar', picard_jar, 'MarkDuplicates', 'I={}'.format(source), 'O={}'.format(destination), 'M={}.dedup_stats'.format(destination), 'REMOVE_DUPLICATES=true', 'BARCODE_TAG=XD', 'ADD_PG_TAG_TO_READS=false', 'MAX_FILE_HANDLES=1000'], check=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="For command line deduplication of presorted bams")
	parser.add_argument("-n", "--num_threads", help="", type=int, default=1)
	parser.add_argument("-f", "--file", help="source bam and destination bam on each line", required=True)
	parser.add_argument("--picard", help="Broad Picard jar file", default='/n/groups/reich/matt/pipeline/static/picard-v2.17.10.jar')
	args = parser.parse_args()

	pool = Pool(processes=args.num_threads)
	with open(args.file) as f:
		results = []
		for line in f:
			fields = line.split()
			result = [pool.apply_async(deduplicate_bam, args=(bam, args.directory, args.picard)) for bam in args.bams]
			results.append(result)
	pool.close()
	pool.join()

	for result in results:
		result.get()
