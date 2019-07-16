import argparse
from multiprocessing import Pool
from os.path import basename, splitext
import subprocess
		
def deduplicate_bam(bam, directory, picard_jar):
	sample_id_filename = basename(bam)
	sample_id_filename_no_extension, extension = splitext(sample_id_filename)
	
	sorted_bam = sample_id_filename_no_extension + ".sorted_coordinate.bam"
	
	subprocess.run(['java', '-Xmx3700m', '-jar', picard_jar, 'SortSam', 'I={}'.format(bam), 'O={}'.format(sorted_bam), 'SORT_ORDER=coordinate'], check=True)
	subprocess.run(['java', '-Xmx3700m', '-jar', picard_jar, 'MarkDuplicates', 'I={}'.format(sorted_bam), 'O={}/{}'.format(directory, sample_id_filename), 'M={}.dedup_stats'.format(sample_id_filename), 'REMOVE_DUPLICATES=true', 'BARCODE_TAG=XD', 'ADD_PG_TAG_TO_READS=false', 'MAX_FILE_HANDLES=1000'], check=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="For command line deduplication of a bam")
	parser.add_argument("-n", "--num_threads", help="", type=int, default=1)
	parser.add_argument("-d", "--directory", help="directory to move files", default='deduplicated')
	parser.add_argument("--picard", help="Broad Picard jar file", default='/n/groups/reich/matt/pipeline/static/picard-v2.17.10.jar')
	parser.add_argument("bams", nargs='+', help="BAM files to deduplicate")
	args = parser.parse_args()

	pool = Pool(processes=args.num_threads)
	results = [pool.apply_async(deduplicate_bam, args=(bam, args.directory, args.picard)) for bam in args.bams]
	pool.close()
	pool.join()

	for result in results:
		result.get()
