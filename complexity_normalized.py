import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool

picard_jar = None

def normalized_unique_reads(number_of_reads, bam):
	try:
		downsampled_library = downsample(number_of_reads, bam)
		unique_bam = remove_duplicates(downsampled_library)
		unique_reads = count_bam_reads(unique_bam)
		return (bam, number_of_reads, unique_reads)
	except:
		return (bam, number_of_reads, -1)

def count_bam_reads(filename):
	result = subprocess.run(['samtools', 'view', '-c', filename], check=True, stdout=subprocess.PIPE)
	count = int(result.stdout.strip().decode('utf-8'))
	return count

def downsample(number_of_reads, bam):
	# compute the fraction of reads to retain
	read_count = count_bam_reads(bam)
	fraction = number_of_reads / read_count
	if fraction > 1.0:
		raise ValueError('Cannot upsample {} from {:d} to {:d}'.format(bam, read_count, number_of_reads))
	bam_path = Path(bam)
	if bam_path.suffix != '.bam':
		raise ValueError('Not a BAM {}'.format(bam))
	output = '{}_{:d}.bam'.format(bam_path.stem, number_of_reads) 
	# run Picard to downsample
	subprocess.run(['java', '-Xmx4500m', '-jar', picard_jar, 'DownsampleSam', 'PROBABILITY={:f}'.format(fraction), 'I={}'.format(bam), 'O={}'.format(output)], check=True)
	return output

def remove_duplicates(bam):
	bam_path = Path(bam)
	output = '{}.noduplicates.bam'.format(bam_path.stem)
	# run Picard MarkDuplicates to remove duplicates to get unique read count
	subprocess.run(['java', '-Xmx4500m', '-jar', picard_jar, 'MarkDuplicates', 'I={}'.format(bam), 'O={}'.format(output), 'M={}.{}'.format(output, 'dedup_stats'), 'REMOVE_DUPLICATES=true', 'BARCODE_TAG=XD', 'ADD_PG_TAG_TO_READS=false', 'MAX_FILE_HANDLES=1000'], check=True)
	return output

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Normalize complexity using number of unique reads given a number of reads (hitting targets)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-n', "--num_threads", help="size of thread pool", type=int, default=1)
	parser.add_argument("--picard", help="Broad picard jar", default='/n/groups/reich/matt/pipeline/static/picard-v2.17.10.jar')
	parser.add_argument("bams", help="bam files, filtered to 1240k targets and sorted", nargs='+')
	args = parser.parse_args()

	picard_jar = args.picard
	
	pool = Pool(processes=args.num_threads)
	results = []
	
	for number_of_reads in [5e5, 1e6, 2e6, 4e6]:
		for library in args.bams:
			results.append(pool.apply_async(normalized_unique_reads, args=(int(number_of_reads), library) ) )
			
	pool.close()
	pool.join()
	for result in results:
		values = result.get()
		library_path = Path(values[0])
		print('{}\t{:d}\t{:d}'.format(library_path.stem, int(values[1]), values[2]))
