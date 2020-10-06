import argparse
from multiprocessing import Pool
import subprocess
		
def samtools_index_bam(bam):
	subprocess.run(['samtools', 'index', bam], check=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="index existing bam files")
	parser.add_argument("-n", "--num_threads", help="", type=int, default=1)
	parser.add_argument("bams", nargs='+', help="BAM files to index")
	args = parser.parse_args()

	pool = Pool(processes=args.num_threads)
	results = [pool.apply_async(samtools_index_bam, args=(bam,)) for bam in args.bams]
	pool.close()
	pool.join()

	for result in results:
		result.get()

