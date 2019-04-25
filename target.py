import subprocess
from os.path import basename
import argparse
from multiprocessing import Pool
# This should be made obsolete by loops in WDL/Cromwell

def statistics_for_bam(fullpath_bam, adna_jar_filename, targets, minimum_mapping_quality):
	filename_only = basename(fullpath_bam)
	histogram = filename_only + ".histogram"
	statistics = filename_only + ".stats"
	with open(statistics, "w") as statistics_file:
		subprocess.run(["java", "-Xmx2500m", "-jar", adna_jar_filename, "SAMStats", "-f", bam, "-t", targets, "-l", histogram, "-q", minimum_mapping_quality], stdout=statistics_file, check=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("adna_jar_filename", help="Reich lab jar file")
	parser.add_argument("targets", help="JSON string for how to map reference sequences")
	parser.add_argument("minimum_mapping_quality", help="minimum mapping quality to count sequences", type=int)
	parser.add_argument("-n", "--num_threads", help="minimum mapping quality to count sequences", type=int, default=2)
	parser.add_argument("filenames", help="bam files to analyze", nargs='+')
	args = parser.parse_args()

	adna_jar_filename = args.adna_jar_filename
	targets = args.targets
	minimum_mapping_quality = args.minimum_mapping_quality

	pool = Pool(processes=args.num_threads)
	results = [pool.apply_async(statistics_for_bam, args=(bam, adna_jar_filename, targets, minimum_mapping_quality)) for bam in args.filenames]
	pool.close()
	[result.get() for result in results] # check for exceptions
	pool.join()
