from subprocess import call
import sys
from os.path import basename
# This should be made obsolete by loops in WDL/Cromwell

adna_screen_jar = sys.argv[1]
targets = sys.argv[2]
minimum_mapping_quality = sys.argv[3]
filenames = sys.argv[4:len(sys.argv)]

# TODO convert this into a thread pool
for bam in filenames:
	filename_only = basename(bam)
	histogram = filename_only + ".histogram"
	statistics = filename_only + ".stats"
	with open(statistics, "w") as statistics_file:
		call(["java", "-Xmx2500m", "-jar", adna_screen_jar, "SAMStats", "-f", bam, "-t", targets, "-l", histogram, "-q", minimum_mapping_quality], stdout=statistics_file)
