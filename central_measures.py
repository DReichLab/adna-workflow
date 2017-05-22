# read the median from a histogram file produced by org.apache.commons.math3.stat.Frequency

import sys
import os

def medianFromHistogram(filename):
	with open(filename, "r") as f:
		# ignore the first line
		f.readline()

		for line in f:
			try:
				# remove percent signs
				noPercentLine = line.replace("%", "")
				# read fields
				values = noPercentLine.split()
				length = values[0]
				count = values[1]
				percentage = values[2]
				cumulativePercent = int(values[3])
				if cumulativePercent >= 50:
					return length
			except:
				return 0
	return 0

def meanFromHistogram(filename):
	totalReads = 0
	totalLength = 0
	with open(filename, "r") as f:
		# ignore the first line
		f.readline()
		

		for line in f:
			try:
				# remove percent signs
				noPercentLine = line.replace("%", "")
				# read fields
				values = noPercentLine.split()
				length = int(values[0])
				count = int(values[1])
				percentage = values[2]
				cumulativePercent = int(values[3])
				
				totalReads += count
				totalLength += count * length
			except:
				break;
	return float(totalLength) / float(totalReads)

# first argument is label to be applied to median output, for example: median_rsrs
# second argument is label to be applied to mean output, for example: mean_rsrs
median_label = sys.argv[1]
mean_label = sys.argv[2]
# iterate over passed files
filenames = sys.argv[3: len(sys.argv)]

for filename in filenames:
	key = os.path.splitext(filename)[0] # filename without extension
	median = medianFromHistogram(filename)
	mean = meanFromHistogram(filename)
	print("%s\t%s\t%s\t%s\t%.1f" % (key, median_label, median, mean_label, mean) )
