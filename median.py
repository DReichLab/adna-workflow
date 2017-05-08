# read the median from a histogram file produced by org.apache.commons.math3.stat.Frequency

import sys

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

filename = sys.argv[1]
print (medianFromHistogram(filename))
