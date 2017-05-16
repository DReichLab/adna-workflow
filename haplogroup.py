# read the output of haplogrep to determine a mt haplogroup

import fileinput
import os
import sys

# iterate over passed files
filenames = sys.argv[1: len(sys.argv)]

for filename in filenames:
	with open(filename) as f:
		# first line is haplogrep header
		header = f.readline()
		headerFields = header.split("\t")

		sampleIndex = headerFields.index("SampleID")
		haplogroupIndex = headerFields.index("Haplogroup")
		rankIndex = headerFields.index("Overall_Rank")

		# second line haplogrep data
		line = f.readline()

		dataFields = line.split("\t")

		key = os.path.splitext(dataFields[sampleIndex])[0] # filename without extension
		haplogroup = dataFields[haplogroupIndex]
		overallRank = float(dataFields[rankIndex])

		print("%s\tHaplogroup\t%s\tHaplogroup_rank\t%.3f" % (key, haplogroup, overallRank) )
