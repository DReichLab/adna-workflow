# read the output of haplogrep to determine a mt haplogroup

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
		notFoundPolysIndex = headerFields.index("Not_Found_Polys")
		foundPolysIndex = headerFields.index("Found_Polys")
		remainingPolysIndex = headerFields.index("Remaining_Polys")

		# second line haplogrep data
		line = f.readline()

		dataFields = line.split("\t")

		key = os.path.splitext(os.path.basename(dataFields[sampleIndex]))[0] # filename without extension
		haplogroup = dataFields[haplogroupIndex]
		overallRank = float(dataFields[rankIndex])
		notFoundPolys = dataFields[notFoundPolysIndex]
		foundPolys = dataFields[foundPolysIndex]
		remainingPolys = dataFields[remainingPolysIndex]

		print(key, end='\t')
		print("MT_Haplogroup\t{}".format(haplogroup), end='\t')
		print("MT_Haplogroup_rank\t{:.3f}".format(overallRank), end='\t')
		print("MT_Haplogroup_NotFoundPolys\t{}".format(notFoundPolys), end='\t')
		print("MT_Haplogroup_FoundPolys\t{}".format(foundPolys), end='\t')
		print("MT_Haplogroup_RemainingPolys\t{}".format(remainingPolys))
