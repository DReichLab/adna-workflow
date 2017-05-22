# This flattens a WDL Array[Array[String]] to a Array[String]
# WDL does not have a built-in way to do this easily, so we do this in python
# The elements of each array are assumed to be non-empty filenames.
# Filenames should NOT include "[", "]", or ", "

#sample input:
# [/path/to/A1, /path/to/A2]
# [/path/to/B1, /path/to/B2]

import sys
import re

for line in sys.stdin:
	s = line
	# remove characters [ ], and remove newlines following, if any
	s = re.sub("([\[\]])(\n)?", "", s)
	# elements in array are separated by ", "
	s = s.replace(", ", "\n")
	print (s)
