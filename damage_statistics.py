import fileinput
import os
# read in pmdtools input --first damage statistics from stdin
# 1-2-Q3-Q4.bam
# C>T at first position and SE: 0.1   0.001


# first line is key (not from pmdtools)
f = fileinput.input()
filename = f.readline()
key = os.path.splitext(filename)[0] # filename without extension
# second line is damage from pmdtools
line = f.readline()
# right fields are damage and error
s = line.split(":")
right = s[1]
values = right.split()
damage = float(values[0])
stderr = float(values[1])

print("%s\t%s\t%.3f" % (key, "damage", damage) )
