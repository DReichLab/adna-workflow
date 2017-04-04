import sys
import re

# take a command line argument that uses a Illumina machine output
# read lane name from this filename
filename = sys.argv[1]
result = re.search("L([0-9]{3})", filename)
print(result.group(0))
