from subprocess import call
import sys

# Retrieve and cat the requested filename from each of the shard directories of a cromwell run
# This is useful for collecting return code values (rc files), for example
# This is not used directly for the production workflow

filename = sys.argv[1]
numShards = int(sys.argv[2])

for i in range(0,numShards):
	path = "shard-" + str(i) + "/execution/" + filename
	call(["cat", path])
