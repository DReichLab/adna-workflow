from multiprocessing import Process, Queue
import subprocess
import sys
from threading import Lock
import os

# We run multiple copies of contammix to retrieve one set of results. 
# We do this because contammix sometimes does not complete and exceeds cluster runtime limits. 
# We continue waiting until we get a result that parses correctly, or all copies have finished. 
numCopies = int(sys.argv[1])
key = sys.argv[2]
contammix_estimate_script = sys.argv[3]
bamFilename = sys.argv[4]
multipleAlignmentFASTA = sys.argv[5]
numChains = sys.argv[6]
minimum_base_quality = sys.argv[7]
deamination_bases_to_clip = sys.argv[8]
seed = int(sys.argv[9])

resultsQueue = Queue(numCopies)

# Failure returns a result, upon which we will stop. 
def contammix(contammix_estimate_script, bamFilename, multipleAlignmentFASTA, numChains, minimum_base_quality, deamination_bases_to_clip, seed, resultsQueue):
	# setup R environment variable to garbage collect more aggressively to restrict memory usage 
	contammix_env = os.environ.copy()
	contammix_env["R_GC_MEM_GROW"] = "0"
	
	result = None
	try:
		filename = "stderr-seed-" + str(seed)
		with open(filename, "w") as f: 
			result = subprocess.run(["Rscript", contammix_estimate_script, "--samFn", bamFilename, "--malnFn", multipleAlignmentFASTA, "--consId", "MT", "--nChains", numChains, "--figure", "diagnostic_plots.pdf", "--baseq", minimum_base_quality, "--trimBases", deamination_bases_to_clip, "--seed", seed, "--tabOutput"], stderr=f, env=contammix_env, check=True, stdout=subprocess.PIPE).stdout
	except:
		pass
	if result == None:
		sys.stderr.write('None')
		resultsQueue.put('\n')
	else:
		result_string = result.decode('utf-8').strip()
		sys.stderr.write(result_string)
		resultsQueue.put(result_string)
	sys.stderr.write('\n')
	
def parseContammixOutput(result):
	fields = result.split('\t')
	inferred_error_rate = float(fields[0])
	map_authentic = float(fields[1])
	lower = float(fields[2])
	upper = float(fields[3])
	gelman_diagnostic = float(fields[4])
	gelman_diagnostic_upper = float(fields[5])
	seed = int(fields[6])
	
	return inferred_error_rate, map_authentic, lower, upper, gelman_diagnostic, gelman_diagnostic_upper, seed

#
processes = [Process(target=contammix,args=(contammix_estimate_script, bamFilename, multipleAlignmentFASTA, numChains, minimum_base_quality, deamination_bases_to_clip, str(seed + i), resultsQueue) ) for i in range(numCopies)]

copyCount = 0
copyCountLock = Lock()
for p in processes:
	p.start()
	with copyCountLock:
		copyCount += 1

inferred_error_rate = float('nan')
map_authentic = float('nan')
lower = float('nan')
upper = float('nan')
gelman_diagnostic = float('nan')
gelman_diagnostic_upper = float('nan')
seed = int(-1)

while copyCount > 0:
	# Queue.get() blocks until next value is available
	firstResult = resultsQueue.get()
	with copyCountLock:
		copyCount -= 1
	
	sys.stderr.write('copyCount: {:d}\n'.format(copyCount))
	sys.stderr.flush()
	try:
		inferred_error_rate, map_authentic, lower, upper, gelman_diagnostic, gelman_diagnostic_upper, seed = parseContammixOutput(firstResult)
		# if this value parses, retrieve values and cancel all other processes
		for p in processes:
			if p.is_alive():
				p.terminate()
		with copyCountLock:
			copyCount = 0
	except:
		pass

print('%s\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%d' % (key, 'contamination_contammix', map_authentic,
					'contamination_contammix_lower', lower, 
					'contamination_contammix_upper', upper,
					'contamination_contammix_gelman', gelman_diagnostic,
					'contamination_contammix_inferred_error', inferred_error_rate,
					'contamination_contammix_seed', seed
					)
)
