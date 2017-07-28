from multiprocessing import Process, Queue
from subprocess import call
import sys

# We run multiple copies of contammix to retrieve one set of results. 
# We do this because contammix sometimes does not complete and exceeds cluster runtime limits. 
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
	result = None
	try:
		result = check_output(["Rscript", contammix_estimate_script, "--samFn", bamFilename, "--malnFn", multipleAlignmentFASTA, "--consId", "MT", "--nChains", numChains, "--figure", "diagnostic_plots.pdf", "--baseq", minimum_base_quality, "--trimBases", deamination_bases_to_clip, "--seed", seed, "--tabOutput"])
	except:
		pass
	sys.stderr.write(result)
	resultsQueue.put(result)

#
processes = [Process(target=contammix,args=(contammix_estimate_script, bamFilename, multipleAlignmentFASTA, numChains, minimum_base_quality, deamination_bases_to_clip, str(seed + i), resultsQueue) ) for i in range(numCopies)]
	
for p in processes:
	p.start()

# Queue.get() blocks until first value is available
# retrieve this value and cancel all other processes
firstResult = resultsQueue.get()
for p in processes:
	if p.is_alive():
		p.terminate()
		
inferred_error_rate = float('nan')
map_authentic = float('nan')
lower = float('nan')
upper = float('nan')
gelman_diagnostic = float('nan')
gelman_diagnostic_upper = float('nan')
seed = int(-1)

try:
	fields = firstResult.split('\t')
	inferred_error_rate = float(fields[0])
	map_authentic = float(fields[1])
	lower = float(fields[2])
	upper = float(fields[3])
	gelman_diagnostic = float(fields[4])
	gelman_diagnostic_upper = float(fields[5])
	seed = int(fields[6])
except:
	pass

finally:
	print('%s\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%d' % (key, 'contamination_contammix', map_authentic,
					   'contamination_contammix_lower', lower, 
					   'contamination_contammix_upper', upper,
					   'contamination_contammix_gelman', gelman_diagnostic,
					   'contamination_contammix_inferred_error', inferred_error_rate,
					   'contamination_contammix_seed', seed
					   )
	)
