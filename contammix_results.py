import sys

key = sys.argv[1]
filename = sys.argv[2]

inferred_error_rate = float('nan')
map_authentic = float('nan')
lower = float('nan')
upper = float('nan')
gelman_diagnostic = float('nan')
gelman_diagnostic_upper = float('nan')

try:
	f = open(filename)
	line = f.readline()
	fields = line.split('\t')
	
	inferred_error_rate = float(fields[0])
	map_authentic = float(fields[1])
	lower = float(fields[2])
	upper = float(fields[3])
	gelman_diagnostic = float(fields[4])
	gelman_diagnostic_upper = float(fields[5])
	f.close()
except:
	pass
finally:
	print('%s\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%f' % (key, 'contamination_contammix', map_authentic,
					   'contamination_contammix_lower', lower, 
					   'contamination_contammix_upper', upper,
					   'contamination_contammix_gelman', gelman_diagnostic,
					   'contamination_contammix_inferred_error', inferred_error_rate,
					   )
	)
