import argparse
import pysam

# parse the read group strings from a bam/sam header
# return array of read group strings
def read_groups_from_bam(bam_filename, use_libraries):
	bam = pysam.AlignmentFile(bam_filename, "rb")
	header = bam.header
	
	read_groups = header['RG']
	
	if use_libraries:
		field = 'LB'
	else:
		field = 'ID'
		
	#print(read_groups)
	results = {}
	for read_group in read_groups:
		results[read_group[field]] = 1
		#read_group['SM'] = sample
		#print(read_group)
	results_without_duplicates = [key for (key, ignored) in results.items()]
	
	sorted_read_groups = sorted(results_without_duplicates)
	return sorted_read_groups

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Prepare pulldown input files for by-sample pulldown batch.")

	parser.add_argument('-l', "--libraries", help="report libraries instead of read groups", action='store_true')
	parser.add_argument("bam", help="bam for read groups")
	
	args = parser.parse_args()
	
	bam_filename = args.bam
	read_groups = read_groups_from_bam(bam_filename, args.libraries)
	for read_group in read_groups:
		print(read_group)
