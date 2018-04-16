import subprocess
import sys
import re

# parse the read group strings from a bam/sam header
# return array of read group strings
def read_groups_from_bam(bam_filename):
	bam_header = subprocess.run(["samtools", "view", "-H", bam_filename], stdout=subprocess.PIPE, check=True)
	#print(bam_header)
	bam_header_lines = bam_header.stdout.decode('utf-8').split('\n')
	read_groups_lines = [line for line in bam_header_lines if line.startswith('@RG')]
	#print(read_groups_lines)
	read_groups = []
	for line in read_groups_lines:
		pattern = re.compile('@RG\s+ID:(\S+)')
		match_obj = pattern.match(line)
		read_group = match_obj.group(1)
		read_groups.append(read_group)
	sorted_read_groups = sorted(read_groups)
	return sorted_read_groups

if __name__ == "__main__":
	bam_filename = sys.argv[1]
	read_groups = read_groups_from_bam(bam_filename)
	for read_group in read_groups:
		print(read_group)
