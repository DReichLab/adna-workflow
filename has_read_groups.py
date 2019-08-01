import pysam
import argparse

def bam_has_read_groups(filename):
	has_read_groups, has_real_library_name = read_group_checks(filename)
	return has_read_groups

def read_group_checks(filename, library_length=2):
	has_read_groups = False
	has_real_library_name = False
	tag_contents = None
	samfile = pysam.AlignmentFile(filename, "rb")
	try:
		read_groups = samfile.header['RG']
		if len(read_groups) > 0:
			has_read_groups = True
			for read_group in read_groups:
				if read_group['LB'] != 'LB':
					has_real_library_name = True
					break
	except:
		pass
	samfile.close()
	#print(tag_contents)
	return has_read_groups, has_real_library_name

# does the read group have 

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Test whether a BAM file has any read groups")
	
	parser.add_argument("bam", help="BAM file to inspect for read groups")
	args = parser.parse_args()
	has_read_groups, has_real_library_name = read_group_checks(args.bam)
	print('Has read groups: {}'.format(has_read_groups))
	print('Has real library name: {}'.format(has_real_library_name))
