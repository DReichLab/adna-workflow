import pysam
import argparse
import dateutil.parser

def bam_has_read_groups(filename):
	has_read_groups, has_real_library_name, date = read_group_checks(filename)
	return has_read_groups

def read_group_checks(filename):
	has_read_groups = False
	has_real_library_name = False
	date = None

	tag_contents = None
	samfile = pysam.AlignmentFile(filename, "rb")
	try:
		read_groups = samfile.header['RG']
		if len(read_groups) > 0:
			has_read_groups = True
			for read_group in read_groups:
				if read_group['LB'] != 'LB' and len(read_group['LB']) > 0:
					has_real_library_name = True
					break
			for read_group in read_groups:
				if 'DT' in read_group:
					date_string = read_group['DT']
					date = dateutil.parser.parse(date_string).date().strftime("%Y%m%d")
					break
	except:
		pass
	samfile.close()
	#print(tag_contents)
	return has_read_groups, has_real_library_name, date

# does the read group have 

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Test whether a BAM file has any read groups")
	
	parser.add_argument("bam", help="BAM file to inspect for read groups")
	args = parser.parse_args()
	has_read_groups, has_real_library_name, date = read_group_checks(args.bam)
	print('Has read groups: {}'.format(has_read_groups))
	print('Has real library name: {}'.format(has_real_library_name))
	print('Date: {}'.format(date))
