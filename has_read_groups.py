import pysam
import argparse

# test first read to see whether it has an XD tag
def bam_has_read_groups(filename):
	has_read_groups = False
	tag_contents = None
	samfile = pysam.AlignmentFile(args.bam, "rb")
	try:
		read_groups = samfile.header['RG']
		if len(read_groups) > 0:
			has_read_groups = True
	except:
		pass
	samfile.close()
	#print(tag_contents)
	return has_read_groups

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Test whether a BAM file has any read groups")
	
	parser.add_argument("bam", help="BAM file to inspect for read groups")
	args = parser.parse_args()
	result = bam_has_read_groups(args.bam)
	print(result)
