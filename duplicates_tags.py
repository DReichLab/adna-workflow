import pysam
import argparse

# test first read to see whether it has an XD tag
def bam_has_XD_tag(filename):
	tag_contents = None
	samfile = pysam.AlignmentFile(filename, "rb")
	for read in samfile:
		if read.has_tag('XD'):
			tag_contents = read.get_tag('XD')
		break
	samfile.close()
	#print(tag_contents)
	return tag_contents is not None

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Test whether a BAM file has XD tags used for barcode-aware deduplication by the Reich Lab pipeline software")
	
	parser.add_argument("bam", help="BAM file to inspect for XD tag")
	args = parser.parse_args()
	result = bam_has_XD_tag(args.bam)
	print(result)
