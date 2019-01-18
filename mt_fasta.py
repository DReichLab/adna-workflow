import subprocess
import argparse
from pathlib import Path
import sys

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Generate the fasta file corresponding to an MT bam.")
	
	parser.add_argument("--htsbox", help="Heng's htsbox program for computing fasta from bam", default='/home/mym11/bin/htsbox')
	parser.add_argument('-Q', "--minimum_base_quality", help="Minimum base quality to be included", type=int, default=20)
	parser.add_argument('-q', "--minimum_mapping_quality", help="Minimum mapping quality to be included", type=int, default=30)
	parser.add_argument('-r', "--reference", help="MT reference file", default='/n/groups/reich/matt/pipeline/static/mtdna_rsrs.fa')
	parser.add_argument('-b', "--bam", help="MT bam file", required=True, nargs='*')
	
	args = parser.parse_args()
	
	htsbox_filename = args.htsbox
	
	for bam in args.bam:
		print(bam, file=sys.stderr)
		output_filename = Path(bam).stem +'.fa'
		with open(output_filename, 'w') as out:
			subprocess.run([htsbox_filename, 'pileup', '-f', args.reference, '-Q', str(args.minimum_base_quality), '-q', str(args.minimum_mapping_quality), '-M', bam], check=True, stdout=out)
