import argparse
import subprocess

#
def picard_read_group(picard_jar, bam, output_bam, library, platform, flowcell, sample, date, memory):
	read_group_id_string = '{}_{}_{}'.format(library, date, flowcell)
	subprocess.run(['java', 
				 memory,
				 '-jar', 
				 picard_jar, 
				 'AddOrReplaceReadGroups', 
				 'I={}'.format(bam),
				 'O={}'.format(output_bam),
				 'RGID={}'.format(read_group_id_string),
				 'RGLB={}'.format(library),
				 'RGPL={}'.format(platform),
				 'RGPU={}'.format(flowcell),
				 'RGSM={}'.format(sample),
				 'VALIDATION_STRINGENCY=LENIENT',
				 ], check=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Use Picard to write one read group for bam file. Read group ID is based on library name.")

	parser.add_argument('-o', "--output", help="output bam file name")
	parser.add_argument('-l', "--library", help="library for read group", required=True)
	parser.add_argument('-p', "--platform", help="platform for read group", default='ILLUMINA')
	parser.add_argument('-f', "--flowcell", help="flowcell ID for read group", required=True)
	parser.add_argument('-s', "--sample", help="sample ID for read group", required=True)
	parser.add_argument('-d', "--date", help="date for read group ID", required=True)
	parser.add_argument('-j', "--picard_jar", help="Picard jar file", default='/n/groups/reich/matt/pipeline/static/picard-v2.17.10.jar')
	parser.add_argument('-m', "--memory", help="java memory", default='-Xmx7G')
	parser.add_argument("bam", help="bam for read group")
	
	args = parser.parse_args()
	
	bam_filename = args.bam
	picard_read_group(args.picard_jar, args.bam, args.output, args.library, args.platform, args.flowcell, args.sample, args.date, memory=args.memory)
	
