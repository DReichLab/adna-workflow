from bam_finder import getBamPath, library_default_dir, MT_default_dir, ShopVersion
import argparse
import re
import sys
from library_id import LibraryID

def readAnnoFile(anno_filename):
	libraries_by_master_id = {}
	remaps = {}
	with open(anno_filename, errors='surrogateescape') as anno_file:
		headers = anno_file.readline().split('\t')
		num_headers = len(headers)
		# find header indices, and map them to correct IDs
		master_id_index = headers.index('Master ID')
		libraries_index = headers.index('LibraryID(s)')

		for line in anno_file:
			try:
				fields = re.split('\t|\n', line)
				master_id = fields[master_id_index]
				libraries = fields[libraries_index].split(',')
				
				master_id_number = int(master_id[1:]) # remove leading 'I'
				libraries_by_master_id[master_id_number] = libraries
				for library in libraries:
					library_id = LibraryID(library)
					if library_id.sample in remaps and master_id_number != remaps[library_id.sample]:
						raise ValueError('{} maps to {} and {}'.format(library_id.sample, master_id_number, remaps[library_id.sample]))
					remaps[library_id.sample] = master_id_number
			except:
				print(line, file=sys.stderr)
			
	return libraries_by_master_id, remaps

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Augment the bam list for a release with a prior existing version of the library")
	parser.add_argument("-a", "--anno", help="Use anno file for bam hints and read groups", required=True)
	parser.add_argument("libraries", help="Use anno file for bam hints and read groups", nargs='+')
	args = parser.parse_args()
	
	libraries_by_master_id, remaps = readAnnoFile(args.anno)

	master_ids = {}
	for libraries_file in args.libraries:
		with open(libraries_file) as f:
			f.readline() # skip header
			for line in f:
				fields = re.split('\t|\n', line)
				library_id = LibraryID(fields[0])
				master_id_number = library_id.sample
				if master_id_number in remaps:
					master_id_number = remaps[master_id_number]
				if master_id_number not in libraries_by_master_id:
					libraries_by_master_id[master_id_number] = []
				if str(library_id) not in libraries_by_master_id[master_id_number]:
					libraries_by_master_id[master_id_number].append(str(library_id))
				master_ids[master_id_number] = len(libraries_by_master_id[master_id_number])
				
	for master_id, count in master_ids.items():
		if count > 1:
			print('I{}_v40\t{}'.format(master_id, '\t'.join(libraries_by_master_id[master_id])))
