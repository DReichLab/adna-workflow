from release_libraries import LibraryParameters
from bam_finder import getBamPath, library_default_dir, MT_default_dir
import argparse
import re

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Augment the bam list for a release with a prior existing version of the library")
	parser.add_argument("bam_list", help="Each line contains the parameters to build a library bam for release. This includes the library ID, the individual ID, experiment, read group description (sequencing run name with experiment type and udg treatment), experiment, and (bam, sequencing run date) pairs ")
	args = parser.parse_args()

	with open(args.bam_list) as f:
		library_parameters = [LibraryParameters(line) for line in f]

	for x in library_parameters:
		experiment = x.experiment
		if '1240k' in experiment:
			experiment = '1240k'
		search_directory = MT_default_dir if x.reference == 'rsrs' else library_default_dir
		existingBAM = getBamPath(x.library_id, experiment=experiment, reference=x.reference, version_policy='latest', shop_parent_directory=search_directory)
		bam = str(existingBAM)
		if len(bam) > 0:
			try: # this will match a new pipeline bam
				match = re.search('v([0-9]+).bam', bam)
				new_version = int(match.group(1)) + 1
			except: # if the existing version is Shop's
				new_version = 1
			#print('{}\t{}\t{:d}'.format(x.library_id, bam, new_version))
			x.version = new_version
			x.bam_filenames.append(str(existingBAM))
			x.bam_date_strings.append('noDate') # the bam date string is used for generating read groups, which the existing bam does not need
		#print('{}\t{}'.format(x.library_id, bam))
		print(x)
