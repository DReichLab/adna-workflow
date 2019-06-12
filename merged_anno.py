import argparse
import sys
from library_results import read_library_file, read_sample_file
from bam_finder import find_pipeline_bam
from library_id import LibraryID
from pathlib import Path
import re

class InstanceAnnoEntry:
	def __init__(self, instance_id):
		self.index = ''
		self.instance_id = instance_id
		self.individual_id = ''
		self.skeletal_code = ''
		self.skeletal_element = ''
		self.library_ids = [] = ''
		self.data_type = 'bam'
		self.pulldown_logfile = ''
		self.pulldown_1_sample_id = ''
		self.pulldown_2_alt_sample_id = ''
		self.pulldown_3_bam = ''
		self.pulldown_4_hefta = ''
		self.pulldown_5_diploid_source = ''
		self.publication = ''
		self.collaborator = ''
		self.date_formatted_check = ''
		self.date = ''
		self.date_notes = ''
		self.group_label = ''
		self.location = ''
		self.site = '..'
		self.country = ''
		self.latitude = ''
		self.longitude = ''
		self.sex = ''
		self.mtdna_haplogroup = ''
		self.mtdna_coverage = 0
		self.y_haplogroup = ''
		self.y_curation = ''
		self.best_shotgun_endogenous = ''
		self.coverage = 0.0 # mean depth
		self.autosome_snps = 0
		self.udg = ''
		self.damage_restricted = 'All' # manual
		self.assessment = ''
		self.carbon14_data_status = ''
		self.more_libraries = ''
		self.publication_plan = '' # similar to, but distinct from publication
		self.mt_bam = ''
		self.mt_fasta = ''
		self.x_contam_point = ''
		self.x_contam_z = ''
		self.endogenous_by_library = []
		self.damage_by_library = []
		self.mt_haplogroup_by_library = []
		self.mt_contamination_by_library = []
		self.year_published = ''
		
	def __repr__(self):
		fields = []
		fields.append(self.index)
		fields.append(self.instance_id)
		fields.append(self.individual_id)
		fields.append(self.skeletal_code)
		fields.append(self.skeletal_element)
		fields.append(','.join(self.library_ids))
		fields.append('{:d}'.format(len(self.library_ids)))
		fields.append(self.data_type)
		if self.pulldown_logfile != None:
			fields.append(self.pulldown_logfile)
		else:
			fields.append('')
		fields.append(self.pulldown_1_sample_id)
		fields.append(self.pulldown_2_alt_sample_id)
		fields.append(self.pulldown_3_bam)
		fields.append(self.pulldown_4_hefta)
		fields.append(self.pulldown_5_diploid_source)
		fields.append(self.publication)
		fields.append(self.collaborator)
		fields.append(self.date_formatted_check)
		fields.append(self.date)
		fields.append(self.date_notes)
		fields.append(self.group_label)
		fields.append(self.location)
		fields.append(self.site)
		fields.append(self.country)
		fields.append(self.latitude)
		fields.append(self.longitude)
		fields.append(self.sex)
		if self.mtdna_coverage >= 2.0:
			fields.append(self.mtdna_haplogroup)
		else:
			fields.append('')
		fields.append(self.y_haplogroup)
		fields.append(self.y_curation)
		if self.best_shotgun_endogenous != '':
			fields.append('{:.3f}'.format(float(self.best_shotgun_endogenous)))
		else:
			fields.append('')
		fields.append('{:.3f}'.format(self.coverage)) # coverage = mean depth
		fields.append('{:d}'.format(self.autosome_snps))
		fields.append(','.join(self.udg))
		fields.append(self.damage_restricted)
		fields.append(self.assessment)
		fields.append(self.carbon14_data_status)
		fields.append(self.more_libraries)
		fields.append(self.publication_plan)
		fields.append(self.mt_bam)
		fields.append(self.mt_fasta)
		fields.append(self.x_contam_point)
		fields.append(self.x_contam_z)
		fields.append(','.join(self.endogenous_by_library))
		fields.append(','.join(self.damage_by_library))
		fields.append(','.join(self.mt_haplogroup_by_library))
		fields.append(','.join(self.mt_contamination_by_library))
		fields.append(self.year_published)
		
		return '\t'.join(fields)
		
def library_id_in_pulldown(library_id, names, pulldown_dir_root):
	possible_udg = ['half', 'minus', 'plus']
	for name in names:
		for udg in possible_udg:
			candidate = '{0}/{0}.{1}.normal.ind'.format(name, udg)
			with open(candidate) as f:
				for line in f:
					if library_id in line:
						logfile = '{2}/{0}/{0}.{1}.normal.parameters.stdout'.format(name, udg, pulldown_dir_root)
						bamfile = find_pipeline_bam(library_id, 'hg19', '1240k', version_policy='latest')
						return logfile, bamfile
	return None, None

def read_merge_summary(filename):
	info = {}
	with open(filename) as f:
		header_line = f.readline()
		headers = re.split('\t|\n', header_line)
		for line in f:
			fields = re.split('\t|\n', line)
			instance_id = fields[headers.index('ID')]
			info[instance_id] = fields
	return headers, info

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Combine pipeline analysis results with Rebecca's library file.")
	parser.add_argument('-l', "--library_file", help="Rebecca's library file", required=True)
	parser.add_argument('-m', "--merge_lists", help="", nargs='*')
	parser.add_argument('-r', "--reports", help="Tab-delimited software pipeline analysis results", nargs='*')
	parser.add_argument('-n', "--names", help="Capture names to look for pulldown logs", nargs='*')
	parser.add_argument('-p', "--pulldown_dir_root", help="Pipeline pulldown directory", required=True)
	parser.add_argument('-s', "--sample_file", help="Sample file describing samples", required=True)
	parser.add_argument('-x', "--delimiter", help="Sample file field delimiter", default='\t')
	parser.add_argument('-z', "--merge_analysis", help="File containing tab-delimited merge analysis")
	args = parser.parse_args()
	
	# read library info into memory
	library_headers, library_ids, library_info = read_library_file(args.library_file)
	# read sample info into memory
	sample_headers, sample_info = read_sample_file(args.sample_file, args.delimiter)
	
	if args.merge_analysis is not None:
		merge_analysis_headers, merge_analysis_info = read_merge_summary(args.merge_analysis)
	
	released_library_list = {}
	instance_list = {}
	# Each separate library gets its own anno entry
	if args.reports is not None:
		for report in args.reports:
			with open(report) as f:
				f.readline() # skip read count
				header_line = f.readline()
				headers = re.split('\t|\n', header_line)
				for line in f:
					fields = re.split('\t|\n', line)
					library_id = fields[headers.index('library_id')]
					experiment = fields[headers.index('experiment')]
					if '1240k' in experiment and 'Contl' not in library_id:
						if library_id in released_library_list:
							#raise ValueError('duplicate library: {}'.format(library_id))
							print('duplicate library: {}'.format(library_id), file=sys.stderr)
						else:
							# store dictionary of fields for library
							library_fields = {}
							for (header, field) in zip(headers, fields):
								library_fields[header] = field
							released_library_list[library_id] = library_fields
							instance = InstanceAnnoEntry(library_id)
							instance.library_ids = [library_id]
							instance.mtdna_haplogroup = library_info[library_id][library_headers.index('mtDNA_Haplogroup')]
							mtdna_coverage = library_info[library_id][library_headers.index('mtDNA_Coverage')]
							if mtdna_coverage != '':
								instance.mtdna_coverage = float(library_info[library_id][library_headers.index('mtDNA_Coverage')])
							instance.best_shotgun_endogenous = library_info[library_id][library_headers.index('Shotgun_Percent_HG19')]
							#instance.udg = [library_info[library_id][library_headers.index('UDG_Treatment')]]
							instance.x_contam_point = library_info[library_id][library_headers.index('Nuclear_ANGSD_Mean')]
							instance.x_contam_z = library_info[library_id][library_headers.index('Nuclear_ANGSD_Z')]
							instance.sex = library_info[library_id][library_headers.index('Nuclear_Sex')]
							
							coverage_string = library_info[library_id][library_headers.index('Nuclear_Coverage_Targeted_Positions')]
							if coverage_string != '':
								instance.coverage = float(coverage_string)
							autosome_snps_string = library_info[library_id][library_headers.index('Nuclear_Unique_SNPS_Hit')]
							if autosome_snps_string != '':
								instance.autosome_snps = int(autosome_snps_string)
							
							instance.endogenous_by_library = [library_info[library_id][library_headers.index('Shotgun_Percent_HG19')]]
							instance.damage_by_library = [library_info[library_id][library_headers.index('Nuclear_Damage_Last_Base')]]
							instance.mt_haplogroup_by_library = [library_info[library_id][library_headers.index('mtDNA_Haplogroup')]]
							instance.mt_contamination_by_library = [library_info[library_id][library_headers.index('mtDNA_Consensus_Match')]]
							
							logfile, bamfile = library_id_in_pulldown(library_id, args.names, args.pulldown_dir_root)
							instance.pulldown_logfile = logfile
							instance.pulldown_1_sample_id = library_id
							instance.pulldown_3_bam = str(bamfile)
							
							sample_id = 'S{:d}'.format(LibraryID(library_id).sample)
							instance.sample_id = sample_id
							
							instance_list[library_id] = instance

	# Read in merge lists. Each merges gets an anno entry
	if args.merge_lists is not None:
		for merge_list in args.merge_lists:
			with open(merge_list) as f:
				for line in f:
					fields = re.split('\t|\n', line.rstrip())
					instance_id = fields[0]
					individual_id = fields[1]
					libraries = fields[2:]
					if '' in libraries:
						print(line, file=sys.stderr)
						raise ValueError('empty library')
					
					instance = InstanceAnnoEntry(instance_id)
					instance.library_ids = libraries
					if instance_id in instance_list:
						raise ValueError('duplicate instance: {}'.format(instance_id))
					else:
						instance_list[instance_id] = instance
						instance.library_ids = libraries
						
						#instance.pulldown_logfile # TODO external
						instance.pulldown_1_sample_id = instance_id
						individual_id = instance_id.split('_')[0]
						instance.sample_id = individual_id.replace('I', 'S')
						instance.pulldown_3_bam = '/n/groups/reich/matt/pipeline/sample_merge/{}/{}.hg19.bam'.format(individual_id, instance_id)
						exist_check = Path(instance.pulldown_3_bam)
						if not exist_check.exists():
							print('{} does not exist'.format(exist_check), file=sys.stderr)
						
						instance.mtdna_haplogroup = merge_analysis_info[instance_id][merge_analysis_headers.index('MT_Haplogroup')]
						mtdna_coverage = merge_analysis_info[instance_id][merge_analysis_headers.index('MT_coverage')]
						if mtdna_coverage != '':
							instance.mtdna_coverage = float(mtdna_coverage)
							
						instance.udg = [library_info[library_id][library_headers.index('UDG_Treatment')] for library_id in libraries]
						
						instance.x_contam_point = merge_analysis_info[instance_id][merge_analysis_headers.index('angsd_MoM')]
						instance.x_contam_z = merge_analysis_info[instance_id][merge_analysis_headers.index('angsd_SE(MoM)')]
						instance.sex = merge_analysis_info[instance_id][merge_analysis_headers.index('1240k_post_sex')]
						
						#autosome_snps_string = library_info[library_id][library_headers.index('Nuclear_Unique_SNPS_Hit')] # TODO get from logs external
						#if autosome_snps_string != '':
							#instance.autosome_snps = int(autosome_snps_string)
						
						instance.endogenous_by_library = [library_info[library_id][library_headers.index('Shotgun_Percent_HG19')] for library_id in libraries]
						instance.damage_by_library = [library_info[library_id][library_headers.index('Nuclear_Damage_Last_Base')] for library_id in libraries]
						instance.mt_haplogroup_by_library = [library_info[library_id][library_headers.index('mtDNA_Haplogroup')] for library_id in libraries]
						instance.mt_contamination_by_library = [library_info[library_id][library_headers.index('mtDNA_Consensus_Match')] for library_id in libraries]
						
						try:
							instance.best_shotgun_endogenous = max([float(x) for x in instance.endogenous_by_library if x != ''])
						except:
							print(instance.endogenous_by_library, file=sys.stderr)
	
	for instance_id, instance in instance_list.items():
		sample_id = instance.sample_id # lookup sample based on instance_id, latest sample number?
		try:
			sample = sample_info[sample_id]
			# fill in sample fields
			instance.skeletal_code = sample[sample_headers.index('Skeletal Code')]
			instance.skeletal_element = sample[sample_headers.index('Skeletal Element')]
			instance.collaborator = sample[sample_headers.index('Collaborator')]
			instance.date_formatted_check = sample[sample_headers.index('Date Fix Flag')]
			instance.date = sample[sample_headers.index('Average BP Date')]
			instance.date_notes = sample[sample_headers.index('Date')]
			instance.group_label = sample[sample_headers.index('Population Label')]
			instance.location = sample[sample_headers.index('Locality')]
			instance.country = sample[sample_headers.index('Country')]
			instance.latitude = sample[sample_headers.index('Latitude')]
			instance.longitude = sample[sample_headers.index('Longitude')]
			
			instance.carbon14_data_status = sample[sample_headers.index('14C Status')].strip()
			#instance.year_published =
		except Exception as exception:
			print(exception, file=sys.stderr)
		print(instance)
