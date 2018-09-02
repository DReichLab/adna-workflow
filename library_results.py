import sys
import os
from pathlib import Path
import argparse
from library_id import LibraryID

# Map headers in Rebecca's sample spreadsheet to the library spreadsheet
header_mapping_sample = {
	'sample_id' : 'Sample_ID',
	'individual_id' : 'individual_id',
	'shipment_id' : 'shipment_id',
	'Skeletal codes' : 'Skeletal codes',
	'Skeletal element' : 'Skeletal element',
	'Collaborator' : 'Collaborator',
	'Date Fix Flag' : 'Date Fix Flag',
	'  Average of 95.4% date range in calBP (defined as 1950 CE)' : '  Average of 95.4% date range in calBP (defined as 1950 CE)',
	'Date' : 'Date: One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age, Lab number) e.g. 5983-5747 calBCE (6980±50 BP, Beta-226472). (Format 2) Archaeological context date, e.g. 2500-1700 BCE',
	'Culture' : 'Culture',
	'Location' : 'Location',
	'Country' : 'Country',
	'  Lat.' : '  Lat.',
	'  Long.' : '  Long.',
	'notes' : 'Notes',
	'notes_2' : 'notes_2',
	'14C_Status' : '14C_Status',
	'AMS_Ship_ID' : '14C Ship-ID',
	'AMS_Ship_Date' : '14C Ship Date',
	'sampling_tech' : 'Sampling_Tech'
}

def read_library_file(filename):
	with open(filename) as f:
		# read header line
		header_line = f.readline()
		headers = header_line.split('\t')
		library_ids = []
		library_info = {}
		id_index = headers.index('Library_id')
		
		for line in f:
			fields = line.split('\t')
			library_id = fields[id_index]
			
			if len(fields) != len(headers):
				raise ValueError('mismatch between headers and fields')
			
			if library_id in library_info:
				raise ValueError('{} appears more than once in library file'.format(library_id))
			else:
				library_ids.append(library_id)
				library_info[library_id] = fields
	return headers, library_ids, library_info

def read_sample_file(filename):
	sample_info = {}
	with open(filename) as f:
		# read header line
		header_line = f.readline()
		headers = header_line.split('\t')
		sample_id_index = headers.index('sample_id')
		
		for line in f:
			fields = line.split('\t')
			sample_id = fields[sample_id_index]
			sample_info[sample_id] = fields
			
	return headers, sample_info

def replace_field(current_library, library_headers, report_entry_values, report_headers, library_header_to_replace, report_field_replacing):
	value = report_entry_values[report_headers.index(report_field_replacing)]
	to_replace_index = library_headers.index(library_header_to_replace)
	if value != '..':
		current_library[to_replace_index] = str(value)
		
def replace_capture_fields(current_library, library_headers, report_entry_values, report_headers):
	header_mapping_1240k = {
		'mt-experiment raw sequences' : 'raw',
		'mt-experiment sequences going into alignment (passing barcode and other filters)' : 'merged',
		'mt sequences aligning to target' : 'MT_pre',
		'mt sequences aligning to target post-dup' : 'MT_post',
		'mt coverage' : 'MT_post-coverage',
		'mt mean length (median for old-style analyses)' : 'mean_rsrs',
		'mt match to consensus (only if mtcov>2)' : 'contamination_contammix',
		'mt haplogroup (only if mtcov>2)' : 'MT_Haplogroup',
		'mt haplogroup confidence (only if mtcov>2)' : 'MT_Haplogroup_rank',

		' 1240K sequences de-indexing for Matt (raw reads for Shop) ' : 'raw',
		' 1240K sequences that merge and pass barcode check ' : 'merged',
		'1240K coverage on targeted positions (from reporting)' : '1240k_post_autosome',
		'1240K Expected Coverage at 10% Marginal Uniqueness' : 'preseq_coverage_at_marginal_uniqueness_0.10',
		'1240K Expected Coverage at 37% Marginal Uniqueness' : 'preseq_coverage_at_marginal_uniqueness_0.368',
		'1240K marginal uniqueness' : 'preseq_marginal_uniqueness',
		'1240K mean sequence length (median for Shop processing)' : 'mean_nuclear',
		'1240K X hits' : '1240k_post_x',
		'1240K Y hits' : '1240k_post_y',
		'1240K Sex' : '1240k_post_sex',
		'1240K ANGSD-SNPs' : 'angsd_nsites',
		'1240K ANGSD-mean' : 'angsd_MoM'
	}

	for library_header_to_replace, report_field_replacing in header_mapping_1240k.items():
		replace_field(current_library, library_headers, report_entry_values, report_headers, library_header_to_replace, report_field_replacing)
	
	try:
		damage_rsrs_ct1 = float(report_entry_values[report_headers.index('damage_rsrs_ct1')])
		damage_rsrs_ga1 = float(report_entry_values[report_headers.index('damage_rsrs_ga1')])
		current_library[library_headers.index('mt fraction damaged in last base')] = '{:.3f}'.format((damage_rsrs_ct1 + damage_rsrs_ga1) / 2)
	except:
		pass
	
	try:
		contamination_contammix_lower = float(report_entry_values[report_headers.index('contamination_contammix_lower')])
		contamination_contammix_upper = float(report_entry_values[report_headers.index('contamination_contammix_upper')])
		current_library[library_headers.index('mt match to consensus 95CI (only if mtcov>2)')] = '[{:.3f}, {:.3f}]'.format(contamination_contammix_lower, contamination_contammix_upper)
	except:
		pass
	
	try:
		z1240k_pre_x = int(report_entry_values[report_headers.index('1240k_pre_x')])
		z1240k_pre_y = int(report_entry_values[report_headers.index('1240k_pre_y')])
		z1240k_pre_autosome = int(report_entry_values[report_headers.index('1240k_pre_autosome')])
		current_library[library_headers.index('1240K sequences passing QC hitting targets pre-dedup')] = '{:d}'.format(z1240k_pre_x + z1240k_pre_y + z1240k_pre_autosome)
	except:
		pass
	
	try:
		z1240k_post_x = int(report_entry_values[report_headers.index('1240k_post_x')])
		z1240k_post_y = int(report_entry_values[report_headers.index('1240k_post_y')])
		z1240k_post_autosome = int(report_entry_values[report_headers.index('1240k_post_autosome')])
		current_library[library_headers.index('1240K sequences passing QC hitting targets post-dedup')] = '{:d}'.format(z1240k_post_x + z1240k_post_y + z1240k_post_autosome)
	except:
		pass
	
	try:
		damage_nuclear_ct1 = float(report_entry_values[report_headers.index('damage_nuclear_ct1')])
		damage_nuclear_ga1 = float(report_entry_values[report_headers.index('damage_nuclear_ga1')])
		current_library[library_headers.index('1240K damage in last base')] = '{:.3f}'.format((damage_nuclear_ct1 + damage_nuclear_ga1) / 2)
	except:
		pass
	
	try:
		angsd_MoM = float(report_entry_values[report_headers.index('angsd_MoM')])
		angsd_SE = float(report_entry_values[report_headers.index('angsd_SE(MoM)')])
		current_library[library_headers.index('1240K ANGSD-Z')] = '{:.3f}'.format(angsd_MoM / angsd_SE)
	except:
		pass
	'1240K Y chromsome haplogroup'

def replace_shotgun_fields(current_library, library_headers, report_entry_values, report_headers):
	header_mapping_shotgun = {
		'Shotgun Raw Sequences' : 'raw',
		'Shotgun sequences going into alignment (passing barcode and other filters)' : 'merged',
		"Shotgun mean length (median length for Shop's processing)" : 'mean_nuclear',
		'Shotgun Fraction Raw Sequences Mapping to hg19 (% Human)' : 'endogenous_pre',
	}
	
	for library_header_to_replace, report_field_replacing in header_mapping_shotgun.items():
		replace_field(current_library, library_headers, report_entry_values, report_headers, library_header_to_replace, report_field_replacing)

	try:
		autosome_pre = int(report_entry_values[report_headers.index('autosome_pre')])
		X_pre = int(report_entry_values[report_headers.index('X_pre')])
		Y_pre = int(report_entry_values[report_headers.index('Y_pre')])
		MT_pre = int(report_entry_values[report_headers.index('MT_pre')])
		current_library[library_headers.index('Shotgun reads mapping to hg19')] = '{:d}'.format(autosome_pre + X_pre + Y_pre + MT_pre)
	except:
		pass
	
	try:
		damage_nuclear_ct1 = float(report_entry_values[report_headers.index('damage_nuclear_ct1')])
		damage_nuclear_ga1 = float(report_entry_values[report_headers.index('damage_nuclear_ga1')])
		current_library[library_headers.index('Shotgun damage rate')] = '{:.3f}'.format((damage_nuclear_ct1 + damage_nuclear_ga1) / 2)
	except:
		pass
	
	try:
		current_library[library_headers.index('Shotgun Fraction Aligning to hg19 that Also Hit mtDNA for ones with at least 1,000 shotgun reads matchin hg19')] = '{:.3f}'.format(MT_pre / (autosome_pre + X_pre + Y_pre + MT_pre))
	except:
		pass

# combine Rebecca's library data with results of pipeline analysis
def read_pipeline_analysis_report(pipeline_report_filename, library_headers, library_info, sample_headers, sample_info):
	with open(pipeline_report_filename) as f:
		f.readline() # first line is read count
		header_line = f.readline() # second line is header fields
		headers = header_line.split('\t')
		report_library_id_index = headers.index('library_id')
		experiment_index = headers.index('experiment')
		
		# each line is one library
		# iterate through report libraries and update corresponding library info
		for line in f:
			fields = line.split('\t')
			
			library_id = fields[report_library_id_index]
			experiment = fields[experiment_index]
			if library_id.startswith('S'): # is not '' and library_id is not 'Contl.Capture':
				current_library = library_info[library_id]
				# sample file data
				try:
					sample_id = LibraryID(library_id).sample
					for sample_header in header_mapping_sample:
						index = sample_headers.index(sample_header)
						value = sample_info['S' + str(sample_id)][index]
						if value != '..':
							mapped_header = header_mapping_sample[sample_header]
							#print(mapped_header, library_headers.index(mapped_header), value)
							current_library[library_headers.index(mapped_header)] = value
				except Exception as exception:
					print(exception, file=sys.stderr)
					#raise
				
				if len(fields) == len(headers): # no data will have fewer fields than headers
					if '1240k' in experiment:
						replace_capture_fields(current_library, library_headers, fields, headers)
					elif 'Raw' in experiment:
						replace_shotgun_fields(current_library, library_headers, fields, headers)
						
# read a pulldown log file and return a map of ids to SNPs
def pulldown_snp_stats(filename):
	library_targets = {}
	with open(filename) as f:
		for line in f:
			if 'coverage' in line:
				fields = line.split()
				library_id = fields[0]
				targets = int(fields[-1])
				library_targets[library_id] = targets
	return library_targets

# run this from release directory where pulldown directories are subdirectories
def logfile_and_dblist(name, library_headers, library_info):
	parent_path = Path(name)
	#logfile_name_only = '{}.half.normal.parameters.stdout'.format(name)
	#print(logfile_name_only, file=sys.stderr)
	logfile_fullpath = str(parent_path.resolve()) + '/{}.half.normal.parameters.stdout'.format(name)
	#print(logfile_fullpath, file=sys.stderr)
	library_targets = pulldown_snp_stats(logfile_fullpath)
	
	logfile_index = library_headers.index('pulldown: logfile location')
	
	dblist_filename = parent_path / ('{}.dblist'.format(name))
	#print(dblist_filename, file=sys.stderr)
	with open(dblist_filename) as dblist:
		for line in dblist:
			fields = line.split()
			if len(fields) == 4:
				library_id = fields[0]
				second = fields[1]
				bam = fields[2]
				read_groups = fields[3]
				
				if library_id in library_info:
					current_library = library_info[library_id]
					
					current_library[logfile_index] = logfile_fullpath
					current_library[library_headers.index('pulldown: 1st column of nickdb (sample ID used in logfile)')] = library_id
					current_library[library_headers.index('pulldown: 2nd column of nickdb (alt sample ID)')] = library_id
					current_library[library_headers.index('pulldown: 3rd column of nickdb (bam)')] = bam
					#'pulldown: 4th column of nickdb (or hetfa)'
					current_library[library_headers.index('pulldown: 5th column of nickdb (readgroup list or diploid source)')] = read_groups
					
					try:
						targets = library_targets[library_id]
						current_library[library_headers.index('1240K number of unique autosomal targets hit (from pulldown)')] = '{:d}'.format(targets)
					except:
						print('{} not in library_targets'.format(library_id), file=sys.stderr)
					#current_library[library_headers.index('1240K fraction of unique targets hit')] = '{:.3f}'.format(targets / )
				else:
					print('{} not in libraries'.format(library_id), file=sys.stderr)
				

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Combine pipeline analysis results with Rebecca's library file.")
	parser.add_argument('-l', "--library_file", help="Library file describing each library from Rebecca", required=True)
	parser.add_argument('-s', "--sample_file", help="Sample file describing samples", required=True)
	parser.add_argument('-r', "--reports", help="Pipeline report files", nargs='*')
	parser.add_argument('-d', "--directory_pulldown", help="Directory containing logfile and dblist", nargs='*')
	args = parser.parse_args()
	
	# read library info into memory
	library_headers, library_ids, library_info = read_library_file(args.library_file)
	# read sample info into memory
	sample_headers, sample_info = read_sample_file(args.sample_file)
	# add pipeline analysis reports into library entries
	reports = args.reports if args.reports is not None else []
	for report_filename in reports:
		read_pipeline_analysis_report(report_filename, library_headers, library_info, sample_headers, sample_info)
		
	pulldown_directories = args.directory_pulldown if args.directory_pulldown is not None else []
	for pulldown_directory in pulldown_directories:
		logfile_and_dblist(pulldown_directory, library_headers, library_info)
	# print back out library info
	# print header
	print('\t'.join(library_headers), end='')
	# print library data
	for library_id in library_ids:
		print('\t'.join(library_info[library_id]), end='')
	
