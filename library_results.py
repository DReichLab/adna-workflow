import sys
import os
from pathlib import Path
import argparse
from library_id import LibraryID

# Map headers in Rebecca's sample spreadsheet to the library spreadsheet
# key is sample spreadsheet column name, value is library spreadsheet column name
header_mapping_sample = {
	'Sample-ID' : 'Sample_ID',
	'Individual_ID' : 'Individual_ID',
	'Shipment_ID' : 'Shipment_ID',
	'Skeletal_Code' : 'Skeletal_Code',
	'Skeletal_Element' : 'Skeletal_Element',
	'Collaborator' : 'Collaborator',
	'Date_Fix_Flag' : 'Date_Fix_Flag',
	'Average_BP_Date' : 'Average_BP_Date',
	'Date' : 'Date',
	'Culture_Period' : 'Culture_Period',
	'Locality' : 'Locality',
	'Country' : 'Country',
	'Latitude' : 'Latitude',
	'Longitude' : 'Longitude',
	'Notes' : 'Notes',
	'Notes_2' : 'Notes_2',
	#'14C_Status' : '14C_Status',
	#'AMS_Ship_ID' : '14C Ship-ID',
	#'AMS_Ship_Date' : '14C Ship Date',
	'Sampling_Techique' : 'Sampling_Tech'
}

def read_library_file(filename):
	with open(filename) as f:
		# read header line
		header_line = f.readline()
		headers = header_line.split('\t')
		library_ids = []
		library_info = {}
		id_index = headers.index('Library_ID')
		
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

def read_sample_file(filename, delimiter):
	sample_info = {}
	with open(filename) as f:
		# read header line
		header_line = f.readline()
		headers = header_line.split(delimiter)
		sample_id_index = headers.index('Sample-ID')
		
		for line in f:
			fields = line.split(delimiter)
			sample_id = fields[sample_id_index]
			sample_info[sample_id] = fields
			
	return headers, sample_info

def replace_field(current_library, library_headers, report_entry_values, report_headers, library_header_to_replace, report_field_replacing):
	value = report_entry_values[report_headers.index(report_field_replacing)]
	to_replace_index = library_headers.index(library_header_to_replace)
	if value != '..':
		current_library[to_replace_index] = str(value)
		
def replace_capture_fields(current_library, library_headers, report_entry_values, report_headers):
	# key is column name in Rebecca's library spreadsheet, value is field name from pipeline report
	header_mapping_1240k = {
		'mtDNA_Raw_Sequences' : 'raw',
		'mtDNA_Sequences_Passing_Filters' : 'merged',
		'mtDNA_Sequences_Target_Alignment' : 'MT_pre',
		'mtDNA_Sequences_Target_Alignment_PostDedup' : 'MT_post',
		'mtDNA_Coverage' : 'MT_post-coverage',
		'mtDNA_Mean_Median_Seq_Length' : 'mean_rsrs',
		'mtDNA_Consensus_Match' : 'contamination_contammix',
		'mtDNA_Haplogroup' : 'MT_Haplogroup',
		'mtDNA_Haplogroup_Confidence' : 'MT_Haplogroup_rank',

		'Nuclear_Raw_Reads_OR_Deindexing' : 'raw',
		'Nuclear_Sequences_Merge_Pass_Barcode' : 'merged',
		'Nuclear_Expected_Coverage_10%_Marginal_Uniqueness' : 'preseq_coverage_at_marginal_uniqueness_0.10',
		'Nuclear_Expected_Coverage_37%_Marginal_Uniqueness' : 'preseq_coverage_at_marginal_uniqueness_0.368',
		'Nuclear_Marginal_Uniqueness' : 'preseq_marginal_uniqueness',
		'Nuclear_Mean_Median_Seq_Length' : 'mean_nuclear',
		'Nuclear_X_Hits' : '1240k_post_x',
		'Nuclear_Y_Hits' : '1240k_post_y',
		'Nuclear_Sex' : '1240k_post_sex',
		'Nuclear_ANGSD_SNPs' : 'angsd_nsites',
		'Nuclear_ANGSD_Mean' : 'angsd_MoM',
		#' 1240K ANGSD-Z ' : 'angsd_MoM_z'
	}

	for library_header_to_replace, report_field_replacing in header_mapping_1240k.items():
		replace_field(current_library, library_headers, report_entry_values, report_headers, library_header_to_replace, report_field_replacing)
	
	try:
		damage_rsrs_ct1 = float(report_entry_values[report_headers.index('damage_rsrs_ct1')])
		damage_rsrs_ga1 = float(report_entry_values[report_headers.index('damage_rsrs_ga1')])
		current_library[library_headers.index('mtDNA_Damage_Last_Base')] = '{:.3f}'.format((damage_rsrs_ct1 + damage_rsrs_ga1) / 2)
	except:
		pass
	
	try:
		contamination_contammix_lower = float(report_entry_values[report_headers.index('contamination_contammix_lower')])
		contamination_contammix_upper = float(report_entry_values[report_headers.index('contamination_contammix_upper')])
		current_library[library_headers.index('mtDNA_Consensus_Match_95CI')] = '[{:.3f}, {:.3f}]'.format(contamination_contammix_lower, contamination_contammix_upper)
	except:
		pass
	
	try:
		z1240k_pre_x = int(report_entry_values[report_headers.index('1240k_pre_x')])
		z1240k_pre_y = int(report_entry_values[report_headers.index('1240k_pre_y')])
		z1240k_pre_autosome = int(report_entry_values[report_headers.index('1240k_pre_autosome')])
		current_library[library_headers.index('Nuclear_Target_Sequences_Pass_QC_PreDedup')] = '{:d}'.format(z1240k_pre_x + z1240k_pre_y + z1240k_pre_autosome)
	except:
		pass
	
	try:
		z1240k_post_x = int(report_entry_values[report_headers.index('1240k_post_x')])
		z1240k_post_y = int(report_entry_values[report_headers.index('1240k_post_y')])
		z1240k_post_autosome = int(report_entry_values[report_headers.index('1240k_post_autosome')])
		current_library[library_headers.index('Nuclear_Target_Sequences_Pass_QC_PostDedup')] = '{:d}'.format(z1240k_post_x + z1240k_post_y + z1240k_post_autosome)
	except:
		pass
	
	try:
		num_1240k_autosome_targets = 1150639
		current_library[library_headers.index('Nuclear_Coverage_Targeted_Positions')] = '{:.3f}'.format(z1240k_post_autosome / num_1240k_autosome_targets)
	except:
		pass
	
	try:
		damage_nuclear_ct1 = float(report_entry_values[report_headers.index('damage_nuclear_ct1')])
		damage_nuclear_ga1 = float(report_entry_values[report_headers.index('damage_nuclear_ga1')])
		current_library[library_headers.index('Nuclear_Damage_Last_Base')] = '{:.3f}'.format((damage_nuclear_ct1 + damage_nuclear_ga1) / 2)
	except:
		pass
	
	try:
		angsd_MoM = float(report_entry_values[report_headers.index('angsd_MoM')])
		angsd_SE = float(report_entry_values[report_headers.index('angsd_SE(MoM)')])
		current_library[library_headers.index('Nuclear_ANGSD_Z')] = '{:.3f}'.format(angsd_MoM / angsd_SE)
	except:
		pass
	'1240K Y chromsome haplogroup'

def replace_shotgun_fields(current_library, library_headers, report_entry_values, report_headers):
	header_mapping_shotgun = {
		'Shotgun_Raw_Sequences' : 'raw',
		'Shotgun_Sequences_Passing_Filters' : 'merged',
		"Shotgun_Mean_Median_Seq_Length" : 'mean_nuclear',
		'Shotgun_Percent_HG19' : 'endogenous_pre',
	}
	
	for library_header_to_replace, report_field_replacing in header_mapping_shotgun.items():
		replace_field(current_library, library_headers, report_entry_values, report_headers, library_header_to_replace, report_field_replacing)

	try:
		autosome_pre = int(report_entry_values[report_headers.index('autosome_pre')])
		X_pre = int(report_entry_values[report_headers.index('X_pre')])
		Y_pre = int(report_entry_values[report_headers.index('Y_pre')])
		MT_pre = int(report_entry_values[report_headers.index('MT_pre')])
		current_library[library_headers.index('Shotgun_Reads_Mapped_HG19')] = '{:d}'.format(autosome_pre + X_pre + Y_pre + MT_pre)
	except:
		pass
	
	try:
		damage_nuclear_ct1 = float(report_entry_values[report_headers.index('damage_nuclear_ct1')])
		damage_nuclear_ga1 = float(report_entry_values[report_headers.index('damage_nuclear_ga1')])
		current_library[library_headers.index('Shotgun_Damage_Rate')] = '{:.3f}'.format((damage_nuclear_ct1 + damage_nuclear_ga1) / 2)
	except:
		pass
	
	try:
		current_library[library_headers.index('Shotgun_Fraction_HG19_Hit_mtDNA')] = '{:.3f}'.format(MT_pre / (autosome_pre + X_pre + Y_pre + MT_pre))
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
				'''
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
				'''
				
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
	
	logfiles = ['{}.half.normal.parameters.stdout'.format(name), '{}.minus.normal.parameters.stdout'.format(name)]
	for logfile in logfiles:
		logfile_fullpath = parent_path.resolve() / logfile
		print(logfile_fullpath, file=sys.stderr)
		if logfile_fullpath.is_file() and  logfile_fullpath.exists():
			library_targets = pulldown_snp_stats(logfile_fullpath)
			
			logfile_index = library_headers.index('Pulldown_Logfile_Location')
			
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
							
							current_library[logfile_index] = str(logfile_fullpath)
							current_library[library_headers.index('Pulldown_1st_Column_NickDB')] = library_id
							current_library[library_headers.index('Pulldown_2nd_Column_NickDB_alt_sample')] = library_id
							current_library[library_headers.index('Pulldown_3rd_Column_NickDB_bam')] = bam
							#'pulldown: 4th column of nickdb (or hetfa)'
							current_library[library_headers.index('Pulldown_5th_Column_NickDB_readgroup_diploid_source')] = read_groups
							
							try:
								targets = library_targets[library_id]
								current_library[library_headers.index('Nuclear_Unique_SNPS_Hit')] = '{:d}'.format(targets)
							except:
								print('{} not in library_targets'.format(library_id), file=sys.stderr)
							#current_library[library_headers.index('1240K fraction of unique targets hit')] = '{:.3f}'.format(targets / )
						else:
							print('{} not in libraries'.format(library_id), file=sys.stderr)
				

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Combine pipeline analysis results with Rebecca's library file.")
	parser.add_argument('-l', "--library_file", help="Library file describing each library from Rebecca", required=True)
	parser.add_argument('-s', "--sample_file", help="Sample file describing samples", required=True)
	parser.add_argument('-x', "--delimiter", help="Sample file field delimiter", default='\t')
	parser.add_argument('-r', "--reports", help="Pipeline report files", nargs='*')
	parser.add_argument('-d', "--directory_pulldown", help="Directory containing logfile and dblist", nargs='*')
	args = parser.parse_args()
	
	# read library info into memory
	library_headers, library_ids, library_info = read_library_file(args.library_file)
	# read sample info into memory
	sample_headers, sample_info = read_sample_file(args.sample_file, args.delimiter)
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
	
