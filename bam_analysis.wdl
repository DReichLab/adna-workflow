import "demultiplex.wdl" as demultiplex_align_bams
import "analysis.wdl" as analysis
import "analysis_clipping.wdl" as analysis_clipping

workflow bam_analysis{
	String label
	String genome_reference_string
	String mt_reference_string
	Array[File] input_nuclear_bams
	Array[File] input_mt_bams

	File adna_screen_jar
	File picard_jar
	File htsbox
	File pmdtools
	File haplogrep_jar
	
	File mt_reference_rsrs_in
	File mt_reference_rcrs_in

	File python_damage_two_bases
	File python_angsd_results
	File python_target
	File python_read_groups_from_bam
	File python_central_measures

	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	Int minimum_mapping_quality
	Int minimum_base_quality
	
	Int deamination_bases_to_clip_half
	Int deamination_bases_to_clip_minus
	Int deamination_bases_to_clip_plus
	File udg_minus_libraries_file
	File udg_plus_libraries_file
	Array[String] udg_minus_libraries = read_lines(udg_minus_libraries_file)
	Array[String] udg_plus_libraries = read_lines(udg_plus_libraries_file)
	
	call demultiplex_align_bams.prepare_reference as prepare_reference_rsrs{ input:
		reference = mt_reference_rsrs_in
	}
	call demultiplex_align_bams.prepare_reference as prepare_reference_rcrs{ input:
		reference = mt_reference_rcrs_in
	}
	
	call remove_marked_duplicates as remove_marked_duplicates_nuclear{ input:
		bams = input_nuclear_bams,
		references = [genome_reference_string, mt_reference_string],
	}
	call analysis_clipping.clip_deamination as clip_nuclear { input:
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		bams = remove_marked_duplicates_nuclear.no_duplicates_bams,
		deamination_bases_to_clip_half = deamination_bases_to_clip_half,
		deamination_bases_to_clip_minus = deamination_bases_to_clip_minus,
		deamination_bases_to_clip_plus = deamination_bases_to_clip_plus,
		udg_minus_libraries = udg_minus_libraries,
		udg_plus_libraries = udg_plus_libraries,
		python_read_groups_from_bam = python_read_groups_from_bam
	}
	
	call remove_marked_duplicates as remove_marked_duplicates_mt{ input:
		bams = input_mt_bams,
		references = [genome_reference_string, mt_reference_string],
	}
	call analysis_clipping.clip_deamination as clip_mt { input:
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		bams = remove_marked_duplicates_mt.no_duplicates_bams,
		deamination_bases_to_clip_half = deamination_bases_to_clip_half,
		deamination_bases_to_clip_minus = deamination_bases_to_clip_minus,
		deamination_bases_to_clip_plus = deamination_bases_to_clip_plus,
		udg_minus_libraries = udg_minus_libraries,
		udg_plus_libraries = udg_plus_libraries,
		python_read_groups_from_bam = python_read_groups_from_bam
	}
	call analysis_clipping.clip_deamination as clip_mt_hard { input:
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		bams = remove_marked_duplicates_mt.no_duplicates_bams,
		deamination_bases_to_clip_half = deamination_bases_to_clip_half,
		deamination_bases_to_clip_minus = deamination_bases_to_clip_minus,
		deamination_bases_to_clip_plus = deamination_bases_to_clip_plus,
		udg_minus_libraries = udg_minus_libraries,
		udg_plus_libraries = udg_plus_libraries,
		python_read_groups_from_bam = python_read_groups_from_bam,
		hard_clip = true
	}
	
	call analysis.damage_loop as damage_nuclear{ input :
		pmdtools = pmdtools,
		python_damage_two_bases = python_damage_two_bases,
		bams = remove_marked_duplicates_nuclear.no_duplicates_bams,
		damage_label = "damage_nuclear",
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		processes = 12
	}
	call analysis.damage_loop as damage_mt{ input :
		pmdtools = pmdtools,
		python_damage_two_bases = python_damage_two_bases,
		bams = remove_marked_duplicates_mt.no_duplicates_bams,
		damage_label = "damage_mt",
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		processes = 4
	}
	call analysis_clipping.angsd_contamination{ input:
		bams = clip_nuclear.clipped_bams,
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		python_angsd_results = python_angsd_results,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		processes = 10
	}
	call analysis_clipping.haplogrep as haplogrep_rcrs{ input:
		missing_alignments_fraction = missing_alignments_fraction,
		max_open_gaps = max_open_gaps,
		seed_length = seed_length,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		bams = clip_mt_hard.clipped_bams,
		reference = prepare_reference_rcrs.reference_fa,
		reference_amb = prepare_reference_rcrs.reference_amb,
		reference_ann = prepare_reference_rcrs.reference_ann,
		reference_bwt = prepare_reference_rcrs.reference_bwt,
		reference_pac = prepare_reference_rcrs.reference_pac,
		reference_sa = prepare_reference_rcrs.reference_sa,
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		haplogrep_jar = haplogrep_jar,
		processes = 3
	}
	call analysis.summarize_haplogroups{ input:
		haplogrep_output = haplogrep_rcrs.haplogroup_report
	}
	call analysis_clipping.chromosome_target as nuclear_chromosome_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = clip_nuclear.clipped_bams,
		targets="\"{'autosome_post':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X_post':'X','Y_post':'Y','human_post':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	call analysis_clipping.chromosome_target as rsrs_chromosome_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = clip_mt.clipped_bams,
		targets="\"{'MT_post':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	call coverage_without_index_barcode_key as rsrs_coverage{ input:
		bam_stats = rsrs_chromosome_target_post.target_stats,
		reference_length = 16569,
		coverage_field = "MT_post-coverageLength"
	}
	scatter(bam in clip_mt_hard.clipped_bams){
		call analysis_clipping.contammix{ input:
			bam = bam,
			picard_jar = picard_jar,
			htsbox = htsbox,
			missing_alignments_fraction = missing_alignments_fraction,
			max_open_gaps = max_open_gaps,
			seed_length = seed_length,
			minimum_mapping_quality = minimum_mapping_quality,
			minimum_base_quality = minimum_base_quality,
			reference = prepare_reference_rsrs.reference_fa,
			reference_amb = prepare_reference_rsrs.reference_amb,
			reference_ann = prepare_reference_rsrs.reference_ann,
			reference_bwt = prepare_reference_rsrs.reference_bwt,
			reference_pac = prepare_reference_rsrs.reference_pac,
			reference_sa = prepare_reference_rsrs.reference_sa,
			reference_fai = prepare_reference_rsrs.reference_fai,
			coverages = rsrs_coverage.coverages
		}
	}
	call analysis.concatenate as concatenate_contammix{ input:
		to_concatenate = contammix.contamination_estimate
	}
	call analysis_clipping.snp_target_bed as count_1240k_post { input:
#		coordinates_autosome = coordinates_1240k_autosome,
#		coordinates_x = coordinates_1240k_x,
#		coordinates_y = coordinates_1240k_y,
		bams = clip_nuclear.clipped_bams,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		label = "1240k_post",
#		python_snp_target_bed = python_snp_target_bed,
		processes = 12
	}
	call analysis.concatenate as concatenate_count_1240k_post{ input:
		to_concatenate = count_1240k_post.snp_target_stats
	}
	
	call analysis_clipping.central_measures as central_measures_nuclear{ input:
		python_central_measures = python_central_measures,
		mean_label = "mean_nuclear",
		median_label = "median_nuclear",
		histograms = nuclear_chromosome_target_post.length_histogram
	}
	call analysis_clipping.central_measures as central_measures_rsrs{ input:
		python_central_measures = python_central_measures,
		mean_label = "mean_rsrs",
		median_label = "median_rsrs",
		histograms = rsrs_chromosome_target_post.length_histogram
	}
	
	call analysis_results{ input:
		keyed_value_results = [
			damage_nuclear.damage_all_samples_two_bases,
			damage_mt.damage_all_samples_two_bases,
			central_measures_nuclear.central_measures_output,
			central_measures_rsrs.central_measures_output,
			angsd_contamination.contamination,
			concatenate_count_1240k_post.concatenated,
			rsrs_coverage.coverage_statistics,
			concatenate_contammix.concatenated,
			summarize_haplogroups.haplogroups
		]
	}
	
	output{
		Array[File] clipped_mt_bams = clip_mt.clipped_bams
		Array[File] mt_bams = remove_marked_duplicates_mt.no_duplicates_bams
		Array[File] clipped_nuclear_bams = clip_nuclear.clipped_bams
		Array[File] nuclear_bams = remove_marked_duplicates_nuclear.no_duplicates_bams
	}
}

# take a list of instance ids with nuclear and MT bams
# rename them using symbolic links so the analysis returns intelligible results using the filenames
task prepare_bam_links{
	# each line looks like:
	# [id] [nuclear bam] [mt bam]library ids
	# does each line also need version number?
	File bams_to_analyze
	
	command{
		set -e
		mkdir -p nuclear
		mkdir -p mt
		python3 <<CODE
			import os
			import re
			
			with open("${bams_to_analyze}", 'r') as f:
				for line in f:
					fields = re.split('\t|\n', line)
					instance_id = fields[0]
					nuclear_bam = fields[1]
					mt_bam = fields[2]
					
					os.symlink(nuclear_bam, 'nuclear/' + instance_id + '.bam')
					os.symlink(mt_bam, 'mt/' + instance_id + '.bam')
		CODE
	}
	output{
		Array[File] nuclear_bams = glob("nuclear/*.bam")
		Array[File] mt_bams = glob("mt/*.bam")
	}
	runtime{
		runtime_minutes: 10
		requested_memory_mb_per_core: 100
	}
}

# Remove duplicates that are already marked
# We do not deduplicate across libraries
task remove_marked_duplicates{
	Array[File] bams
	Array[String] references
	
	Int processes = 4
	command{
		python3 <<CODE
		from multiprocessing import Pool
		from os.path import basename
		import subprocess
		
		reference_string = "${sep=',' references}"
		references = reference_string.split(',')
		
		def remove_marked_duplicates_bam(input_bam):
			bam = basename(input_bam) # same name, in working directory
			# remove reference specific part of filename because downstream tools do not expect
			for reference in references:
				search_for = '.' + reference
				if search_for in bam:
					bam = bam.replace(search_for, '', 1)
			subprocess.run(['samtools', 'view', '-h', '-b', '-F', '0x400', '-o', bam, input_bam], check=True)
		
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(remove_marked_duplicates_bam, args=(bam, )) for bam in bams]
		pool.close()
		pool.join()
		for result in results:
			result.get()
		CODE
	}
	output{
		Array[File] no_duplicates_bams = glob("*.bam")
	}
	runtime{
		cpus: if length(bams) < processes then length(bams) else processes
		runtime_minutes: 100
		requested_memory_mb_per_core: 1000
	}
}

task coverage_without_index_barcode_key{
	Array[File] bam_stats
	Int reference_length
	String coverage_field
	
	command{
		python3 <<CODE
		from os.path import basename
		
		# adna SamStats assumes that the identifier is an index-barcode key 
		# trim trailing _ characters that are added when processing the filename as a key
		def cleaned_coverage(stats_file, coverage_field, reference_length):
			sample_id_from_filename = basename(stats_file)
			if sample_id_from_filename.find('.bam') >= 0:
				sample_id_from_filename = sample_id_from_filename[0:sample_id_from_filename.find('.bam')]
			elif sample_id_from_filename.find('.sam') >= 0:
				sample_id_from_filename = sample_id_from_filename[0:sample_id_from_filename.find('.sam')]
			
			with open(stats_file) as f:
				f.readline() # skip first line with read total
				for line in f:
					fields = line.strip().split('\t')
					sample_id = fields[0].rstrip('_') # remove key _ artifacts
					labels = fields[1:len(fields):2]
					values = fields[2:len(fields):2]
					if sample_id_from_filename == sample_id:
						for label, value in zip(labels, values):
							if coverage_field == label:
								coverage = float(value) / reference_length
								return (sample_id, coverage)
						return None
					raise ValueError('no ID match for %s' % (stats_file))
		
		bam_stats_string = "${sep=',' bam_stats}"
		stats_files = bam_stats_string.split(',')
		
		results = [cleaned_coverage(stats_file, "${coverage_field}", int(${reference_length}) ) for stats_file in stats_files]
		with open('coverage_statistics', 'w') as f:
			for result in results:
				if result is not None:
					print("%s\t%f" % result)
					print("%s\tMT_coverage\t%f" % result, file=f)
		
		CODE
	}
	output{
		File coverages = stdout()
		File coverage_statistics = "coverage_statistics"
	}
	runtime{
		runtime_minutes: 3
		requested_memory_mb_per_core: 100
	}
	
}

task analysis_results{
	File python_combine_dictionary_results
	Array[File] keyed_value_results
	
	command{
		python3 ${python_combine_dictionary_results} ${sep=' ' keyed_value_results} > results
	}
	output{
		File results = "results"
	}
	runtime{
		runtime_minutes: 5
		requested_memory_mb_per_core: 100
	}
}
