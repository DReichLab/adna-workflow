import "demultiplex.wdl" as demultiplex_master

workflow demultiplex_align_bams{
	# exclusive to Broad processing
	# end of exclusive section

	String blc_input_directory
	String dataset_label
	String date
	
	File i5_indices
	File i7_indices
	File barcodeSets
	
	File barcodes_q_only
	File labeled_i5
	File labeled_i7
	
	File adna_screen_jar
	File picard_jar
	
	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	
	Int samples_to_demultiplex
	
	File index_barcode_keys
	
	File python_lane_name
	File python_index_pairs_without_barcodes
	File python_common_unknown_barcodes
	File python_kmer_analysis
	File python_prepare_report
	
	# the references need to appear in the same directory as the derived files
	# in the prepare_reference, we put all of these into the same directory
	# all subsequent uses of the reference need to use that copy
	File reference_in
	File mt_reference_rsrs_in
	
	String output_path_parent
	String output_path = output_path_parent + "/" + date + "_" + dataset_label
	String output_path_nuclear_aligned_filtered = output_path + "/nuclear_aligned_filtered"
	String output_path_rsrs_aligned_filtered = output_path + "/rsrs_aligned_filtered"
	
	call demultiplex_master.prepare_reference as prepare_reference_nuclear{ input:
		reference = reference_in
	}
	call demultiplex_master.prepare_reference as prepare_reference_rsrs{ input:
		reference = mt_reference_rsrs_in
	}
	call demultiplex_master.versions{ input:
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		index_barcode_keys_to_stop_call_caching = index_barcode_keys
	}
	
	# exclusive to broad input
	call intake_fastq{
		# configured in inputs json
	}
	# 
	
	scatter(lane in intake_fastq.read_files_by_lane){
		call demultiplex_master.barcode_count_check{ input:
			adna_screen_jar = adna_screen_jar,
			i5_indices = i5_indices,
			i7_indices = i7_indices,
			barcodeSets = barcodeSets,
			read_files_by_lane = lane,
			reverse_complement_i5 = true
		}
	}
	call demultiplex_master.aggregate_statistics as aggregate_barcode_count_statistics{ input :
		adna_screen_jar=adna_screen_jar,
		statistics_by_group=barcode_count_check.barcode_count_statistics
	}
	scatter(lane in intake_fastq.read_files_by_lane){
		call discover_lane_name_from_filename_broad as discover_lane_name_from_filename{ input:
			filename = lane[0]
		}
		call demultiplex_master.merge_and_trim_lane { input : 
			adna_screen_jar = adna_screen_jar,
			i5_indices = i5_indices,
			i7_indices = i7_indices,
			barcodeSets = barcodeSets,
			read_files_by_lane = lane,
			label = discover_lane_name_from_filename.lane,
			barcode_count_statistics = aggregate_barcode_count_statistics.statistics,
			index_barcode_keys = index_barcode_keys,
			reverse_complement_i5 = true
		}
	}
	call demultiplex_master.collect_read_group_info{ input:
		read_groups_by_lane = merge_and_trim_lane.read_group
	}
	call demultiplex_master.aggregate_statistics as aggregate_lane_statistics{ input :
		adna_screen_jar=adna_screen_jar,
		statistics_by_group=merge_and_trim_lane.statistics
	}
	call demultiplex_master.prepare_demultiplex_report{ input:
		python_prepare_report = python_prepare_report,
		demultiplex_statistics = aggregate_lane_statistics.statistics,
		index_barcode_keys = index_barcode_keys,
		dataset_label = dataset_label,
		date = date
	}
	call demultiplex_master.collect_filenames{ input:
		filename_arrays = merge_and_trim_lane.fastq_to_align
	}
	String read_group = dataset_label
	scatter(fastq_to_align in collect_filenames.filenames){
		call demultiplex_master.align as align_nuclear{ input:
			missing_alignments_fraction = missing_alignments_fraction,
			max_open_gaps = max_open_gaps,
			seed_length = seed_length,
			fastq_to_align = fastq_to_align,
			reference = prepare_reference_nuclear.reference_fa,
			reference_amb = prepare_reference_nuclear.reference_amb,
			reference_ann = prepare_reference_nuclear.reference_ann,
			reference_bwt = prepare_reference_nuclear.reference_bwt,
			reference_pac = prepare_reference_nuclear.reference_pac,
			reference_sa = prepare_reference_nuclear.reference_sa
		}
	}
	call demultiplex_master.align_pool as align_rsrs{ input:
		missing_alignments_fraction = missing_alignments_fraction,
		max_open_gaps = max_open_gaps,
		seed_length = seed_length,
		fastq_to_align = collect_filenames.filenames,
		reference = prepare_reference_rsrs.reference_fa,
		reference_amb = prepare_reference_rsrs.reference_amb,
		reference_ann = prepare_reference_rsrs.reference_ann,
		reference_bwt = prepare_reference_rsrs.reference_bwt,
		reference_pac = prepare_reference_rsrs.reference_pac,
		reference_sa = prepare_reference_rsrs.reference_sa
	}
	call demultiplex_master.demultiplex as demultiplex_nuclear {input:
		adna_screen_jar = adna_screen_jar,
		prealignment_statistics = aggregate_lane_statistics.statistics,
		aligned_bam_files = align_nuclear.bam,
		samples_to_demultiplex = samples_to_demultiplex,
		index_barcode_keys = index_barcode_keys
	}
	call demultiplex_master.filter_aligned_only as filter_aligned_only_nuclear { input:
		picard_jar = picard_jar,
		bams = demultiplex_nuclear.demultiplexed_bam,
	}
	call demultiplex_master.demultiplex as demultiplex_rsrs {input:
		adna_screen_jar = adna_screen_jar,
		prealignment_statistics = aggregate_lane_statistics.statistics,
		aligned_bam_files = align_rsrs.bam,
		samples_to_demultiplex = samples_to_demultiplex,
		index_barcode_keys = index_barcode_keys
	}
	call demultiplex_master.filter_aligned_only as filter_aligned_only_rsrs{ input:
		picard_jar = picard_jar,
		bams = demultiplex_rsrs.demultiplexed_bam,
		minutes = 20
	}
	
	call demultiplex_master.index_pairs_without_barcodes{ input:
		python_index_pairs_without_barcodes = python_index_pairs_without_barcodes,
		barcode_count_statistics = aggregate_barcode_count_statistics.statistics
	}
	call demultiplex_master.demultiplex as demultiplex_for_unknown_barcodes{ input:
		adna_screen_jar = adna_screen_jar,
		prealignment_statistics = aggregate_lane_statistics.statistics,
		aligned_bam_files = align_rsrs.bam,
		samples_to_demultiplex = 0,
		index_barcode_keys = index_pairs_without_barcodes.index_pairs
	}
	call demultiplex_master.common_unknown_barcodes{ input:
		python_common_unknown_barcodes = python_common_unknown_barcodes,
		bams_without_known_barcodes = demultiplex_for_unknown_barcodes.demultiplexed_bam
	}
	call demultiplex_master.kmer_analysis{ input :
		python_kmer_analysis = python_kmer_analysis,
		python_prepare_report = python_prepare_report,
		barcodes_q_only = barcodes_q_only,
		labeled_i5 = labeled_i5,
		labeled_i7 = labeled_i7,
		counts_by_index_barcode_key = aggregate_lane_statistics.statistics,
		index_barcode_keys = index_barcode_keys,
		unknown_barcodes = common_unknown_barcodes.unknown_barcodes,
		dataset_label = dataset_label,
		date = date,
	}
	
	# output
	call demultiplex_master.copy_output as copy_nuclear_aligned_filtered{ input:
		files = filter_aligned_only_nuclear.filtered,
		output_path = output_path_nuclear_aligned_filtered
	}
	call demultiplex_master.copy_output as copy_rsrs_aligned_filtered{ input:
		files = filter_aligned_only_rsrs.filtered,
		output_path = output_path_rsrs_aligned_filtered
	}
	call demultiplex_master.copy_and_rename as copy_and_rename_demultiplex_nuclear_statistics{ input:
		source_file = demultiplex_nuclear.statistics,
		output_path = output_path,
		output_filename_no_path = "nuclear_statistics"
	}
	call demultiplex_master.copy_and_rename as copy_and_rename_demultiplex_mt_statistics{ input:
		source_file = demultiplex_rsrs.statistics,
		output_path = output_path,
		output_filename_no_path = "mt_statistics"
	}
	
	Array[File] misc_output_files = [collect_read_group_info.read_groups, kmer_analysis.analysis, versions.versions, prepare_demultiplex_report.report]
	call demultiplex_master.copy_output as copy_misc_output_files{input :
		files = misc_output_files,
		output_path = output_path
	}
	call demultiplex_master.copy_and_rename as copy_and_rename_lane_statistics{ input:
		source_file = aggregate_lane_statistics.statistics,
		output_path = output_path,
		output_filename_no_path = date + "_" + dataset_label + ".demultiplex_statistics"
	}
	
	call demultiplex_master.update_database_with_demultiplexed{ input:
		date_string = date,
		name = dataset_label,
		flowcell_by_lane = true,
		unused = (copy_nuclear_aligned_filtered.copied + copy_rsrs_aligned_filtered.copied + copy_and_rename_demultiplex_nuclear_statistics.copied + copy_and_rename_demultiplex_mt_statistics.copied + copy_misc_output_files.copied + copy_and_rename_lane_statistics.copied)
	}
	
	output{
		Array[File] nuclear_bams = filter_aligned_only_nuclear.filtered
		Array[File] rsrs_bams = filter_aligned_only_rsrs.filtered
		File kmer_analysis_report = kmer_analysis.analysis
		File aggregated_statistics = aggregate_lane_statistics.statistics
		File demultiplex_report = prepare_demultiplex_report.report
	}
}

task intake_fastq{
	String fastq_directory
	File python_arrange_lane_fastq
	
	command{
		python3 ${python_arrange_lane_fastq} ${fastq_directory + "/*.fastq.gz"} > files_by_lane
	}
	output{
		Array[Array[File]] read_files_by_lane = read_tsv('files_by_lane')
	}
	runtime{
		runtime_minutes: 2
		requested_memory_mb_per_core: 100
	}
}

task discover_lane_name_from_filename_broad{
	String filename
	
	command{
		python3 <<CODE
			from pathlib import Path
			with open('lane_name', 'w') as f:
				s = '${filename}'
				print(s[0], file=f)
		CODE
	}
	output{
		String lane = read_string('lane_name')
	}
	runtime{
		runtime_minutes: 2
		requested_memory_mb_per_core: 100
	}
}
