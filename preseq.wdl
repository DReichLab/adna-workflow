import "demultiplex.wdl" as demultiplex_align_bams
import "analysis.wdl" as analysis

workflow adna_analysis{
	File nuclear_bam_lists_to_merge
	File mt_bam_lists_to_merge
	File demultiplex_statistics_file_list
	Array[File] demultiplex_statistics_files = read_lines(demultiplex_statistics_file_list)
	File index_barcode_keys
	
	File labeled_barcodes
	File labeled_i5
	File labeled_i7
	
	String dataset_label
	String date
	
	String output_path_parent
	String output_path = output_path_parent + "/" + date + "_" + dataset_label

	File adna_screen_jar
	File picard_jar
	File pmdtools
	File htsbox
	File haplogrep_jar
	
	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	
	Int minimum_mapping_quality
	Int minimum_base_quality
	
	Int deamination_bases_to_clip
	
	File python_damage
	File python_damage_two_bases
	File python_target
	File python_central_measures
	File python_snp_target_bed
	File python_coverage
	File python_depth_histogram
	File python_angsd_results
	File python_prepare_report
	
	File coordinates_1240k_autosome
	File coordinates_1240k_x
	File coordinates_1240k_y
	
	# the references need to appear in the same directory as the derived files
	# in the prepare_reference, we put all of these into the same directory
	# all subsequent uses of the reference need to use that copy
	File mt_reference_rsrs_in
	File mt_reference_rcrs_in
	
	call demultiplex_align_bams.prepare_reference as prepare_reference_rsrs{ input:
		reference = mt_reference_rsrs_in
	}
	call demultiplex_align_bams.prepare_reference as prepare_reference_rcrs{ input:
		reference = mt_reference_rcrs_in
	}
	
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_across_sequencing_runs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = demultiplex_statistics_files
	}
	
	call analysis.combine_bams_into_libraries as combine_nuclear_libraries{ input:
		bam_lists = nuclear_bam_lists_to_merge,
		picard_jar = picard_jar
	}
	
	call analysis.preseq as preseq{ input:
		bams = combine_nuclear_libraries.library_bams,
		targets_bed = coordinates_1240k_autosome,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip,
		statistics = aggregate_statistics_across_sequencing_runs.statistics,	
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		python_depth_histogram = python_depth_histogram
	}	
}
