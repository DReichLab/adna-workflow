workflow ancientDNA_screen{
	String blc_input_directory
	String dataset_label
	String date
	
	File i5_indices
	File i7_indices
	File barcodeSets
	
	File barcodes_q_only
	
	File adna_screen_jar
	File picard_jar
	File pmdtools
	
	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	Int samples_to_demultiplex
	
	File python_lane_name
	File python_damage
	File python_damage_two_bases
	File python_target
	File python_central_measures
	File python_snp_target
	File python_snp_target_bed
	File python_coverage
	File python_floor
	File python_kmer_analysis
	
	File htsbox
	
	File spike3k_coordinates_autosome
	File spike3k_coordinates_x
	File spike3k_coordinates_y
	
	# the references need to appear in the same directory as the derived files
	# in the prepare_reference, we put all of these into the same directory
	# all subsequent uses of the reference need to use that copy
	File reference_in
	File mt_reference_rsrs_in
	File mt_reference_rcrs_in
	File mt_reference_human_95_consensus
	
	String output_path_parent
	String output_path = output_path_parent + "/" + date + "_" + dataset_label
	String output_path_hs37d5_aligned_unfiltered = output_path + "/hs37d5_aligned_unfiltered"
	String output_path_hs37d5_aligned_filtered = output_path + "/hs37d5_aligned_filtered"
	String output_path_rsrs_aligned_filtered = output_path + "/rsrs_aligned_filtered"

	call prepare_reference as prepare_reference_hs37d5{ input:
		reference = reference_in
	}
	call prepare_reference as prepare_reference_rsrs{ input:
		reference = mt_reference_rsrs_in
	}
	call prepare_reference as prepare_reference_rcrs{ input:
		reference = mt_reference_rcrs_in
	}
	call prepare_reference as prepare_reference_human_95_consensus{ input:
		reference = mt_reference_human_95_consensus
	}
	call bcl2fastq { input : blc_input_directory=blc_input_directory} 
	scatter(lane in bcl2fastq.read_files_by_lane){
		call barcode_count_check{ input:
			adna_screen_jar = adna_screen_jar,
			i5_indices = i5_indices,
			i7_indices = i7_indices,
			barcodeSets = barcodeSets,
			read_files_by_lane = lane
		}
	}
	call aggregate_statistics as aggregate_barcode_count_statistics{ input :
		adna_screen_jar=adna_screen_jar,
		statistics_by_group=barcode_count_check.barcode_count_statistics
	}
	scatter(lane in bcl2fastq.read_files_by_lane){
		call discover_lane_name_from_filename{ input:
			python_lane_name = python_lane_name,
			filename = lane[0]
		}
		call merge_and_trim_lane { input : 
			adna_screen_jar = adna_screen_jar,
			i5_indices = i5_indices,
			i7_indices = i7_indices,
			barcodeSets = barcodeSets,
			read_files_by_lane = lane,
			label = discover_lane_name_from_filename.lane,
			barcode_count_statistics = aggregate_barcode_count_statistics.statistics
		}
	}
	call aggregate_statistics as aggregate_lane_statistics{ input :
		adna_screen_jar=adna_screen_jar,
		statistics_by_group=merge_and_trim_lane.statistics
	}
	call kmer_analysis{ input :
		python_kmer_analysis = python_kmer_analysis,
		barcodes_q_only = barcodes_q_only,
		counts_by_index_barcode_key = aggregate_lane_statistics.statistics,
		dataset_label = dataset_label,
		date = date
	}
	call collect_filenames{ input:
		filename_arrays = merge_and_trim_lane.fastq_to_align
	}
	String read_group = dataset_label
	scatter(fastq_to_align in collect_filenames.filenames){
		call align as align_hs37d5{ input:
			missing_alignments_fraction = missing_alignments_fraction,
			max_open_gaps = max_open_gaps,
			seed_length = seed_length,
			fastq_to_align = fastq_to_align,
			reference = prepare_reference_hs37d5.reference_fa,
			reference_amb = prepare_reference_hs37d5.reference_amb,
			reference_ann = prepare_reference_hs37d5.reference_ann,
			reference_bwt = prepare_reference_hs37d5.reference_bwt,
			reference_pac = prepare_reference_hs37d5.reference_pac,
			reference_sa = prepare_reference_hs37d5.reference_sa
		}
		call align as align_rsrs{ input:
			missing_alignments_fraction = missing_alignments_fraction,
			max_open_gaps = max_open_gaps,
			seed_length = seed_length,
			fastq_to_align = fastq_to_align,
			reference = prepare_reference_rsrs.reference_fa,
			reference_amb = prepare_reference_rsrs.reference_amb,
			reference_ann = prepare_reference_rsrs.reference_ann,
			reference_bwt = prepare_reference_rsrs.reference_bwt,
			reference_pac = prepare_reference_rsrs.reference_pac,
			reference_sa = prepare_reference_rsrs.reference_sa
		}
	}
	call demultiplex as demultiplex_hs37d5 {input:
		adna_screen_jar = adna_screen_jar,
		prealignment_statistics = aggregate_lane_statistics.statistics,
		aligned_bam_files = align_hs37d5.bam,
		samples_to_demultiplex = samples_to_demultiplex
	}
	call chromosome_target as hs37d5_chromosome_target{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = demultiplex_hs37d5.demultiplexed_bam,
		targets="\"{'autosome_pre':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X_pre':'X','Y_pre':'Y','human_pre':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	call filter_aligned_only as filter_aligned_only_hs37d5 { input:
		bams = demultiplex_hs37d5.demultiplexed_bam
	}
	call snp_target_bed as spike3k_pre{ input:
		coordinates_autosome = spike3k_coordinates_autosome,
		coordinates_x = spike3k_coordinates_x,
		coordinates_y = spike3k_coordinates_y,
		bams = filter_aligned_only_hs37d5.filtered,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip,
		label = "spike3k_pre",
		picard_jar = picard_jar,
		adna_screen_jar = adna_screen_jar,
		python_snp_target_bed = python_snp_target_bed
	}
	call duplicates as duplicates_hs37d5 { input: 
		picard_jar = picard_jar,
		adna_screen_jar = adna_screen_jar,
		pmdtools = pmdtools,
		unsorted = filter_aligned_only_hs37d5.filtered,
		duplicates_label = "duplicates_hs37d5"
	}
	call snp_target_bed as spike3k_post{ input:
		coordinates_autosome = spike3k_coordinates_autosome,
		coordinates_x = spike3k_coordinates_x,
		coordinates_y = spike3k_coordinates_y,
		bams = duplicates_hs37d5.aligned_deduplicated,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip,
		label = "spike3k_post",
		picard_jar = picard_jar,
		adna_screen_jar = adna_screen_jar,
		python_snp_target_bed = python_snp_target_bed
	}
	call damage_loop as damage_loop_hs37d5{ input:
		pmdtools = pmdtools,
		python_damage_two_bases = python_damage_two_bases,
		bams = duplicates_hs37d5.aligned_deduplicated,
		damage_label = "damage_hs37d5",
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality
	}
	call chromosome_target as hs37d5_chromosome_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = duplicates_hs37d5.aligned_deduplicated,
		targets="\"{'autosome_post':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X_post':'X','Y_post':'Y','human_post':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	
	call demultiplex as demultiplex_rsrs {input:
		adna_screen_jar = adna_screen_jar,
		prealignment_statistics = aggregate_lane_statistics.statistics,
		aligned_bam_files = align_rsrs.bam,
		samples_to_demultiplex = samples_to_demultiplex
	}
	call chromosome_target as rsrs_chromosome_target{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = demultiplex_rsrs.demultiplexed_bam,
		targets="\"{'MT_pre':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	call filter_aligned_only as filter_aligned_only_rsrs{ input:
		bams = demultiplex_rsrs.demultiplexed_bam
	}
	call duplicates as duplicates_rsrs{ input: 
		picard_jar = picard_jar,
		adna_screen_jar = adna_screen_jar,
		pmdtools = pmdtools,
		unsorted = filter_aligned_only_rsrs.filtered,
		duplicates_label = "duplicates_rsrs"
	}
	call haplogrep as haplogrep_rcrs{ input:
		missing_alignments_fraction = missing_alignments_fraction,
		max_open_gaps = max_open_gaps,
		seed_length = seed_length,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip,
		bams = duplicates_rsrs.aligned_deduplicated,
		reference = prepare_reference_rcrs.reference_fa,
		reference_amb = prepare_reference_rcrs.reference_amb,
		reference_ann = prepare_reference_rcrs.reference_ann,
		reference_bwt = prepare_reference_rcrs.reference_bwt,
		reference_pac = prepare_reference_rcrs.reference_pac,
		reference_sa = prepare_reference_rcrs.reference_sa,
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar
	}
	call damage_loop as damage_loop_rsrs{ input:
		pmdtools = pmdtools,
		python_damage_two_bases = python_damage_two_bases,
		bams = duplicates_rsrs.aligned_deduplicated,
		damage_label = "damage_rsrs",
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality
	}
	call chromosome_target as rsrs_chromosome_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = duplicates_rsrs.aligned_deduplicated,
		targets="\"{'MT_post':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	call chromosome_coverage as rsrs_coverage{ input:
		bam_stats = rsrs_chromosome_target_post.target_stats,
		python_coverage = python_coverage,
		reference_length = 16569,
		coverage_field = "MT_post-coverageLength"
	}
	scatter(bam in duplicates_rsrs.aligned_deduplicated){
#		call schmutzi{ input:
#			bam = duplicates_rsrs.aligned_deduplicated,
#			picard_jar = picard_jar,
#			reference = prepare_reference_rsrs.reference_fa,
#			reference_amb = prepare_reference_rsrs.reference_amb,
#			reference_ann = prepare_reference_rsrs.reference_ann,
#			reference_bwt = prepare_reference_rsrs.reference_bwt,
#			reference_pac = prepare_reference_rsrs.reference_pac,
#			reference_sa = prepare_reference_rsrs.reference_sa,
#			reference_fai = prepare_reference_rsrs.reference_fai,
# TODO when using schmutzi, need to read coverage in from file as in contammix
#			coverage = chromosome_target_single_rsrs.coverage
#		}
#		call contamination_rare_variant{ input:
#			bam = duplicates_rsrs.aligned_deduplicated,
#			picard_jar = picard_jar,
#			adna_screen_jar = adna_screen_jar,
#			missing_alignments_fraction = missing_alignments_fraction,
#			max_open_gaps = max_open_gaps,
#			seed_length = seed_length,
#			minimum_mapping_quality = minimum_mapping_quality,
#			minimum_base_quality = minimum_base_quality,
#			deamination_bases_to_clip = deamination_bases_to_clip,
#			reference = prepare_reference_human_95_consensus.reference_fa,
#			reference_amb = prepare_reference_human_95_consensus.reference_amb,
#			reference_ann = prepare_reference_human_95_consensus.reference_ann,
#			reference_bwt = prepare_reference_human_95_consensus.reference_bwt,
#			reference_pac = prepare_reference_human_95_consensus.reference_pac,
#			reference_sa = prepare_reference_human_95_consensus.reference_sa,
#			reference_fai = prepare_reference_human_95_consensus.reference_fai
#		}
		call contammix{ input:
			bam = duplicates_rsrs.aligned_deduplicated,
			picard_jar = picard_jar,
			htsbox = htsbox,
			missing_alignments_fraction = missing_alignments_fraction,
			max_open_gaps = max_open_gaps,
			seed_length = seed_length,
			minimum_mapping_quality = minimum_mapping_quality,
			minimum_base_quality = minimum_base_quality,
			deamination_bases_to_clip = deamination_bases_to_clip,
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
	
	call central_measures as central_measures_hs37d5{ input:
		python_central_measures = python_central_measures,
		mean_label = "mean_hs37d5",
		median_label = "median_hs37d5",
		histograms = hs37d5_chromosome_target_post.length_histogram
	}
	call central_measures as central_measures_rsrs{ input:
		python_central_measures = python_central_measures,
		mean_label = "mean_rsrs",
		median_label = "median_rsrs",
		histograms = rsrs_chromosome_target_post.length_histogram
	}
	call summarize_haplogroups{ input:
		haplogrep_output = haplogrep_rcrs.haplogroup_report
	}

	call aggregate_statistics as aggregate_statistics_duplicates_hs37d5{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = duplicates_hs37d5.duplicates_statistics
	}
	call aggregate_statistics as aggregate_statistics_duplicates_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = duplicates_rsrs.duplicates_statistics
	}
	
	call aggregate_statistics as aggregate_statistics_pre_hs37d5{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = hs37d5_chromosome_target.target_stats
	}
	call aggregate_statistics as aggregate_statistics_post_hs37d5{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = hs37d5_chromosome_target_post.target_stats
	}
	call aggregate_statistics as aggregate_statistics_pre_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = rsrs_chromosome_target.target_stats
	}
	call aggregate_statistics as aggregate_statistics_post_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = rsrs_chromosome_target_post.target_stats
	}
	
	Array[File] cumulative_statistics = [
		aggregate_lane_statistics.statistics,
		aggregate_statistics_duplicates_hs37d5.statistics,
		aggregate_statistics_duplicates_rsrs.statistics,
		aggregate_statistics_pre_hs37d5.statistics,
		aggregate_statistics_post_hs37d5.statistics,
		aggregate_statistics_pre_rsrs.statistics,
		aggregate_statistics_post_rsrs.statistics
	]
	call aggregate_statistics as aggregate_statistics_final{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = cumulative_statistics
	}
	
	call copy_output as copy_hs37d5_aligned_unfiltered{ input:
		files = demultiplex_hs37d5.demultiplexed_bam,
		output_path = output_path_hs37d5_aligned_unfiltered
	}
	call copy_output as copy_hs37d5_aligned_filtered{ input:
		files = duplicates_hs37d5.aligned_deduplicated,
		output_path = output_path_hs37d5_aligned_filtered
	}
	call copy_output as copy_hs37d5_histogram{ input:
		files = hs37d5_chromosome_target_post.length_histogram,
		output_path = output_path_hs37d5_aligned_filtered
	}
	call copy_output as copy_rsrs_aligned_filtered{ input:
		files = duplicates_rsrs.aligned_deduplicated,
		output_path = output_path_rsrs_aligned_filtered
	}
	call copy_output as copy_rsrs_histogram{ input:
		files = rsrs_chromosome_target_post.length_histogram,
		output_path = output_path_rsrs_aligned_filtered
	}
	call copy_output as copy_consensus_mt{ input:
		files = contammix.consensus,
		output_path = output_path_rsrs_aligned_filtered
	}
	Array[File] kmer_array = [kmer_analysis.analysis]
	call copy_output as copy_kmer{input :
		files = kmer_array,
		output_path = output_path
	}
	call concatenate as concatenate_spike3k_pre{ input:
		to_concatenate = spike3k_pre.snp_target_stats
	}
	call concatenate as concatenate_spike3k_post{ input:
		to_concatenate = spike3k_post.snp_target_stats
	}
#	call concatenate as concatenate_schmutzi{ input:
#		to_concatenate = schmutzi.contamination_estimate
#	}
#	call concatenate as concatenate_contamination_rare_variant{ input:
#		to_concatenate = contamination_rare_variant.contamination_estimate
#	}
	call concatenate as concatenate_contammix{ input:
		to_concatenate = contammix.contamination_estimate
	}
	call spike3k_complexity{ input:
		spike3k_pre_data = concatenate_spike3k_pre.concatenated,
		spike3k_post_data =concatenate_spike3k_post.concatenated
	}
	
	Array[File] final_keyed_statistics = [
		damage_loop_hs37d5.damage_all_samples_two_bases,
		damage_loop_rsrs.damage_all_samples_two_bases,
		central_measures_hs37d5.central_measures_output,
		central_measures_rsrs.central_measures_output,
		summarize_haplogroups.haplogroups,
		concatenate_spike3k_pre.concatenated,
		concatenate_spike3k_post.concatenated,
		spike3k_complexity.estimates,
#		concatenate_schmutzi.concatenated,
#		concatenate_contamination_rare_variant.concatenated,
		concatenate_contammix.concatenated
	]
	call prepare_report{ input:
		aggregated_statistics = aggregate_statistics_final.statistics,
		keyed_statistics = final_keyed_statistics,
		dataset_label = dataset_label,
		date = date
	}
	Array[File] report_array = [prepare_report.report]
	call copy_output as copy_report{ input:
		files = report_array,
		output_path = output_path
	}
}

# output needs to have all files in the same directory
task prepare_reference{
	File reference
	String filename = basename(reference)
	
	# java -jar ${picard_jar} CreateSequenceDictionary R=${reference} O=${reference}.dict
	command{
		set -e
		bwa index ${reference}
		samtools faidx ${reference}
		
		cp -s ${reference}     ${filename}
		cp -l ${reference}.amb ${filename}.amb
		cp -l ${reference}.ann ${filename}.ann
		cp -l ${reference}.bwt ${filename}.bwt
		cp -l ${reference}.pac ${filename}.pac
		cp -l ${reference}.sa  ${filename}.sa
		cp -l ${reference}.fai ${filename}.fai
	}
	output{
		File reference_fa  = "${filename}"
		File reference_amb = "${filename}.amb"
		File reference_ann = "${filename}.ann"
		File reference_bwt = "${filename}.bwt"
		File reference_pac = "${filename}.pac"
		File reference_sa  = "${filename}.sa"
		File reference_fai = "${filename}.fai"
	}
	runtime{
		cpus: 4
		requested_memory_mb_per_core: 8192
	}
}

task bcl2fastq{
	String blc_input_directory

	command{
		set -e
		touch empty
		bcl2fastq \
			-R ${blc_input_directory} \
			-o ./ \
			--create-fastq-for-index-reads \
			--sample-sheet empty
	}
	
	output{
		Array[Array[File]] read_files_by_lane = [
			[
				"Undetermined_S0_L001_R1_001.fastq.gz", 
				"Undetermined_S0_L001_R2_001.fastq.gz",
				"Undetermined_S0_L001_I1_001.fastq.gz",
				"Undetermined_S0_L001_I2_001.fastq.gz"
			],
			[
				"Undetermined_S0_L002_R1_001.fastq.gz",
				"Undetermined_S0_L002_R2_001.fastq.gz",
				"Undetermined_S0_L002_I1_001.fastq.gz",
				"Undetermined_S0_L002_I2_001.fastq.gz"
			],
			[
				"Undetermined_S0_L003_R1_001.fastq.gz",
				"Undetermined_S0_L003_R2_001.fastq.gz",
				"Undetermined_S0_L003_I1_001.fastq.gz",
				"Undetermined_S0_L003_I2_001.fastq.gz"
			],
			[
				"Undetermined_S0_L004_R1_001.fastq.gz",
				"Undetermined_S0_L004_R2_001.fastq.gz",
				"Undetermined_S0_L004_I1_001.fastq.gz",
				"Undetermined_S0_L004_I2_001.fastq.gz"
			]
		]
	}
	runtime{
		cpus: 4
		requested_memory_mb_per_core: 8192
	}
}

task discover_lane_name_from_filename{
	String filename
	File python_lane_name
	
	command{
		python ${python_lane_name} ${filename} > lane_name
	}
	output{
		String lane = read_string("./lane_name")
	}
	runtime{
		runtime_minutes: 10
		requested_memory_mb_per_core: 1024
	}
}

# Count the number of paired reads that would demultiplex with barcodes, and those without
task barcode_count_check{
	File adna_screen_jar
	File i5_indices
	File i7_indices
	File barcodeSets
	Array[File] read_files_by_lane
	
	command{
		java -Xmx14g -jar ${adna_screen_jar} BarcodeCount --i5-indices ${i5_indices} --i7-indices ${i7_indices} --barcodes ${barcodeSets} ${read_files_by_lane[0]} ${read_files_by_lane[1]} ${read_files_by_lane[2]} ${read_files_by_lane[3]} > barcodeCount.stats
	}
	output{
		File barcode_count_statistics = "barcodeCount.stats"
	}
	runtime{
		requested_memory_mb_per_core: 16384
	}
}

task merge_and_trim_lane{
	File adna_screen_jar
	File i5_indices
	File i7_indices
	File barcodeSets
	Array[File] read_files_by_lane
	String label
	File barcode_count_statistics
	
	command{
		java -Xmx14g -jar ${adna_screen_jar} IndexAndBarcodeScreener --i5-indices ${i5_indices} --i7-indices ${i7_indices} --barcodes ${barcodeSets} --barcode-count ${barcode_count_statistics} ${read_files_by_lane[0]} ${read_files_by_lane[1]} ${read_files_by_lane[2]} ${read_files_by_lane[3]} ${label} > ${label}.stats
	}
	
	output{
		Array[File] fastq_to_align = glob("${label}*.fastq.gz")
		File statistics = "${label}.stats"
	}
	runtime{
		requested_memory_mb_per_core: 16384
	}
}

task aggregate_statistics{
	File adna_screen_jar
	Array [File] statistics_by_group
	
	command{
		java -jar ${adna_screen_jar} AggregateStatistics ${sep=' ' statistics_by_group} > aggregated_statistics
	}
	output{
		File statistics = "aggregated_statistics"
	}
	runtime{
		runtime_minutes: 60
		requested_memory_mb_per_core: 4096
	}
}

task align{
	File fastq_to_align
	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	Int threads
	
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_pac
	File reference_sa
	
	# the bwa -r option for specifying the read group leads to problems
	# bwa uses tab delimiters, but these are illegal in the program group section of a sam file
	# so we leave out the read group here and plan to insert these at the end of processing
	command {
		set -e
		bwa aln -t ${threads} -o ${max_open_gaps} -n ${missing_alignments_fraction} -l ${seed_length} ${reference} ${fastq_to_align} > aligned.sai
		bwa samse ${reference} aligned.sai ${fastq_to_align} | samtools view -bS - > aligned.bam
	}
	output{
		File bam = "aligned.bam"
	}
	runtime{
		cpus: "${threads}"
		requested_memory_mb_per_core: 8192
	}
}

# use String instead of filename to avoid file copying overhead
task collect_filenames{
	Array[Array[String]] filename_arrays
	File python_flatten
	
	command{
		echo "${sep='\n' filename_arrays}" > raw_array
		python ${python_flatten} < raw_array > file_of_filenames
	}
	output{
		Array[String] filenames = read_lines("./file_of_filenames")
	}
	runtime{
		runtime_minutes: 60
		requested_memory_mb_per_core: 2048
	}
}

task demultiplex{
	File adna_screen_jar
	File prealignment_statistics
	Array[File] aligned_bam_files
	Int samples_to_demultiplex
	File? index_barcode_keys
	File? barcodes
	
	command{
		java -Xmx14g -jar ${adna_screen_jar} DemultiplexSAM -b -n ${samples_to_demultiplex} -s ${prealignment_statistics} ${"-e " + index_barcode_keys} ${"--barcodeFile " + barcodes} ${sep=' ' aligned_bam_files} > postalignment_statistics
	}
	output{
		Array[File] demultiplexed_bam = glob("*.bam")
		File statistics = "postalignment_statistics"
	}
	runtime{
		cpus: 2
		requested_memory_mb_per_core: 8000
	}
}

# filter out unaligned reads
task filter_aligned_only{
	Array[File] bams
	Int processes = 10
	
	# picard complains "MAPQ should be 0 for unmapped read." while trying to filter unmapped reads
	#java -jar ${picard_jar} SortSam I=${bam} O=sorted_queryname.bam SORT_ORDER=queryname
	#java -jar ${picard_jar} FilterSamReads I=sorted_queryname.bam O=${filename} FILTER=includeAligned
	command{
		set -e
		python <<CODE
		from multiprocessing import Pool
		from os.path import basename
		import subprocess
		
		def filter_bam_aligned_only(bam):
			output_filename = basename(bam)
			subprocess.check_output("samtools view -h -b -F 4 -o %s %s" % (output_filename, bam), shell=True)
		
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		[pool.apply_async(filter_bam_aligned_only, args=(bam,)) for bam in bams]
		pool.close()
		pool.join()
		CODE
	}
	output{
		Array[File] filtered = glob("*.bam")
	}
	runtime{
		cpus: processes
		requested_memory_mb_per_core: 1000
	}
}

task duplicates{
	File picard_jar
	File adna_screen_jar
	File pmdtools
	Array[File] unsorted
	String duplicates_label
	
	Int processes = 8
	
	command{
		set -e
		
		python <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def deduplicate_bam(bam):
			sample_id_filename = basename(bam)
			sample_id_filename_no_extension, extension = splitext(sample_id_filename)
			
			sorted_bam = sample_id_filename_no_extension + ".sorted_coordinate.bam"
			
			subprocess.check_output("java -Xmx9g -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (bam, sorted_bam), shell=True)
			subprocess.check_output("java -Xmx9g -jar ${picard_jar} MarkDuplicates I=%s O=%s M=%s.dedup_stats REMOVE_DUPLICATES=true BARCODE_TAG=XD ADD_PG_TAG_TO_READS=false MAX_FILE_HANDLES=1000" % (sorted_bam, sample_id_filename, sample_id_filename), shell=True)
			subprocess.check_output("java -Xmx9g -jar ${adna_screen_jar} ReadMarkDuplicatesStatistics -l ${duplicates_label} %s.dedup_stats > %s.stats" % (sample_id_filename, sample_id_filename), shell=True)
		
		bams_string = "${sep=',' unsorted}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		[pool.apply_async(deduplicate_bam, args=(bam,)) for bam in bams]
		pool.close()
		pool.join()
		CODE
	}
	output{
		Array[File] aligned_deduplicated = glob("*.bam")
		Array[File] duplicates_statistics = glob("*.stats")
	}
	runtime{
		cpus: processes
		requested_memory_mb_per_core: 10000
	}
}

task damage_loop{
	File pmdtools
	File python_damage_two_bases
	Array[File] bams
	String damage_label
	Int minimum_mapping_quality
	Int minimum_base_quality
	
	Int processes = 10
	
	command{
		set -e
		python <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def damage_for_bam(bam):
			sample_id_filename = basename(bam)
			sample_id_filename_no_extension, extension = splitext(sample_id_filename)
			
			damage_filename = sample_id_filename_no_extension + ".damage"
			
			subprocess.check_output("samtools view %s | python ${pmdtools} -d --requiremapq=${minimum_mapping_quality} --requirebaseq=${minimum_base_quality} > %s" %(bam, damage_filename), shell=True)
			damage_result = subprocess.check_output("python ${python_damage_two_bases} ${damage_label} %s" % (damage_filename,), shell=True)
			return damage_result.strip()
			
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(damage_for_bam, args=(bam, resultsQueue)) for bam in bams]
		pool.close()
		pool.join()
		with open('damage_all_samples_two_bases', 'w') as f:
			for result in results:
				f.write(result.get())
				f.write('\n')
		CODE
	}
	output{
		File damage_all_samples_two_bases = "damage_all_samples_two_bases"
	}
	runtime{
		cpus: processes
		requested_memory_mb_per_core: 2000
	}
}

# Alternative in place of looping in WDL, run loop in python
task chromosome_target{
	File python_target
	File adna_screen_jar
	Array[File] bams
	String targets
	Int minimum_mapping_quality
	
	#String sample_id_filename = basename(bam, ".bam")

	command{
		python ${python_target} ${adna_screen_jar} ${targets} ${minimum_mapping_quality} ${sep=' ' bams}
	}
	output{
		Array[File] target_stats = glob("*.stats")
		Array[File] length_histogram = glob("*.histogram")
	}
	runtime{
		cpus: 2
		requested_memory_mb_per_core: 8192
	}
}

task chromosome_coverage{
	Array[File] bam_stats
	File python_coverage
	Int reference_length
	String coverage_field
	
	command{
		set -e
		for bam in ${sep=' ' bam_stats}
		do
			sample_id = basename($bam, ".stats")
			python ${python_coverage} $sample_id.stats ${reference_length} $sample_id ${coverage_field} >> coverages
		done
	}
	output{
		File coverages = "coverages"
	}
	runtime{
		requested_memory_mb_per_core: 2000
	}
}

# This counts reads similarly to snp_target, but uses a bed file to look for read overlaps in an interval around the SNP. 
task snp_target_bed{
	File coordinates_autosome
	File coordinates_x
	File coordinates_y
	Array[File] bams
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	String label
	File picard_jar
	File adna_screen_jar
	
	File python_snp_target_bed
	
	Int processes = 10

	command{
		set -e
		
		python <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def count_SNPs(bam):
			sample_id_filename = basename(bam)
			sample_id_filename_no_extension, extension = splitext(sample_id_filename)
			
			clipped_bam = sample_id_filename_no_extension + ".clipped.bam"
			clipped_sorted_bam = sample_id_filename_no_extension + ".clipped.sorted.bam"
			
			subprocess.check_output("java -Xmx2500m -jar ${adna_screen_jar} softclip -b -n ${deamination_bases_to_clip} -i %s -o %s" % (bam, clipped_bam), shell=True)
			subprocess.check_output("java -Xmx2500m -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (clipped_bam, clipped_sorted_bam), shell=True)
			subprocess.check_output("samtools index %s" % (clipped_sorted_bam,), shell=True)
			subprocess.check_output("samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_autosome} %s > %s.autosome" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_x}        %s > %s.x" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_y}        %s > %s.y" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("python ${python_snp_target_bed} ${label} %s.autosome %s.x %s.y > %s.snp_target_stats" % (sample_id_filename, sample_id_filename, sample_id_filename, sample_id_filename), shell=True)
			
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		[pool.apply_async(count_SNPs, args=(bam,)) for bam in bams]
		pool.close()
		pool.join()
		CODE
	}
	output{
		Array[File] snp_target_stats = glob("*.snp_target_stats")
	}
	runtime{
		cpus: processes
		requested_memory_mb_per_core: 3000
	}
}

task spike3k_complexity{
	File python_spike3k_complexity_prep
	File python_spike3k_complexity_results
	File spike3k_complexity_binary
	File spike3k_pre_data
	File spike3k_post_data
	
	command{
		set -e
		python ${python_spike3k_complexity_prep} ${spike3k_pre_data} ${spike3k_post_data} > spike3k_for_complexity
		${spike3k_complexity_binary} -i spike3k_for_complexity -o nick_table
		python ${python_spike3k_complexity_results} nick_table > spike3k_complexity_estimates
	}
	output{
		File estimates = "spike3k_complexity_estimates"
	}
	runtime{
		runtime_minutes: 240
		requested_memory_mb_per_core: 4096
	}
}

task concatenate{
	Array[File] to_concatenate
	
	command{
		for file in ${sep=' ' to_concatenate}  ; do 
			cat $file >> concatenated
		done
	}
	output{
		File concatenated = "concatenated"
	}
	runtime{
		runtime_minutes: 60
		requested_memory_mb_per_core: 4096
	}
}

task copy_output{
	Array[File] files
	String output_path
	
	command{
		mkdir -p ${output_path};
		for file in ${sep=' ' files}  ; do 
			cp -l $file "${output_path}" || cp $file "${output_path}"
		done
	}
	runtime{
		requested_memory_mb_per_core: 2048
	}
}

task haplogrep{
	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	String region = "MT"
	Array[File] bams
	File haplogrep_jar
	Int phylotree_version
	File adna_screen_jar
	File picard_jar
	
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_pac
	File reference_sa
	
	Int processes = 10
	
	command{
		set -e
		
		python <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
				
		def haplogrep_run(bam):
			sample_id_filename = basename(bam)
			sample_id, extension = splitext(sample_id_filename)
			
			subprocess.check_output("java -jar ${picard_jar} SamToFastq I=%s FASTQ=%s.fastq" % (bam, sample_id), shell=True)
			subprocess.check_output("bwa aln -t 2 -o ${max_open_gaps} -n ${missing_alignments_fraction} -l ${seed_length} ${reference} %s.fastq > %s.sai" % (sample_id, sample_id), shell=True)
			subprocess.check_output("bwa samse ${reference} %s.sai %s.fastq | samtools view -bS - > %s.realigned.bam" % (sample_id, sample_id, sample_id), shell=True)
			subprocess.check_output("java -jar ${adna_screen_jar} softclip -b -n ${deamination_bases_to_clip} -i %s.realigned.bam -o %s.clipped_unsorted_realigned.bam" % (sample_id, sample_id), shell=True)
			subprocess.check_output("java -jar ${picard_jar} SortSam I=%s.clipped_unsorted_realigned.bam O=%s.bam SORT_ORDER=coordinate" % (sample_id, sample_id), shell=True)
			subprocess.check_output("samtools index %s.bam" % (sample_id,), shell=True)
			subprocess.check_output("samtools mpileup -q ${minimum_mapping_quality} -Q ${minimum_base_quality} -r ${region} -u -f ${reference} %s.bam | bcftools call -c -v --ploidy 1 > %s.vcf" % (sample_id, sample_id), shell=True)
			subprocess.check_output("java -jar ${haplogrep_jar} --format vcf --phylotree ${phylotree_version} --in %s.vcf --out %s.haplogroup" % (sample_id, sample_id), shell=True)
		
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		[pool.apply_async(haplogrep_run, args=(bam,)) for bam in bams]
		pool.close()
		pool.join()
		CODE
	}
	output{
		Array[File] haplogroup_report = glob("*.haplogroup")
	}
	runtime{
		cpus: processes
		requested_memory_mb_per_core: 2000
	}
}

task summarize_haplogroups{
	File python_haplogroup
	Array[File] haplogrep_output
	
	command{
		python ${python_haplogroup} ${sep=' ' haplogrep_output} > haplogroups
	}
	output{
		File haplogroups = "haplogroups"
	}
	runtime{
		requested_memory_mb_per_core: 4096
	}
}

task central_measures{
	File python_central_measures
	String median_label
	String mean_label
	Array[File] histograms
	
	command{
		python ${python_central_measures} ${median_label} ${mean_label} ${sep=' ' histograms} > central_measures
	}
	output{
		File central_measures_output = "central_measures"
	}
	runtime{
		requested_memory_mb_per_core: 4096
	}
}

task schmutzi{
	File bam
	String sample_id_filename = basename(bam)
	String key = sub(sample_id_filename, ".bam$", "") # remove file extension
	Int deamination_length
	Float coverage
	
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_pac
	File reference_sa
	File reference_fai
	
	File python_schumtzi_output
	File picard_jar
	
	# The schmutzi scripts require many other programs and files, which we need to include
	File schmutzi_contDeam
	File schmutzi_pl
	String path_to_eurasion_freqs
	File schmutzi_approxDist
	File schmutzi_bam2prof
	File schmutzi_contDeam_pl
	File schmutzi_countRecords
	File schmutzi_endoCaller
	File schmutzi_filterlog
	File schmutzi_insertSize
	File schmutzi_jointFreqDeaminated
	File schmutzi_jointFreqDeaminatedDouble
	File schmutzi_log2ConsensusLog
	File schmutzi_log2fasta
	File schmutzi_log2freq
	File schmutzi_logdiff
	File schmutzi_logmask
	File schmutzi_logs2pos
	File schmutzi_msa2freq
	File schmutzi_msa2log
	File schmutzi_msa2singlefreq
	File schmutzi_mtCont
	File schmutzi_parseSchmutzi
	File schmutzi_posteriorDeam
	File schmutzi_wrapper_pl
	File schmutzi_wrapperRMDUP_pl
	
	File schmutzi_illuminaProf_error
	File schmutzi_illuminaProf_null
	
	File schmutzi_splitEndoVsCont_denisovaHuman
	File schmutzi_splitEndoVsCont_neandertalHuman
	File schmutzi_splitEndoVsCont_poshap2splitbam
	
	# We downsample to 300x MT coverage
	Float threshold = 300.0
	Float retain_probability = if (coverage > threshold) then (threshold / coverage) else 1.0
	Int threads = if (coverage >= 50) then 8 else 4
	Int iterations

	# some of these commands may fail
	# the python command will report nan in this case
	command{
		java -jar ${picard_jar} DownsampleSam I=${bam} O=downsampled.bam PROBABILITY=${retain_probability}
		samtools calmd -b downsampled.bam ${reference} > schmutzi.bam
		samtools index schmutzi.bam
		${schmutzi_contDeam_pl} --lengthDeam ${deamination_length} --library single --out ${key} ${reference} schmutzi.bam
		${schmutzi_pl} --iterations ${iterations} -t ${threads} --notusepredC --uselength --ref ${reference} --out ${key}_npred ${key} ${path_to_eurasion_freqs} schmutzi.bam
		${schmutzi_pl} --iterations ${iterations} -t ${threads}               --uselength --ref ${reference} --out ${key}_wpred ${key} ${path_to_eurasion_freqs} schmutzi.bam
		python ${python_schumtzi_output} ${key} ${key}_wpred_final.cont.est > contamination_estimate
	}
	output{
		File contamination_estimate = "contamination_estimate"
	}
	runtime{
		cpus: threads
		requested_memory_mb_per_core: if (threads <= 4) then 8000 else 4000
	}
}

task contamination_rare_variant{
	File bam
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	
	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	Int threads = 1

	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_pac
	File reference_sa
	File reference_fai
	
	File picard_jar
	File adna_screen_jar
	File python_calico
	File python_contamination_rare_variant_results
	
	String sample_id = basename(bam, ".bam")
	
	command{
		java -jar ${picard_jar} SamToFastq I=${bam} FASTQ=sample.fastq
		bwa aln -t ${threads} -o ${max_open_gaps} -n ${missing_alignments_fraction} -l ${seed_length} ${reference} sample.fastq > realigned.sai
		bwa samse ${reference} realigned.sai sample.fastq | samtools view -bS - > realigned.bam
		java -jar ${adna_screen_jar} softclip -b -n ${deamination_bases_to_clip} -i realigned.bam -o clipped_unsorted.bam
		java -jar ${picard_jar} SortSam I=clipped_unsorted.bam O=sorted.bam SORT_ORDER=coordinate
		samtools index sorted.bam
		samtools mpileup -q ${minimum_mapping_quality} -Q ${minimum_base_quality} -f ${reference} sorted.bam | python ${python_calico} --maxdepth 100000 --indels > contamination_rare_variant_results
		python ${python_contamination_rare_variant_results} ${sample_id} contamination_rare_variant_results > contamination_estimate
	}
	output{
		File contamination_estimate = "contamination_estimate"
	}
	runtime{
		cpus: threads
		runtime_minutes: 60
		requested_memory_mb_per_core: 4096
	}
}

task contammix{
	File bam
	File picard_jar
	File htsbox
	File potential_contaminants_fa
	File contammix_estimate
	File python_contammix_multiprocess
	
	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	
	Int threads
	Int copies
	Int chains
	Int seed
	
	# this reference is only for computing the consensus sequence
	# all later alignment for contammix is relative to this consensus sequence
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_pac
	File reference_sa
	File reference_fai
	
	File coverages
	# We downsample to max MT coverage
	Float max_coverage
	
	String sample_id = basename(bam, ".bam")
	
	#samtools mpileup -u -Q ${minimum_base_quality} -q ${minimum_mapping_quality} -f ${reference} ${bam} | bcftools call -c -O z --ploidy 1 -o calls.vcf.gz
	#tabix calls.vcf.gz
	#cat ${reference} | bcftools consensus calls.vcf.gz > consensus.fa
	
	command{
		${htsbox} pileup -f ${reference} -Q ${minimum_base_quality} -q ${minimum_mapping_quality} -M ${bam} > ${sample_id}.consensus.fa
		java -jar ${picard_jar} SamToFastq I=${bam} FASTQ=for_alignment_to_consensus.fastq 
		bwa index ${sample_id}.consensus.fa
		bwa aln -t ${threads} -o ${max_open_gaps} -n ${missing_alignments_fraction} -l ${seed_length} ${sample_id}.consensus.fa for_alignment_to_consensus.fastq > realigned.sai
		bwa samse ${sample_id}.consensus.fa realigned.sai for_alignment_to_consensus.fastq | samtools view -bS - > realigned.bam
		
		python <<CODE
		# read in coverages from file
		coverage_by_sampleID = dict()
		filename = ${coverages}
		with open(filename) as f:
			for line in f:
				sampleID, coverage = line.split('\t')
				coverage_float = float(coverage)
				coverage_by_sampleID[sampleID] = coverage_float
		
		coverage = coverage_by_sampleID[${sample_id}]
		retain_probability = (${max_coverage} / coverage) if (coverage > ${max_coverage}) else 1.0
		subprocess.check_output("java -jar ${picard_jar} DownsampleSam I=realigned.bam O=${sample_id}.downsampled.bam PROBABILITY=%s" % (retain_probability,), shell=True)
		CODE
		
		cat ${sample_id}.consensus.fa ${potential_contaminants_fa} > all_fasta
		mafft all_fasta > multiple_alignment.fa
		python ${python_contammix_multiprocess} ${copies} ${sample_id} ${contammix_estimate} ${sample_id}.downsampled.bam multiple_alignment.fa ${chains} ${minimum_base_quality} ${deamination_bases_to_clip} ${seed} > contamination_estimate
	}
	output{
		File contamination_estimate = "contamination_estimate"
		File consensus = "${sample_id}.consensus.fa"
	}
	runtime{
		cpus: threads
		requested_memory_mb_per_core: 7000
	}
}

task kmer_analysis{
	File python_kmer_analysis
	File barcodes_q_only
	File counts_by_index_barcode_key
	
	String dataset_label
	String date

	command{
		python ${python_kmer_analysis} ${barcodes_q_only} ${counts_by_index_barcode_key} > ${date}_${dataset_label}.kmer
	}
	output{
		File analysis = "${date}_${dataset_label}.kmer"
	}
	runtime{
		runtime_minutes: 60
		requested_memory_mb_per_core: 1000
	}
	
}

task prepare_report{
	File python_prepare_report
	File aggregated_statistics
	File index_barcode_keys
	# these files contain statistics that are not necessarily count based
	# they do not contain a leading number of reads
	Array[File] keyed_statistics
	
	String dataset_label
	String date
	
	command{
		python ${python_prepare_report} ${aggregated_statistics} ${index_barcode_keys} ${sep=' ' keyed_statistics} > ${date}_${dataset_label}.report
	}
	output{
		File report = "${date}_${dataset_label}.report"
	}
	runtime{
		requested_memory_mb_per_core: 2000
	}
}
