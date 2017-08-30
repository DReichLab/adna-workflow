workflow ancientDNA_screen{
	String blc_input_directory
	String dataset_label
	String date
	
	File i5_indices
	File i7_indices
	File barcodeSets
	
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
	File python_target
	File python_central_measures
	File python_snp_target
	File python_snp_target_bed
	File python_coverage
	File python_floor
	
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
			label = discover_lane_name_from_filename.lane
		}
	}
	call aggregate_statistics as aggregate_lane_statistics{ input :
		adna_screen_jar=adna_screen_jar,
		statistics_by_group=merge_and_trim_lane.statistics
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
	scatter(bam in demultiplex_hs37d5.demultiplexed_bam){
		call filter_aligned_only as filter_aligned_only_hs37d5 { input:
			#picard_jar = picard_jar,
			bam = bam			
		}
		call snp_target_bed as spike3k_pre{ input:
			coordinates_autosome = spike3k_coordinates_autosome,
			coordinates_x = spike3k_coordinates_x,
			coordinates_y = spike3k_coordinates_y,
			bam = filter_aligned_only_hs37d5.filtered,
			minimum_mapping_quality = minimum_mapping_quality,
			minimum_base_quality = minimum_base_quality,
			deamination_bases_to_clip = deamination_bases_to_clip,
			label = "spike3k_pre",
			picard_jar = picard_jar,
			adna_screen_jar = adna_screen_jar,
			python_snp_target_bed = python_snp_target_bed,
			reference = prepare_reference_hs37d5.reference_fa,
			reference_amb = prepare_reference_hs37d5.reference_amb,
			reference_ann = prepare_reference_hs37d5.reference_ann,
			reference_bwt = prepare_reference_hs37d5.reference_bwt,
			reference_pac = prepare_reference_hs37d5.reference_pac,
			reference_sa = prepare_reference_hs37d5.reference_sa,
			reference_fai = prepare_reference_hs37d5.reference_fai
		}
		call duplicates_and_damage as duplicates_and_damage_hs37d5 { input: 
			picard_jar = picard_jar,
			adna_screen_jar = adna_screen_jar,
			pmdtools = pmdtools,
			unsorted = filter_aligned_only_hs37d5.filtered,
			python_damage = python_damage,
			duplicates_label = "duplicates_hs37d5",
			damage_label = "damage_hs37d5"
		}
		call snp_target_bed as spike3k_post{ input:
			coordinates_autosome = spike3k_coordinates_autosome,
			coordinates_x = spike3k_coordinates_x,
			coordinates_y = spike3k_coordinates_y,
			bam = duplicates_and_damage_hs37d5.aligned_deduplicated,
			minimum_mapping_quality = minimum_mapping_quality,
			minimum_base_quality = minimum_base_quality,
			deamination_bases_to_clip = deamination_bases_to_clip,
			label = "spike3k_post",
			picard_jar = picard_jar,
			adna_screen_jar = adna_screen_jar,
			python_snp_target_bed = python_snp_target_bed,
			reference = prepare_reference_hs37d5.reference_fa,
			reference_amb = prepare_reference_hs37d5.reference_amb,
			reference_ann = prepare_reference_hs37d5.reference_ann,
			reference_bwt = prepare_reference_hs37d5.reference_bwt,
			reference_pac = prepare_reference_hs37d5.reference_pac,
			reference_sa = prepare_reference_hs37d5.reference_sa,
			reference_fai = prepare_reference_hs37d5.reference_fai
		} 
	}
	call chromosome_target as hs37d5_chromosome_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = duplicates_and_damage_hs37d5.aligned_deduplicated,
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
	scatter(bam in demultiplex_rsrs.demultiplexed_bam){
		call filter_aligned_only as filter_aligned_only_rsrs{ input:
			#picard_jar = picard_jar,
			bam = bam
		}
		call duplicates_and_damage as duplicates_and_damage_rsrs{ input: 
			picard_jar = picard_jar,
			adna_screen_jar = adna_screen_jar,
			pmdtools = pmdtools,
			unsorted = filter_aligned_only_rsrs.filtered,
			python_damage = python_damage,
			duplicates_label = "duplicates_rsrs",
			damage_label = "damage_rsrs"
		}
		call haplogrep as haplogrep_rcrs{ input:
			missing_alignments_fraction = missing_alignments_fraction,
			max_open_gaps = max_open_gaps,
			seed_length = seed_length,
			minimum_mapping_quality = minimum_mapping_quality,
			minimum_base_quality = minimum_base_quality,
			deamination_bases_to_clip = deamination_bases_to_clip,
			bam = duplicates_and_damage_rsrs.aligned_deduplicated,
			reference = prepare_reference_rcrs.reference_fa,
			reference_amb = prepare_reference_rcrs.reference_amb,
			reference_ann = prepare_reference_rcrs.reference_ann,
			reference_bwt = prepare_reference_rcrs.reference_bwt,
			reference_pac = prepare_reference_rcrs.reference_pac,
			reference_sa = prepare_reference_rcrs.reference_sa,
			adna_screen_jar = adna_screen_jar,
			picard_jar = picard_jar,
			coverage_int = chromosome_target_single_rsrs.coverage_int
		}
		call chromosome_target_single as chromosome_target_single_rsrs{ input:
			adna_screen_jar = adna_screen_jar,
			python_coverage = python_coverage,
			python_floor = python_floor,
			bam = duplicates_and_damage_rsrs.aligned_deduplicated,
			targets="\"{'MT_post':'MT'}\"",
			minimum_mapping_quality = minimum_mapping_quality,
			reference_length = 16569,
			coverage_field = "MT_post-coverageLength"
		}
#		call schmutzi{ input:
#			bam = duplicates_and_damage_rsrs.aligned_deduplicated,
#			picard_jar = picard_jar,
#			reference = prepare_reference_rsrs.reference_fa,
#			reference_amb = prepare_reference_rsrs.reference_amb,
#			reference_ann = prepare_reference_rsrs.reference_ann,
#			reference_bwt = prepare_reference_rsrs.reference_bwt,
#			reference_pac = prepare_reference_rsrs.reference_pac,
#			reference_sa = prepare_reference_rsrs.reference_sa,
#			reference_fai = prepare_reference_rsrs.reference_fai,
#			coverage = chromosome_target_single_rsrs.coverage
#		}
		call contamination_rare_variant{ input:
			bam = duplicates_and_damage_rsrs.aligned_deduplicated,
			picard_jar = picard_jar,
			adna_screen_jar = adna_screen_jar,
			missing_alignments_fraction = missing_alignments_fraction,
			max_open_gaps = max_open_gaps,
			seed_length = seed_length,
			minimum_mapping_quality = minimum_mapping_quality,
			minimum_base_quality = minimum_base_quality,
			deamination_bases_to_clip = deamination_bases_to_clip,
			reference = prepare_reference_human_95_consensus.reference_fa,
			reference_amb = prepare_reference_human_95_consensus.reference_amb,
			reference_ann = prepare_reference_human_95_consensus.reference_ann,
			reference_bwt = prepare_reference_human_95_consensus.reference_bwt,
			reference_pac = prepare_reference_human_95_consensus.reference_pac,
			reference_sa = prepare_reference_human_95_consensus.reference_sa,
			reference_fai = prepare_reference_human_95_consensus.reference_fai
		}
		call contammix{ input:
			bam = duplicates_and_damage_rsrs.aligned_deduplicated,
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
			coverage = chromosome_target_single_rsrs.coverage
		}
	}
	call chromosome_target as rsrs_chromosome_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = duplicates_and_damage_rsrs.aligned_deduplicated,
		targets="\"{'MT_post':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
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
		statistics_by_group = duplicates_and_damage_hs37d5.duplicates_statistics
	}
	call aggregate_statistics as aggregate_statistics_duplicates_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = duplicates_and_damage_rsrs.duplicates_statistics
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
		files = duplicates_and_damage_hs37d5.aligned_deduplicated,
		output_path = output_path_hs37d5_aligned_filtered
	}
	call copy_output as copy_hs37d5_histogram{ input:
		files = hs37d5_chromosome_target_post.length_histogram,
		output_path = output_path_hs37d5_aligned_filtered
	}
	call copy_output as copy_rsrs_aligned_filtered{ input:
		files = duplicates_and_damage_rsrs.aligned_deduplicated,
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
	call concatenate as concatenate_rsrs_damage{ input:
		to_concatenate = duplicates_and_damage_rsrs.damage
	}
	call concatenate as concatenate_hs37d5_damage{ input:
		to_concatenate = duplicates_and_damage_hs37d5.damage
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
	call concatenate as concatenate_contamination_rare_variant{ input:
		to_concatenate = contamination_rare_variant.contamination_estimate
	}
	call concatenate as concatenate_contammix{ input:
		to_concatenate = contammix.contamination_estimate
	}
	call spike3k_complexity{ input:
		spike3k_pre_data = concatenate_spike3k_pre.concatenated,
		spike3k_post_data =concatenate_spike3k_post.concatenated
	}
	
	Array[File] final_keyed_statistics = [
		concatenate_hs37d5_damage.concatenated,
		concatenate_rsrs_damage.concatenated,
		central_measures_hs37d5.central_measures_output,
		central_measures_rsrs.central_measures_output,
		summarize_haplogroups.haplogroups,
		concatenate_spike3k_pre.concatenated,
		concatenate_spike3k_post.concatenated,
		spike3k_complexity.estimates,
#		concatenate_schmutzi.concatenated,
		concatenate_contamination_rare_variant.concatenated,
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
		bcl2fastq \
			-R ${blc_input_directory} \
			-o ./ \
			--create-fastq-for-index-reads
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
		runtime_minutes: 240
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
		requested_memory_mb_per_core: 2048
	}
}

task merge_and_trim_lane{
	File adna_screen_jar
	File i5_indices
	File i7_indices
	File barcodeSets
	Array[File] read_files_by_lane
	String label
	command{
		java -Xmx14g -jar ${adna_screen_jar} IndexAndBarcodeScreener --i5-indices ${i5_indices} --i7-indices ${i7_indices} --barcodes ${barcodeSets} ${read_files_by_lane[0]} ${read_files_by_lane[1]} ${read_files_by_lane[2]} ${read_files_by_lane[3]} ${label} > ${label}.stats
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
		runtime_minutes: 20
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
		runtime_minutes: 10
		requested_memory_mb_per_core: 2048
	}
}

task demultiplex{
	File adna_screen_jar
	File prealignment_statistics
	Array[File] aligned_bam_files
	Int samples_to_demultiplex
	
	command{
		java -Xmx14g -jar ${adna_screen_jar} DemultiplexSAM -b -n ${samples_to_demultiplex} -s ${prealignment_statistics} ${sep=' ' aligned_bam_files} > postalignment_statistics
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
	#File picard_jar
	File bam
	
	String filename = basename(bam)
	
	# picard complains "MAPQ should be 0 for unmapped read." while trying to filter unmapped reads
	#java -jar ${picard_jar} SortSam I=${bam} O=sorted_queryname.bam SORT_ORDER=queryname
	#java -jar ${picard_jar} FilterSamReads I=sorted_queryname.bam O=${filename} FILTER=includeAligned
	command{
		samtools view -h -b -F 4 -o ${filename} ${bam}
	}
	output{
		File filtered = filename
	}
}

task duplicates_and_damage{
	File picard_jar
	File adna_screen_jar
	File pmdtools
	File unsorted
	File python_damage
	String duplicates_label
	String damage_label
	
	String sample_id_filename = basename(unsorted)
	
	command{
		set -e
		java -jar ${picard_jar} SortSam I=${unsorted} O=sorted_coordinate.bam SORT_ORDER=coordinate
		java -jar ${picard_jar} MarkDuplicates I=sorted_coordinate.bam O=${sample_id_filename} M=${sample_id_filename}.dedup_stats REMOVE_DUPLICATES=true BARCODE_TAG=XD ADD_PG_TAG_TO_READS=false
		java -jar ${adna_screen_jar} ReadMarkDuplicatesStatistics -l ${duplicates_label} ${sample_id_filename}.dedup_stats > ${sample_id_filename}.stats
		
		echo "${sample_id_filename}" > damage
		samtools view ${sample_id_filename} | python ${pmdtools} --first >> damage
		python ${python_damage} ${damage_label} damage > damage_statistics
	}
	output{
		String id = sample_id_filename
		File aligned_deduplicated = "${sample_id_filename}"
		File duplicates_statistics = "${sample_id_filename}.stats"
		File damage = "damage_statistics"
	}
	runtime{
		cpus: 2
	}
}

# This runs too quickly to run on Orchestra short queue
# The best solution would be to run a loop within this task, but this is not yet supported
task chromosome_target_single{
	File adna_screen_jar
	File python_coverage
	File python_floor
	File bam
	String targets
	Int minimum_mapping_quality
	Int reference_length
	String coverage_field
	
	String sample_id = basename(bam, ".bam")

	command{
		set -e
		java -jar ${adna_screen_jar} SAMStats -f ${bam} -t ${targets} -l ${sample_id}.histogram -q ${minimum_mapping_quality} > ${sample_id}.stats
		python ${python_coverage} ${sample_id}.stats ${reference_length} ${sample_id} ${coverage_field} > coverage
		python ${python_floor} coverage > coverage_int
	}
	output{
		File target_stats = "${sample_id}.stats"
		File length_histogram = "${sample_id}.histogram"
		Float coverage = read_float("coverage")
		Int coverage_int = read_int("coverage_int")
	}
	runtime{
		runtime_minutes: 30
		requested_memory_mb_per_core: 4096
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
		runtime_minutes: 300
		requested_memory_mb_per_core: 8192
	}
}

task snp_target{
	File coordinates
	File bam
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	String label
	File picard_jar
	File adna_screen_jar
	
	File python_snp_target
		
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_fai
	File reference_pac
	File reference_sa
	
	String sample_id_filename = basename(bam)

	#${htsbox} pileup -vcf ${reference} -q ${minimum_mapping_quality} -Q ${minimum_base_quality} -T ${deamination_bases_to_clip} -b ${coordinates} ${sample_id_filename}.bam > ${sample_id_filename}.vcf
	
	command{
		set -e
		java -jar ${adna_screen_jar} softclip -b -n ${deamination_bases_to_clip} -i ${bam} -o clipped_unsorted.bam
		java -jar ${picard_jar} SortSam I=clipped_unsorted.bam O=sorted.bam SORT_ORDER=coordinate
		samtools index sorted.bam
		samtools mpileup -q ${minimum_mapping_quality} -Q ${minimum_base_quality} -v -u -f ${reference} -l ${coordinates} sorted.bam > ${sample_id_filename}.vcf
		python ${python_snp_target} ${label} ${sample_id_filename}.vcf > snp_target_stats
	}
	output{
		File snp_target_stats = "snp_target_stats"
	}
}

# This counts reads similarly to snp_target, but uses a bed file to look for read overlaps in an interval around the SNP. 
task snp_target_bed{
	File coordinates_autosome
	File coordinates_x
	File coordinates_y
	File bam
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	String label
	File picard_jar
	File adna_screen_jar
	
	File python_snp_target_bed
		
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_fai
	File reference_pac
	File reference_sa
	
	String sample_id_filename = basename(bam)
	
	command{
		set -e
		java -jar ${adna_screen_jar} softclip -b -n ${deamination_bases_to_clip} -i ${bam} -o clipped_unsorted.bam
		java -jar ${picard_jar} SortSam I=clipped_unsorted.bam O=sorted.bam SORT_ORDER=coordinate
		samtools index sorted.bam
		samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_autosome} sorted.bam > ${sample_id_filename}.autosome
		samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_x}        sorted.bam > ${sample_id_filename}.x
		samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_y}        sorted.bam > ${sample_id_filename}.y
		python ${python_snp_target_bed} ${label} ${sample_id_filename}.autosome ${sample_id_filename}.x ${sample_id_filename}.y > snp_target_stats
	}
	output{
		File snp_target_stats = "snp_target_stats"
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
		runtime_minutes: 300
		requested_memory_mb_per_core: 2048
	}
}

task haplogrep{
	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	
	Int minimum_pileup_depth
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	String region = "MT"
	File bam
	File haplogrep_jar
	Int phylotree_version
	File adna_screen_jar
	File picard_jar
	
	Int coverage_int
	Int pileup_depth_from_coverage = coverage_int / 10
	Int pileup_depth = if(pileup_depth_from_coverage > minimum_pileup_depth) then pileup_depth_from_coverage else minimum_pileup_depth
	
	String sample_id_filename = basename(bam, ".bam")
	
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_pac
	File reference_sa
	
	command{
		set -e
		java -jar ${picard_jar} SamToFastq I=${bam} FASTQ=for_alignment_to_rcrs.fastq 
		bwa aln -t 2 -o ${max_open_gaps} -n ${missing_alignments_fraction} -l ${seed_length} ${reference} for_alignment_to_rcrs.fastq > realigned.sai
		bwa samse ${reference} realigned.sai for_alignment_to_rcrs.fastq | samtools view -bS - > realigned.bam
		java -jar ${adna_screen_jar} softclip -b -n ${deamination_bases_to_clip} -i realigned.bam -o clipped_unsorted_realigned.bam
		java -jar ${picard_jar} SortSam I=clipped_unsorted_realigned.bam O=${sample_id_filename}.bam SORT_ORDER=coordinate
		samtools index ${sample_id_filename}.bam
		samtools mpileup -q ${minimum_mapping_quality} -Q ${minimum_base_quality} -r ${region} -u -f ${reference} ${sample_id_filename}.bam | bcftools call -c -v --ploidy 1 > ${sample_id_filename}.vcf
		java -jar ${haplogrep_jar} --format vcf --phylotree ${phylotree_version} --in ${sample_id_filename}.vcf --out ${sample_id_filename}.haplogroup
	}
	output{
		File haplogroup_report = "${sample_id_filename}.haplogroup"
	}
	runtime{
		cpus: 2
		runtime_minutes: 180
		requested_memory_mb_per_core: 8192
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
		runtime_minutes: 30
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
		runtime_minutes: 30
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
	Int threads = 2

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
	
	Float coverage
	# We downsample to 600x MT coverage
	Float threshold = 600.0
	Float retain_probability = if (coverage > threshold) then (threshold / coverage) else 1.0
	
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
		java -jar ${picard_jar} DownsampleSam I=realigned.bam O=${sample_id}.downsampled.bam PROBABILITY=${retain_probability}
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
	}
}

task prepare_report{
	File python_prepare_report
	File aggregated_statistics
	# these files contain statistics that are not necessarily count based
	# they do not contain a leading number of reads
	Array[File] keyed_statistics
	
	String dataset_label
	String date
	
	command{
		python ${python_prepare_report} ${aggregated_statistics} ${sep=' ' keyed_statistics} > ${date}_${dataset_label}.report
	}
	output{
		File report = "${date}_${dataset_label}.report"
	}
	runtime{
		runtime_minutes: 30
		requested_memory_mb_per_core: 4096
	}
}
