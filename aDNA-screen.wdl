workflow ancientDNA_screen{
	String blc_input_directory
	File i5_indices
	File i7_indices
	File barcodeSets
	
	File adna_screen_jar
	File picard_jar
	File pmdtools
	
	Int minimum_mapping_quality
	Int samples_to_demultiplex
	
	File python_damage
	File python_target
	File python_central_measures
	
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_pac
	File reference_sa
	
	File mt_reference
	File mt_reference_amb
	File mt_reference_ann
	File mt_reference_bwt
	File mt_reference_pac
	File mt_reference_sa
	
	String output_path_hs37d5_aligned_unfiltered
	String output_path_hs37d5_aligned_filtered
	String output_path_rsrs_aligned_filtered

	call bcl2fastq { input : blc_input_directory=blc_input_directory} 
	scatter(lane in bcl2fastq.read_files_by_lane){
		call discover_lane_name_from_filename{ input:
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
	scatter(fastq_to_align in collect_filenames.filenames){
		call align as align_hs37d5{ input:
			fastq_to_align = fastq_to_align,
			reference = reference,
			reference_amb = reference_amb,
			reference_ann = reference_ann,
			reference_bwt = reference_bwt,
			reference_pac = reference_pac,
			reference_sa = reference_sa
		}
		call align as align_rsrs{ input:
			fastq_to_align = fastq_to_align,
			reference = mt_reference,
			reference_amb = mt_reference_amb,
			reference_ann = mt_reference_ann,
			reference_bwt = mt_reference_bwt,
			reference_pac = mt_reference_pac,
			reference_sa = mt_reference_sa
		}
	}
	call demultiplex as demultiplex_hs37d5 {input:
		adna_screen_jar = adna_screen_jar,
		prealignment_statistics = aggregate_lane_statistics.statistics,
		aligned_sam_files = align_hs37d5.sam,
		samples_to_demultiplex = samples_to_demultiplex
	}
	call target as hs37d5_target{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = demultiplex_hs37d5.demultiplexed_bam,
		targets="\"{'autosome_pre':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X_pre':'X','Y_pre':'Y'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	scatter(bam in demultiplex_hs37d5.demultiplexed_bam){
		call process_sample as process_sample_hs37d5 { input: 
			picard_jar = picard_jar,
			adna_screen_jar = adna_screen_jar,
			pmdtools = pmdtools,
			unsorted = bam,
			python_damage = python_damage,
			duplicates_label = "duplicates_hs37d5",
			damage_label = "damage_hs37d5"
		}
	}
	call target as hs37d5_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = process_sample_hs37d5.aligned_deduplicated,
		targets="\"{'autosome_post':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X_post':'X','Y_post':'Y'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	
	call demultiplex as demultiplex_rsrs {input:
		adna_screen_jar = adna_screen_jar,
		prealignment_statistics = aggregate_lane_statistics.statistics,
		aligned_sam_files = align_rsrs.sam,
		samples_to_demultiplex = samples_to_demultiplex
	}
	call target as rsrs_target{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = demultiplex_rsrs.demultiplexed_bam,
		targets="\"{'MT_pre':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	scatter(bam in demultiplex_rsrs.demultiplexed_bam){
		call process_sample as process_sample_rsrs{ input: 
			picard_jar = picard_jar,
			adna_screen_jar = adna_screen_jar,
			pmdtools = pmdtools,
			unsorted = bam,
			python_damage = python_damage,
			duplicates_label = "duplicates_rsrs",
			damage_label = "damage_rsrs"
		}
		call haplogrep as haplogrep_rsrs{ input:
			minimum_mapping_quality = minimum_mapping_quality,
			region = "MT",
			bam = process_sample_rsrs.aligned_deduplicated,
			reference = mt_reference,
			reference_amb = mt_reference_amb,
			reference_ann = mt_reference_ann,
			reference_bwt = mt_reference_bwt,
			reference_pac = mt_reference_pac,
			reference_sa = mt_reference_sa
		}
	}
	call target as rsrs_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = process_sample_rsrs.aligned_deduplicated,
		targets="\"{'MT_post':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	
	call central_measures as central_measures_hs37d5{ input:
		python_central_measures = python_central_measures,
		mean_label = "mean_hs37d5",
		median_label = "median_hs37d5",
		histograms = hs37d5_target_post.length_histogram
	}
	call central_measures as central_measures_rsrs{ input:
		python_central_measures = python_central_measures,
		mean_label = "mean_rsrs",
		median_label = "median_rsrs",
		histograms = rsrs_target_post.length_histogram
	}
	call summarize_haplogroups{ input:
		haplogrep_output = haplogrep_rsrs.haplogroup_report
	}

	call aggregate_statistics as aggregate_statistics_duplicates_hs37d5{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = process_sample_hs37d5.duplicates_statistics
	}
	call aggregate_statistics as aggregate_statistics_duplicates_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = process_sample_rsrs.duplicates_statistics
	}
	
	call aggregate_statistics as aggregate_statistics_pre_hs37d5{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = hs37d5_target.target_stats
	}
	call aggregate_statistics as aggregate_statistics_post_hs37d5{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = hs37d5_target_post.target_stats
	}
	call aggregate_statistics as aggregate_statistics_pre_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = rsrs_target.target_stats
	}
	call aggregate_statistics as aggregate_statistics_post_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = rsrs_target_post.target_stats
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
		files = process_sample_hs37d5.aligned_deduplicated,
		output_path = output_path_hs37d5_aligned_filtered
	}
	call copy_output as copy_hs37d5_histogram{ input:
		files = hs37d5_target_post.length_histogram,
		output_path = output_path_hs37d5_aligned_filtered
	}
	call copy_output as copy_rsrs_aligned_filtered{ input:
		files = process_sample_rsrs.aligned_deduplicated,
		output_path = output_path_rsrs_aligned_filtered
	}
	call copy_output as copy_rsrs_histogram{ input:
		files = rsrs_target_post.length_histogram,
		output_path = output_path_rsrs_aligned_filtered
	}
	call concatenate as concatenate_rsrs_damage{ input:
		to_concatenate = process_sample_rsrs.damage
	}
	call concatenate as concatenate_hs37d5_damage{ input:
		to_concatenate = process_sample_hs37d5.damage
	}
	
	Array[File] final_keyed_statistics = [
		concatenate_hs37d5_damage.concatenated,
		concatenate_rsrs_damage.concatenated,
		central_measures_hs37d5.central_measures_output,
		central_measures_rsrs.central_measures_output,
		summarize_haplogroups.haplogroups
		
	]
	call prepare_report{ input:
		aggregated_statistics = aggregate_statistics_final.statistics,
		keyed_statistics = final_keyed_statistics
	}
}

task bcl2fastq{
	String blc_input_directory

	command{
		bcl2fastq \
			-R ${blc_input_directory} \
			-o ./ \
			--create-fastq-for-index-reads \
			--use-bases-mask Y76,I7,I7,Y76
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
		queue: "mcore"
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
			cpus: 1
			runtime_minutes: 10
			requested_memory_mb_per_core: 2048
			queue: "short"
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
		java -Xmx14g -jar ${adna_screen_jar} ${i5_indices} ${i7_indices} ${barcodeSets} ${read_files_by_lane[0]} ${read_files_by_lane[1]} ${read_files_by_lane[2]} ${read_files_by_lane[3]} ${label} > ${label}.stats
	}
	
	output{
		Array[File] fastq_to_align = glob("${label}*.fastq.gz")
		File statistics = "${label}.stats"
	}
	runtime{
			cpus: 1
			runtime_minutes: 720
			requested_memory_mb_per_core: 16384
			queue: "short"
	}
}

task aggregate_statistics{
	File adna_screen_jar
	Array [File] statistics_by_group
	
	command{
		java -cp ${adna_screen_jar} adnascreen.AggregateStatistics ${sep=' ' statistics_by_group} > aggregated_statistics
	}
	output{
		File statistics = "aggregated_statistics"
	}
	runtime{
			cpus: 1
			runtime_minutes: 20
			requested_memory_mb_per_core: 4096
			queue: "short"
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
	
	command {
		bwa aln -t ${threads} -o ${max_open_gaps} -n ${missing_alignments_fraction} -l ${seed_length} ${reference} ${fastq_to_align} > aligned.sai && \
		bwa samse ${reference} aligned.sai ${fastq_to_align} > aligned.sam
	}
	output{
		File sam = "aligned.sam"
	}
	runtime{
			cpus: "${threads}"
			runtime_minutes: 600
			requested_memory_mb_per_core: 8192
			queue: "mcore"
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
			cpus: 1
			runtime_minutes: 10
			requested_memory_mb_per_core: 2048
			queue: "short"
	}
}

task demultiplex{
	File adna_screen_jar
	File prealignment_statistics
	Array[File] aligned_sam_files
	Int samples_to_demultiplex
	
	command{
		java -Xmx14g -cp ${adna_screen_jar} adnascreen.DemultiplexSAM -b -n ${samples_to_demultiplex} -s ${prealignment_statistics} ${sep=' ' aligned_sam_files} > postalignment_statistics
	}
	output{
		Array[File] demultiplexed_bam = glob("*.bam")
		File statistics = "postalignment_statistics"
	}
	runtime{
			cpus: 2
			runtime_minutes: 600
			requested_memory_mb_per_core: 8000
			queue: "short"
	}
}

# filter out unaligned and duplicate reads
# compute damage
task process_sample{
	File picard_jar
	File adna_screen_jar
	File pmdtools
	File unsorted
	File python_damage
	String duplicates_label
	String damage_label
	
	String sample_id_filename = sub(unsorted, ".*/", "") # remove leading directories from full path to leave only filename
	
	command{
		set -e
		java -jar ${picard_jar} SortSam I=${unsorted} O=sorted_queryname.bam SORT_ORDER=queryname
		java -jar ${picard_jar} FilterSamReads I=sorted_queryname.bam O=filtered.bam FILTER=includeAligned 
		java -jar ${picard_jar} SortSam I=filtered.bam O=sorted_coordinate.bam SORT_ORDER=coordinate
		java -jar ${picard_jar} MarkDuplicates I=sorted_coordinate.bam O=${sample_id_filename} M=${sample_id_filename}.dedup_stats REMOVE_DUPLICATES=true
		java -cp ${adna_screen_jar} adnascreen.ReadMarkDuplicatesStatistics -l ${duplicates_label} ${sample_id_filename}.dedup_stats > ${sample_id_filename}.stats
		
		echo "${sample_id_filename}" > damage
		java -jar ${picard_jar} ViewSam INPUT=${sample_id_filename} ALIGNMENT_STATUS=Aligned | python ${pmdtools} --first >> damage
		python ${python_damage} ${damage_label} damage > damage_statistics
	}
	output{
		String id = sample_id_filename
		File aligned_deduplicated = "${sample_id_filename}"
		File duplicates_statistics = "${sample_id_filename}.stats"
		File damage = "damage_statistics"
	}
}

# This runs too quickly to run on Orchestra short queue
# The best solution would be to run a loop within this task, but this is not yet supported
task target_forloop{
	File adna_screen_jar
	File bam
	String targets
	Int minimum_mapping_quality
	
	String sample_id_filename = sub(bam, ".*/", "") # remove leading directories from full path to leave only filename

	command{
		java -cp ${adna_screen_jar} adnascreen.SAMStats -f ${bam} -t ${targets} -l ${sample_id_filename}.histogram -q ${minimum_mapping_quality} > ${sample_id_filename}.stats
	}
	output{
		File target_stats = "${sample_id_filename}.stats"
		File length_histogram = "${sample_id_filename}.histogram"
	}
	runtime{
			cpus: 1
			runtime_minutes: 30
			requested_memory_mb_per_core: 4096
			queue: "short"
	}
}

# Alternative in place of looping in WDL, run loop in python
task target{
	File python_target
	File adna_screen_jar
	Array[File] bams
	String targets
	Int minimum_mapping_quality
	
	#String sample_id_filename = sub(bam, ".*/", "") # remove leading directories from full path to leave only filename

	command{
		python ${python_target} ${adna_screen_jar} ${targets} ${minimum_mapping_quality} ${sep=' ' bams}
	}
	output{
		Array[File] target_stats = glob("*.stats")
		Array[File] length_histogram = glob("*.histogram")
	}
	runtime{
			cpus: 2
			runtime_minutes: 180
			requested_memory_mb_per_core: 8192
			queue: "short"
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
			cpus: 1
			runtime_minutes: 30
			requested_memory_mb_per_core: 4096
			queue: "short"
	}
}

task copy_output{
	Array[File] files
	String output_path
	
	command{
		mkdir -p ${output_path};
		for file in ${sep=' ' files}  ; do 
			cp -l $file "${output_path}"
		done
	}
	runtime{
			cpus: 1
			runtime_minutes: 300
			requested_memory_mb_per_core: 2048
			queue: "short"
	}
}

task haplogrep{
	Int minimum_mapping_quality
	Int minimum_base_quality
	String region
	File bam
	File haplogrep_jar
	Int phylotree_version
	
	String sample_id_filename = sub(bam, ".*/", "") # remove leading directories from full path to leave only filename
	
	# value from samtools for bwa
	Int excessive_mismatch_penalty = 50
	
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_pac
	File reference_sa

	command{
		set -e
		samtools index ${bam}
samtools mpileup -q ${minimum_mapping_quality} -Q ${minimum_base_quality} -C ${excessive_mismatch_penalty} -r ${region} -u -f ${reference} ${bam} | bcftools call -m -v > ${sample_id_filename}.vcf
		java -jar ${haplogrep_jar} --format vcf --phylotree ${phylotree_version} --in ${sample_id_filename}.vcf --out ${sample_id_filename}.haplogroup
	}
	output{
		File haplogroup_report = "${sample_id_filename}.haplogroup"
	}
	runtime{
		cpus: 1
		runtime_minutes: 30
		requested_memory_mb_per_core: 8192
		queue: "short"
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
		cpus: 1
		runtime_minutes: 30
		requested_memory_mb_per_core: 4096
		queue: "short"
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
		cpus: 1
		runtime_minutes: 30
		requested_memory_mb_per_core: 4096
		queue: "short"
	}
}

task prepare_report{
	File python_prepare_report
	File aggregated_statistics
	# these files contain statistics that are not necessarily count based
	# they do not contain a leading number of reads
	Array[File] keyed_statistics
	
	command{
		python ${python_prepare_report} ${aggregated_statistics} ${sep=' ' keyed_statistics} > report
	}
	output{
		File report = "report"
	}
	runtime{
		cpus: 1
		runtime_minutes: 30
		requested_memory_mb_per_core: 4096
		queue: "short"
	}
}
