import "demultiplex.wdl" as demultiplex_align_bams

workflow adna_analysis{
	File adna_screen_jar
	File picard_jar
	File pmdtools
	File htsbox
	
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
	File python_snp_target
	File python_snp_target_bed
	File python_coverage
	File python_floor
	
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
	
	call demultiplex_align_bams.prepare_reference as prepare_reference_rsrs{ input:
		reference = mt_reference_rsrs_in
	}
	call demultiplex_align_bams.prepare_reference as prepare_reference_rcrs{ input:
		reference = mt_reference_rcrs_in
	}
	call demultiplex_align_bams.prepare_reference as prepare_reference_human_95_consensus{ input:
		reference = mt_reference_human_95_consensus
	}
	
	call chromosome_target as hs37d5_chromosome_target{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = demultiplex_hs37d5.demultiplexed_bam,
		targets="\"{'autosome_pre':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X_pre':'X','Y_pre':'Y','human_pre':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']}\"",
		minimum_mapping_quality = minimum_mapping_quality
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
	call chromosome_target as rsrs_chromosome_target{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = demultiplex_rsrs.demultiplexed_bam,
		targets="\"{'MT_pre':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
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
			bam = bam,
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

	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_duplicates_hs37d5{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = duplicates_hs37d5.duplicates_statistics
	}
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_duplicates_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = duplicates_rsrs.duplicates_statistics
	}
	
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_pre_hs37d5{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = hs37d5_chromosome_target.target_stats
	}
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_post_hs37d5{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = hs37d5_chromosome_target_post.target_stats
	}
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_pre_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = rsrs_chromosome_target.target_stats
	}
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_post_rsrs{ input:
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
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_final{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = cumulative_statistics
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
	
	# same as final statistics, but missing the expensive contammix calculation
	Array[File] preliminary_keyed_statistics = [
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
#		concatenate_contammix.concatenated
	]
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
	call prepare_report as preliminary_report{ input:
		aggregated_statistics = aggregate_statistics_final.statistics,
		keyed_statistics = preliminary_keyed_statistics,
		index_barcode_keys = index_barcode_keys,
		dataset_label = dataset_label,
		date = date
	}
	call prepare_report{ input:
		aggregated_statistics = aggregate_statistics_final.statistics,
		keyed_statistics = final_keyed_statistics,
		index_barcode_keys = index_barcode_keys,
		dataset_label = dataset_label,
		date = date
	}
	Array[File] report_array = [prepare_report.report]
	call demultiplex_align_bams.copy_output as copy_report{ input:
		files = report_array,
		output_path = output_path
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
		mkdir deduplicated
		
		python <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def deduplicate_bam(bam):
			sample_id_filename = basename(bam)
			sample_id_filename_no_extension, extension = splitext(sample_id_filename)
			
			sorted_bam = sample_id_filename_no_extension + ".sorted_coordinate.bam"
			
			subprocess.check_output("java -Xmx9g -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (bam, sorted_bam), shell=True)
			subprocess.check_output("java -Xmx9g -jar ${picard_jar} MarkDuplicates I=%s O=deduplicated/%s M=%s.dedup_stats REMOVE_DUPLICATES=true BARCODE_TAG=XD ADD_PG_TAG_TO_READS=false MAX_FILE_HANDLES=1000" % (sorted_bam, sample_id_filename, sample_id_filename), shell=True)
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
		Array[File] aligned_deduplicated = glob("deduplicated/*.bam")
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
		results = [pool.apply_async(damage_for_bam, args=(bam,)) for bam in bams]
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
		for bam_stat in ${sep=' ' bam_stats}
		do
			sample_id=$(basename $bam_stat .bam.stats)
			python ${python_coverage} $bam_stat ${reference_length} $sample_id ${coverage_field} >> coverages
		done
	}
	output{
		File coverages = "coverages"
	}
	runtime{
		runtime_minutes: 120
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

# measure the number of reads at each 1240k site
task target_depth_bed{
	File coordinates_autosome
	File coordinates_x
	File coordinates_y
	Array[File] bams
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	File picard_jar
	File adna_screen_jar
	
	File python_depth_histogram
	File python_depth_histogram_combine
	
	# depth requires qualities greater than arguments, so we adjust them to preserve external minimum quality usage
	Int base_quality = minimum_base_quality-1
	Int mapping_quality = minimum_mapping_quality-1
	
	Int processes = 4
	
	# samtools depth -b ${targets_bed} -q ${base_quality} -Q ${mapping_quality} $bam > depths
	command{
		set -e
		
		python <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def calculate_SNP_depths(bam):
			sample_id_filename = basename(bam)
			sample_id_filename_no_extension, extension = splitext(sample_id_filename)
			
			clipped_bam = sample_id_filename_no_extension + ".clipped.bam"
			clipped_sorted_bam = sample_id_filename_no_extension + ".clipped.sorted.bam"
			
			subprocess.check_output("java -Xmx2500m -jar ${adna_screen_jar} softclip -b -n ${deamination_bases_to_clip} -i %s -o %s" % (bam, clipped_bam), shell=True)
			subprocess.check_output("java -Xmx2500m -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (clipped_bam, clipped_sorted_bam), shell=True)
			subprocess.check_output("samtools index %s" % (clipped_sorted_bam,), shell=True)
			subprocess.check_output("samtools depth -b ${coordinates_autosome} -q ${base_quality} -Q ${mapping_quality} %s | python ${python_depth_histogram} > %s.autosome" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools depth -b ${coordinates_x}        -q ${base_quality} -Q ${mapping_quality} %s | python ${python_depth_histogram} > %s.x" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools depth -b ${coordinates_y}        -q ${base_quality} -Q ${mapping_quality} %s | python ${python_depth_histogram} > %s.y" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("python ${python_depth_histogram_combine} autosome %s.autosome X %s.x Y %s.y > %s.depth_histogram" % (sample_id_filename, sample_id_filename, sample_id_filename, sample_id_filename), shell=True)
			
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		[pool.apply_async(calculate_SNP_depths, args=(bam,)) for bam in bams]
		pool.close()
		pool.join()
		CODE
	}
	output{
		Array[File] depth_histograms = glob("*.depth_histogram")
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
		runtime_minutes: 60
		requested_memory_mb_per_core: 2000
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
		runtime_minutes: 60
		requested_memory_mb_per_core: 2000
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
		import subprocess
		
		# read in coverages from file
		coverage_by_sampleID = dict()
		filename = '${coverages}'
		with open(filename) as f:
			for line in f:
				sampleID, coverage = line.split('\t')
				coverage_float = float(coverage)
				coverage_by_sampleID[sampleID] = coverage_float
		
		coverage = coverage_by_sampleID['${sample_id}']
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
		runtime_minutes: 480
		requested_memory_mb_per_core: 7000
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
		runtime_minutes: 60
		requested_memory_mb_per_core: 2000
	}
}