import "demultiplex.wdl" as demultiplex_align_bams

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
	File angsd
	
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
	
	call analysis_versions as versions{ input:
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		htsbox = htsbox,
		haplogrep_jar = haplogrep_jar,
		angsd = angsd,
		index_barcode_keys_to_stop_call_caching = index_barcode_keys
	}
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_across_sequencing_runs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = demultiplex_statistics_files
	}
	
	call combine_bams_into_libraries as combine_nuclear_libraries{ input:
		bam_lists = nuclear_bam_lists_to_merge,
		picard_jar = picard_jar
	}
	call combine_bams_into_libraries as combine_mt_libraries{ input:
		bam_lists = mt_bam_lists_to_merge,
		picard_jar = picard_jar
	}
	
	call preseq{ input:
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
	call chromosome_target as nuclear_chromosome_target{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = combine_nuclear_libraries.library_bams,
		targets="\"{'autosome_pre':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X_pre':'X','Y_pre':'Y','human_pre':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	call snp_target_bed as count_1240k_pre{ input:
		coordinates_autosome = coordinates_1240k_autosome,
		coordinates_x = coordinates_1240k_x,
		coordinates_y = coordinates_1240k_y,
		bams = combine_nuclear_libraries.library_bams,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip,
		label = "1240k_pre",
		picard_jar = picard_jar,
		adna_screen_jar = adna_screen_jar,
		python_snp_target_bed = python_snp_target_bed
	}
	call duplicates as duplicates_nuclear { input: 
		picard_jar = picard_jar,
		adna_screen_jar = adna_screen_jar,
		unsorted = combine_nuclear_libraries.library_bams,
		duplicates_label = "duplicates_nuclear"
	}
	call snp_target_bed as count_1240k_post{ input:
		coordinates_autosome = coordinates_1240k_autosome,
		coordinates_x = coordinates_1240k_x,
		coordinates_y = coordinates_1240k_y,
		bams = duplicates_nuclear.aligned_deduplicated,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip,
		label = "1240k_post",
		picard_jar = picard_jar,
		adna_screen_jar = adna_screen_jar,
		python_snp_target_bed = python_snp_target_bed
	}
	call damage_loop as damage_loop_nuclear{ input:
		pmdtools = pmdtools,
		python_damage_two_bases = python_damage_two_bases,
		bams = duplicates_nuclear.aligned_deduplicated,
		damage_label = "damage_nuclear",
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality
	}
	call chromosome_target as nuclear_chromosome_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = duplicates_nuclear.aligned_deduplicated,
		targets="\"{'autosome_post':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X_post':'X','Y_post':'Y','human_post':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	call chromosome_target as rsrs_chromosome_target{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = combine_mt_libraries.library_bams,
		targets="\"{'MT_pre':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	
	call duplicates as duplicates_rsrs{ input: 
		picard_jar = picard_jar,
		adna_screen_jar = adna_screen_jar,
		unsorted = combine_mt_libraries.library_bams,
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
		picard_jar = picard_jar,
		haplogrep_jar = haplogrep_jar
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
			coverages = rsrs_coverage.coverages,
			unused = signal_preliminary_report_ready.finished
		}
	}
	
	call angsd_contamination{ input:
		bams = duplicates_nuclear.aligned_deduplicated,
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		angsd = angsd,
		python_angsd_results = python_angsd_results,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip
	}
	
	call central_measures as central_measures_nuclear{ input:
		python_central_measures = python_central_measures,
		mean_label = "mean_nuclear",
		median_label = "median_nuclear",
		histograms = nuclear_chromosome_target_post.length_histogram
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

	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_duplicates_nuclear{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = duplicates_nuclear.duplicates_statistics
	}
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_duplicates_rsrs{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = duplicates_rsrs.duplicates_statistics
	}
	
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_pre_nuclear{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = nuclear_chromosome_target.target_stats
	}
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_post_nuclear{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = nuclear_chromosome_target_post.target_stats
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
		aggregate_statistics_across_sequencing_runs.statistics,
		aggregate_statistics_duplicates_nuclear.statistics,
		aggregate_statistics_duplicates_rsrs.statistics,
		aggregate_statistics_pre_nuclear.statistics,
		aggregate_statistics_post_nuclear.statistics,
		aggregate_statistics_pre_rsrs.statistics,
		aggregate_statistics_post_rsrs.statistics
	]
	call demultiplex_align_bams.aggregate_statistics as aggregate_statistics_final{ input:
		adna_screen_jar = adna_screen_jar,
		statistics_by_group = cumulative_statistics
	}
	
	call concatenate as concatenate_count_1240k_pre{ input:
		to_concatenate = count_1240k_pre.snp_target_stats
	}
	call concatenate as concatenate_count_1240k_post{ input:
		to_concatenate = count_1240k_post.snp_target_stats
	}
	call concatenate as concatenate_contammix{ input:
		to_concatenate = contammix.contamination_estimate
	}
	
	# same as final statistics, but missing the expensive contammix calculation
	Array[File] preliminary_keyed_statistics = [
		damage_loop_nuclear.damage_all_samples_two_bases,
		damage_loop_rsrs.damage_all_samples_two_bases,
		central_measures_nuclear.central_measures_output,
		central_measures_rsrs.central_measures_output,
		summarize_haplogroups.haplogroups,
		concatenate_count_1240k_pre.concatenated, 
		concatenate_count_1240k_post.concatenated,
		preseq.results,
		angsd_contamination.contamination
	]
	Array[File] final_keyed_statistics = [
		damage_loop_nuclear.damage_all_samples_two_bases,
		damage_loop_rsrs.damage_all_samples_two_bases,
		central_measures_nuclear.central_measures_output,
		central_measures_rsrs.central_measures_output,
		summarize_haplogroups.haplogroups,
		concatenate_count_1240k_pre.concatenated, 
		concatenate_count_1240k_post.concatenated,
		preseq.results,
		angsd_contamination.contamination,
		concatenate_contammix.concatenated
	]
	call prepare_report as preliminary_report{ input:
		python_prepare_report = python_prepare_report,
		aggregated_statistics = aggregate_statistics_final.statistics,
		keyed_statistics = preliminary_keyed_statistics,
		index_barcode_keys = index_barcode_keys,
		labeled_barcodes = labeled_barcodes,
		labeled_i5 = labeled_i5,
		labeled_i7 = labeled_i7,
		dataset_label = dataset_label,
		date = date
	}
	Array[File] preliminary_report_array = [preliminary_report.report, versions.versions]
	call demultiplex_align_bams.copy_output as preliminary_copy_report{ input:
		files = preliminary_report_array,
		output_path = output_path
	}
	call signal_preliminary_report_ready{ input:
		date_string = date,
		name = dataset_label,
		unused = preliminary_copy_report.copied
	}
	call prepare_report{ input:
		python_prepare_report = python_prepare_report,
		aggregated_statistics = aggregate_statistics_final.statistics,
		keyed_statistics = final_keyed_statistics,
		index_barcode_keys = index_barcode_keys,
		labeled_barcodes = labeled_barcodes,
		labeled_i5 = labeled_i5,
		labeled_i7 = labeled_i7,
		dataset_label = dataset_label,
		date = date
	}
	Array[File] report_array = [prepare_report.report, versions.versions]
	call demultiplex_align_bams.copy_output as copy_report{ input:
		files = report_array,
		output_path = output_path
	}
}

task analysis_versions{
	File adna_screen_jar
	File picard_jar
	File htsbox
	File haplogrep_jar
	File angsd
	String python_version_git_hash
	File index_barcode_keys_to_stop_call_caching

	command{
		echo "adna-workflow " >> versions
		python3 ${python_version_git_hash} >> versions
		java -version >> versions 2>&1
		python3 --version >> versions 2>&1
		bwa >> versions 2>&1
		samtools --version >> versions 2>&1
		
		echo "reichlab adna_jar " >> versions
		java -jar ${adna_screen_jar} version >> versions
		
		${htsbox} >> versions 2>&1
		echo "mafft" >> versions
		mafft --version >> versions 2>&1
		java -jar ${haplogrep_jar} >> versions 2>&1
		preseq >> versions 2>&1
		${angsd} >> versions 2>&1
	}
	output{
		File versions = "versions"
	}
	runtime{
		runtime_minutes: 5
		requested_memory_mb_per_core: 1000
		continueOnReturnCode: true
	}
}

# bams entering this task are filtered to aligned reads only and sorted by coordinate
task combine_bams_into_libraries{
	File bam_lists
	File picard_jar
	
	Int processes = 6
	
	Int num_files = length(read_lines(bam_lists))
	
	command{
		python3 <<CODE
		from multiprocessing import Pool
		from os.path import basename
		import subprocess
		from pathlib import Path
		
		def merge_bam(bam_filenames):
			sample_id_filename = basename(bam_filenames[0])
			#print('combine bams ' + sample_id_filename)
			
			# remove unaligned reads from bam components in temp directory
			temp_dir = sample_id_filename + '_dir'
			Path(temp_dir).mkdir(exist_ok=True) # make a temporary directory
			count = 0 # ensure filenames are unique by using count
			filtered_bam_filenames = []
			with open('%s/stdout_build' % (temp_dir,), 'w') as stdout_build, \
				open('%s/stderr_build' % (temp_dir,), 'w') as stderr_build:
				for bam in bam_filenames:
					count += 1
					filtered_bam_filename = '%s/%d_%s' % (temp_dir, count, basename(bam))
					subprocess.run(['samtools', 'view', '-h', '-b', '-F', '4', '-o', filtered_bam_filename, bam], check=True, stdout=stdout_build, stderr=stderr_build)
					filtered_bam_filenames += [filtered_bam_filename]
			
				merge_file_list = 'I=' + ' I='.join(filtered_bam_filenames)
				command = "java -Xmx2000m -jar ${picard_jar} MergeSamFiles %s O=%s SORT_ORDER=coordinate COMPRESSION_LEVEL=9" % (merge_file_list, sample_id_filename)
				#print('combine bam lists ' + command)
				subprocess.check_output(command, shell=True)
		
		with open("${bam_lists}") as f:
			bam_filenames_for_library = [line.split() for line in f]
			pool = Pool(processes=${processes})
			results = [pool.apply_async(merge_bam, args=(bam_filenames,)) for bam_filenames in bam_filenames_for_library]
			pool.close()
			[result.get() for result in results] # check for errors
			pool.join()
		CODE
	}
	output{
		Array[File] library_bams = glob("*.bam")
	}
	runtime{
		cpus: if processes < num_files then processes else num_files
		runtime_minutes: 180
		requested_memory_mb_per_core: 2200
	}
}

task duplicates{
	File picard_jar
	File adna_screen_jar
	Array[File] unsorted
	String duplicates_label
	
	Int processes = 8
	
	command{
		set -e
		mkdir deduplicated
		
		python3 <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def deduplicate_bam(bam):
			sample_id_filename = basename(bam)
			sample_id_filename_no_extension, extension = splitext(sample_id_filename)
			
			sorted_bam = sample_id_filename_no_extension + ".sorted_coordinate.bam"
			
			subprocess.check_output("java -Xmx4500m -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (bam, sorted_bam), shell=True)
			subprocess.check_output("java -Xmx4500m -jar ${picard_jar} MarkDuplicates I=%s O=deduplicated/%s M=%s.dedup_stats REMOVE_DUPLICATES=true BARCODE_TAG=XD ADD_PG_TAG_TO_READS=false MAX_FILE_HANDLES=1000" % (sorted_bam, sample_id_filename, sample_id_filename), shell=True)
			subprocess.check_output("java -Xmx4500m -jar ${adna_screen_jar} ReadMarkDuplicatesStatistics -l ${duplicates_label} %s.dedup_stats > %s.stats" % (sample_id_filename, sample_id_filename), shell=True)
		
		bams_string = "${sep=',' unsorted}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(deduplicate_bam, args=(bam,)) for bam in bams]
		pool.close()
		[result.get() for result in results] # check for errors
		pool.join()
		CODE
	}
	output{
		Array[File] aligned_deduplicated = glob("deduplicated/*.bam")
		Array[File] duplicates_statistics = glob("*.stats")
	}
	runtime{
		cpus: if length(unsorted) < processes then length(unsorted) else processes
		runtime_minutes: 480
		requested_memory_mb_per_core: 5000
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
	Int maximum_reads_for_damage = 5000000
	
	command{
		set -e
		python3 <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def num_reads_in_bam(bam_filename):
			result = subprocess.run(['samtools', 'view', '-c', bam_filename], check=True, universal_newlines=True, stdout=subprocess.PIPE).stdout.strip()
			return int(result)
		
		def damage_for_bam(bam):
			sample_id_filename = basename(bam)
			sample_id_filename_no_extension, extension = splitext(sample_id_filename)
			
			damage_filename = sample_id_filename_no_extension + ".damage"
			
			num_reads = num_reads_in_bam(bam)
			subsample_option = "" if num_reads <= ${maximum_reads_for_damage} else "-s %f" % (${maximum_reads_for_damage} / num_reads)
			
			subprocess.check_output("samtools view %s %s | python3 ${pmdtools} -d --requiremapq=${minimum_mapping_quality} --requirebaseq=${minimum_base_quality} > %s" %(subsample_option, bam, damage_filename), shell=True)
			damage_result = subprocess.check_output("python3 ${python_damage_two_bases} ${damage_label} %s" % (damage_filename,), shell=True)
			return damage_result.strip()
			
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(damage_for_bam, args=(bam,)) for bam in bams]
		pool.close()
		pool.join()
		with open('damage_all_samples_two_bases', 'w') as f:
			for result in results:
				f.write(result.get().decode('utf-8').strip())
				f.write('\n')
		CODE
	}
	output{
		File damage_all_samples_two_bases = "damage_all_samples_two_bases"
	}
	runtime{
		cpus: if length(bams) < processes then length(bams) else processes
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
		python3 ${python_target} ${adna_screen_jar} ${targets} ${minimum_mapping_quality} ${sep=' ' bams}
	}
	output{
		Array[File] target_stats = glob("*.stats")
		Array[File] length_histogram = glob("*.histogram")
	}
	runtime{
		cpus: 2
		runtime_minutes: 180
		requested_memory_mb_per_core: 3000
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
			python3 ${python_coverage} $bam_stat ${reference_length} $sample_id ${coverage_field} >> coverages
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
	
	Int processes = 6

	command{
		set -e
		
		python3 <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def count_SNPs(bam):
			sample_id_filename = basename(bam)
			sample_id_filename_no_extension, extension = splitext(sample_id_filename)
			
			clipped_bam = sample_id_filename_no_extension + ".clipped.bam"
			clipped_sorted_bam = sample_id_filename_no_extension + ".clipped.sorted.bam"
			
			subprocess.check_output("java -Xmx3500m -jar ${adna_screen_jar} softclip -b -n ${deamination_bases_to_clip} -i %s -o %s" % (bam, clipped_bam), shell=True)
			subprocess.check_output("java -Xmx3500m -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (clipped_bam, clipped_sorted_bam), shell=True)
			subprocess.check_output("samtools index %s" % (clipped_sorted_bam,), shell=True)
			subprocess.check_output("samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_autosome} %s > %s.autosome" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_x}        %s > %s.x" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_y}        %s > %s.y" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("python3 ${python_snp_target_bed} ${label} %s.autosome %s.x %s.y > %s.snp_target_stats" % (sample_id_filename, sample_id_filename, sample_id_filename, sample_id_filename), shell=True)
			
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(count_SNPs, args=(bam,)) for bam in bams]
		pool.close()
		[result.get() for result in results] # check for errors
		pool.join()
		CODE
	}
	output{
		Array[File] snp_target_stats = glob("*.snp_target_stats")
	}
	runtime{
		cpus: if length(bams) < processes then length(bams) else processes
		runtime_minutes: 480
		requested_memory_mb_per_core: 4000
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
		
		python3 <<CODE
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
			subprocess.check_output("samtools depth -b ${coordinates_autosome} -q ${base_quality} -Q ${mapping_quality} %s | python3 ${python_depth_histogram} > %s.autosome" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools depth -b ${coordinates_x}        -q ${base_quality} -Q ${mapping_quality} %s | python3 ${python_depth_histogram} > %s.x" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools depth -b ${coordinates_y}        -q ${base_quality} -Q ${mapping_quality} %s | python3 ${python_depth_histogram} > %s.y" % (clipped_sorted_bam, sample_id_filename), shell=True)
			subprocess.check_output("python3 ${python_depth_histogram_combine} autosome %s.autosome X %s.x Y %s.y > %s.depth_histogram" % (sample_id_filename, sample_id_filename, sample_id_filename, sample_id_filename), shell=True)
			
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(calculate_SNP_depths, args=(bam,)) for bam in bams]
		pool.close()
		[result.get() for result in results] # check for errors
		pool.join()
		CODE
	}
	output{
		Array[File] depth_histograms = glob("*.depth_histogram")
	}
	runtime{
		cpus: if length(bams) < processes then length(bams) else processes
		runtime_minutes: 480
		requested_memory_mb_per_core: 3000
	}
}

task concatenate{
	Array[File] to_concatenate
	
	command{
		set -e
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
	
	Int processes = 8
	
	command{
		set -e
		
		python3 <<CODE
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
		results = [pool.apply_async(haplogrep_run, args=(bam,)) for bam in bams]
		pool.close()
		[result.get() for result in results] # check for errors
		pool.join()
		CODE
	}
	output{
		Array[File] haplogroup_report = glob("*.haplogroup")
	}
	runtime{
		cpus: if length(bams) < processes then length(bams) else processes
		runtime_minutes: 360
		requested_memory_mb_per_core: 2000
	}
}

task summarize_haplogroups{
	File python_haplogroup
	Array[File] haplogrep_output
	
	command{
		python3 ${python_haplogroup} ${sep=' ' haplogrep_output} > haplogroups
	}
	output{
		File haplogroups = "haplogroups"
	}
	runtime{
		runtime_minutes: 5
		requested_memory_mb_per_core: 100
	}
}

task central_measures{
	File python_central_measures
	String median_label
	String mean_label
	Array[File] histograms
	
	command{
		python3 ${python_central_measures} ${median_label} ${mean_label} ${sep=' ' histograms} > central_measures
	}
	output{
		File central_measures_output = "central_measures"
	}
	runtime{
		runtime_minutes: 20
		requested_memory_mb_per_core: 2000
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
	
	Float bam_size = size(bam)
	# used to delay contammix until after other tasks
	Int? unused
	
	#samtools mpileup -u -Q ${minimum_base_quality} -q ${minimum_mapping_quality} -f ${reference} ${bam} | bcftools call -c -O z --ploidy 1 -o calls.vcf.gz
	#tabix calls.vcf.gz
	#cat ${reference} | bcftools consensus calls.vcf.gz > consensus.fa
	
	command{
		${htsbox} pileup -f ${reference} -Q ${minimum_base_quality} -q ${minimum_mapping_quality} -M ${bam} > ${sample_id}.consensus.fa
		java -jar ${picard_jar} SamToFastq I=${bam} FASTQ=for_alignment_to_consensus.fastq 
		bwa index ${sample_id}.consensus.fa
		bwa aln -t ${threads} -o ${max_open_gaps} -n ${missing_alignments_fraction} -l ${seed_length} ${sample_id}.consensus.fa for_alignment_to_consensus.fastq > realigned.sai
		bwa samse ${sample_id}.consensus.fa realigned.sai for_alignment_to_consensus.fastq | samtools view -bS - > realigned.bam
		
		python3 <<CODE
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
		python3 ${python_contammix_multiprocess} ${copies} ${sample_id} ${contammix_estimate} ${sample_id}.downsampled.bam multiple_alignment.fa ${chains} ${minimum_base_quality} ${deamination_bases_to_clip} ${seed} > contamination_estimate
	}
	output{
		File contamination_estimate = "contamination_estimate"
		File consensus = "${sample_id}.consensus.fa"
	}
	runtime{
		cpus: if bam_size < 10000 then 1 else threads
		runtime_minutes: if bam_size < 10000 then 30 else 150
		requested_memory_mb_per_core: if bam_size < 10000 then 3000 else 1000
	}
}

task preseq{
	Array[File] bams
	File targets_bed
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	File statistics
	
	File adna_screen_jar
	File picard_jar
	File python_depth_histogram
	File python_preseq_process
	File python_histogram_counts
	
	Int subsampling_iterations = 10
	
	Int processes = 8
	Float model_a
	Float model_b
	
	command{
		python3 <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		MINIMUM_RAW_READS_TO_SEQUENCE = 10e7
		
		def count_unique_reads(filename): 
			unique_count = 0
			total_count = 0
			with open(filename) as f:
				for line in f:
					depth, count = line.split()
					unique_count += int(count)
					total_count += int(depth) * int(count)
			return unique_count, total_count
		
		def preseq_run(bam):
			sample_id_filename = basename(bam)
			sample_id, extension = splitext(sample_id_filename)
			sample_id_key_not_filename = sample_id.replace('-', ':')
			
			# bed locations only
			filtered_filename = sample_id + ".filtered.bam"
			subprocess.check_output("java -Xmx4500m -jar ${adna_screen_jar} FilterSAM -b -i %s -o %s -n ${deamination_bases_to_clip} -m ${minimum_mapping_quality} -q ${minimum_base_quality} -p ${targets_bed}" % (bam, filtered_filename), shell=True)
			# sort
			sorted_filename = sample_id + ".sorted.bam"
			subprocess.check_output("java -Xmx4500m -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (filtered_filename, sorted_filename), shell=True)
			
			# demultiplexing reads
			raw_count_str = subprocess.check_output("java -Xmx4500m -jar ${adna_screen_jar} AggregateStatistics -k %s -l raw ${statistics}" % (sample_id_key_not_filename), shell=True)
			raw_count = int(raw_count_str)
			
			# build histogram
			unique_reads_histogram_filename = sample_id + ".unique_reads_histogram"
			subprocess.check_output("java -Xmx4500m -jar ${adna_screen_jar} DuplicatesHistogram -i %s > %s" % (sorted_filename, unique_reads_histogram_filename), shell=True)
			unique_read_count, total_count = count_unique_reads(unique_reads_histogram_filename)
			
			preseq_table_filename = sample_id + ".preseq_table"
			read_ratio = (raw_count / total_count) if total_count > 0 else float('inf')
			perform_subsampling = (total_count >= 1000) and (total_count < 1e5 or ((unique_read_count / total_count) < 0.20))
			perform_extrapolation = (unique_read_count > 0) and ((unique_read_count / total_count) > 0.01)
			
			if perform_subsampling:
				# For low (but not super low) complexity libraries, we want a marginal uniqueness estimate
				# We compute this using a subsampling approach to guarantee an estimate for libraries where preseq fails
				subsampling_table = sample_id + ".low.preseq_table"
				fractions = [0.01, 0.02, 0.04, 0.08, 0.16]
				if not perform_extrapolation:
					fractions += [0.32, 0.64, 0.8, 1.0]
				for fraction_retained in fractions:
					subprocess.run("python3 ${python_histogram_counts} %s -f %f -n %d >> %s" % (unique_reads_histogram_filename, fraction_retained, ${subsampling_iterations}, subsampling_table), shell=True)
			if perform_extrapolation:
				# This is extrapolation estimate for what we get from more sequencing
				# First step should be after contents of interpolation preseq table
				preseq_table_filename_extrapolation = sample_id + ".high.preseq_table"
				step = int(total_count / 4)
				extrapolation_max = int((total_count * 5) if read_ratio > 0 else 0)
				subprocess.run("preseq lc_extrap -H %s -s %d -e %d > %s" % (unique_reads_histogram_filename, step, extrapolation_max, preseq_table_filename_extrapolation), shell=True)
				
			if unique_read_count > 0:
				# combine low and high coverage preseq tables into one
				with open(preseq_table_filename, 'w') as combined_table:
					data_lines_to_write = []
					if perform_subsampling:
						with open(subsampling_table) as low_table:
							for line in low_table:
								data_lines_to_write.append(line)
					if perform_extrapolation:
						with open(preseq_table_filename_extrapolation) as high_table:
							high_table.readline() # skip header
							high_table.readline() # skip zeros
							for line in high_table:
								data_lines_to_write.append(line)
					if len(data_lines_to_write) > 0:
						header = 'TOTAL_READS\tEXPECTED_DISTINCT\tLOWER_0.95CI\tUPPER_0.95CI'
						zeros = '0\t0\t0\t0'
						print(header, file=combined_table)
						print(zeros, file=combined_table)
						for line in data_lines_to_write:
							print(line, end='', file=combined_table)
			else:
				subprocess.run("touch %s" % (preseq_table_filename), shell=True)
			
			targets_histogram_filename = sample_id + ".targets_histogram"
			subprocess.check_output("samtools depth -b ${targets_bed} -q ${minimum_base_quality} -Q ${minimum_mapping_quality} %s | python3 ${python_depth_histogram} > %s" % (sorted_filename, targets_histogram_filename), shell=True)
			# keyed statistics are written to stdout 
			result = subprocess.check_output("python3 ${python_preseq_process} %s %s -n %d -r %d -a ${model_a} -b ${model_b} -k %s | tee %s" % (targets_histogram_filename, preseq_table_filename, raw_count, total_count, sample_id_key_not_filename, sample_id + '.final_results'), shell=True)
			return result.strip()
		
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(preseq_run, args=(bam,)) for bam in bams]
		pool.close()
		pool.join()
		with open('preseq_results', 'w') as f:
			for result in results:
				f.write(result.get().decode('utf-8').strip())
				f.write('\n')
		CODE
	}
	output{
		File results = "preseq_results"
	}
	runtime{
		cpus: if length(bams) < processes then length(bams) else processes
		runtime_minutes: 400
		requested_memory_mb_per_core: 5000
	}
}

# compute contamination using X chromosome for 1240k data
task angsd_contamination{
	Array[File] bams
	
	File adna_screen_jar
	File picard_jar
	File angsd
	File angsd_contamination_bin
	File python_angsd_results
	
	File HapMap
	
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	
	Int processes = 10
	
	Int angsd_threads
	Int seed
	
	command{
		python3 <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def angsd_run(bam):
			sample_id_filename = basename(bam)
			sample_id, extension = splitext(sample_id_filename)
			sample_id_key_not_filename = sample_id.replace('-', ':')
			
			clipped_bam_filename = sample_id + ".clipped.bam"
			subprocess.run("java -Xmx3500m -jar ${adna_screen_jar} softclip -b -n ${deamination_bases_to_clip} -i %s -o %s" % (bam, clipped_bam_filename), shell=True, check=True)
			sorted_bam_filename = sample_id + ".sorted.bam"
			subprocess.run("java -Xmx3500m -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (clipped_bam_filename, sorted_bam_filename), shell=True, check=True)
			subprocess.run("samtools index %s" % (sorted_bam_filename,), shell=True, check=True)
			subprocess.run("${angsd} -i %s -r X:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ ${minimum_mapping_quality} -minQ ${minimum_base_quality} -out %s" % (sorted_bam_filename, sample_id), shell=True, check=True)
			
			angsd_output_filename = sample_id + ".angsd"
			subprocess.run("${angsd_contamination_bin} -a %s -h ${HapMap} -p ${angsd_threads} -s ${seed} > %s 2>&1" % (sample_id + ".icnts.gz", angsd_output_filename), shell=True, check=False)
			result = subprocess.run("python3 ${python_angsd_results} %s | tee %s" % (angsd_output_filename, sample_id  + ".keyed_angsd"), shell=True, check=True, stdout=subprocess.PIPE, encoding='utf-8')
			return sample_id_key_not_filename + '\t' + result.stdout.strip()
		
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(angsd_run, args=(bam,)) for bam in bams]
		pool.close()
		pool.join()
		with open('angsd_contamination_results', 'w') as f:
			for result in results:
				f.write(result.get())
				f.write('\n')
		CODE
	}
	output{
		File contamination = "angsd_contamination_results"
	}
	runtime{
		cpus: if length(bams) < processes then length(bams) else processes
		runtime_minutes: 480
		requested_memory_mb_per_core: 4000
	}
}

task signal_preliminary_report_ready{
	String django_manage_for_command
	String date_string
	String name
	
	Int unused
	
	command{
		ssh -t o2.hms.harvard.edu ssh app-wsgi-prod01.rc.hms.harvard.edu /opt/rc/python/3.8.1/bin/python3 ${django_manage_for_command} preliminary_report_ready --date_string ${date_string} --name ${name}
	}
	runtime{
		runtime_minutes: 20
		requested_memory_mb_per_core: 2000
	}
	output{
		Int finished = unused
	}
}

task prepare_report{
	File python_prepare_report
	File aggregated_statistics
	File index_barcode_keys
	# these files contain statistics that are not necessarily count based
	# they do not contain a leading number of reads
	Array[File] keyed_statistics
	
	File labeled_barcodes
	File labeled_i5
	File labeled_i7
	
	String dataset_label
	String date
	
	command{
		python3 ${python_prepare_report} --barcode_labels_filename ${labeled_barcodes} --i5_labels_filename ${labeled_i5} --i7_labels_filename ${labeled_i7} ${aggregated_statistics} ${index_barcode_keys} ${sep=' ' keyed_statistics} > ${date}_${dataset_label}.report
	}
	output{
		File report = "${date}_${dataset_label}.report"
	}
	runtime{
		runtime_minutes: 20
		requested_memory_mb_per_core: 2000
	}
}
