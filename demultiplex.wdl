workflow demultiplex_align_bams{
	String blc_input_directory
	String dataset_label
	String date
	
	File i5_indices
	File i7_indices
	File barcodeSets
	
	File barcodes_q_only
	
	File adna_screen_jar
	File picard_jar
	
	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	
	Int samples_to_demultiplex
	
	File index_barcode_keys
	
	File python_lane_name
	File python_kmer_analysis
	
	# the references need to appear in the same directory as the derived files
	# in the prepare_reference, we put all of these into the same directory
	# all subsequent uses of the reference need to use that copy
	File reference_in
	File mt_reference_rsrs_in
	
	String output_path_parent
	String output_path = output_path_parent + "/" + date + "_" + dataset_label
	String output_path_nuclear_aligned_filtered = output_path + "/nuclear_aligned_filtered"
	String output_path_rsrs_aligned_filtered = output_path + "/rsrs_aligned_filtered"
	
	call prepare_reference as prepare_reference_nuclear{ input:
		reference = reference_in
	}
	call prepare_reference as prepare_reference_rsrs{ input:
		reference = mt_reference_rsrs_in
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
	call collect_read_group_info{ input:
		read_groups_by_lane = merge_and_trim_lane.read_group
	}
	call aggregate_statistics as aggregate_lane_statistics{ input :
		adna_screen_jar=adna_screen_jar,
		statistics_by_group=merge_and_trim_lane.statistics
	}
	call kmer_analysis{ input :
		python_kmer_analysis = python_kmer_analysis,
		barcodes_q_only = barcodes_q_only,
		counts_by_index_barcode_key = aggregate_lane_statistics.statistics,
		index_barcode_keys = index_barcode_keys,
		dataset_label = dataset_label,
		date = date
	}
	call collect_filenames{ input:
		filename_arrays = merge_and_trim_lane.fastq_to_align
	}
	String read_group = dataset_label
	scatter(fastq_to_align in collect_filenames.filenames){
		call align as align_nuclear{ input:
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
	call demultiplex as demultiplex_nuclear {input:
		adna_screen_jar = adna_screen_jar,
		prealignment_statistics = aggregate_lane_statistics.statistics,
		aligned_bam_files = align_nuclear.bam,
		samples_to_demultiplex = samples_to_demultiplex,
		index_barcode_keys = index_barcode_keys
	}
	call filter_aligned_only as filter_aligned_only_nuclear { input:
		picard_jar = picard_jar,
		bams = demultiplex_nuclear.demultiplexed_bam,
	}
	call demultiplex as demultiplex_rsrs {input:
		adna_screen_jar = adna_screen_jar,
		prealignment_statistics = aggregate_lane_statistics.statistics,
		aligned_bam_files = align_rsrs.bam,
		samples_to_demultiplex = samples_to_demultiplex,
		index_barcode_keys = index_barcode_keys
	}
	call filter_aligned_only as filter_aligned_only_rsrs{ input:
		bams = demultiplex_rsrs.demultiplexed_bam
	}
	
	# output
	call copy_output as copy_nuclear_aligned_filtered{ input:
		files = filter_aligned_only_nuclear.filtered,
		output_path = output_path_nuclear_aligned_filtered
	}
	call copy_output as copy_rsrs_aligned_filtered{ input:
		files = filter_aligned_only_rsrs.filtered,
		output_path = output_path_rsrs_aligned_filtered
	}
	Array[File] demultiplex_statistics = [demultiplex_nuclear.statistics, demultiplex_rsrs.statistics]
	call aggregate_statistics as aggregate_demultiplex_statistics{ input :
		adna_screen_jar=adna_screen_jar,
		statistics_by_group=demultiplex_statistics
	}
	Array[File] misc_output_files = [collect_read_group_info.read_groups, aggregate_demultiplex_statistics.statistics, kmer_analysis.analysis]
	call copy_output as copy_misc_output_files{input :
		files = misc_output_files,
		output_path = output_path
	}
	call copy_and_rename as copy_and_rename_lane_statistics{ input:
		source_file = aggregate_lane_statistics.statistics,
		output_path = output_path,
		output_filename_no_path = "lane_statistics"
	}
	
	output{
		Array[File] nuclear_bams = filter_aligned_only_nuclear.filtered
		Array[File] rsrs_bams = filter_aligned_only_rsrs.filtered
		File kmer_analysis_report = kmer_analysis.analysis
		File aggregated_statistics = aggregate_lane_statistics.statistics
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
		runtime_minutes: 480
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
		requested_memory_mb_per_core: 500
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
		java -Xmx3500m -jar ${adna_screen_jar} BarcodeCount --i5-indices ${i5_indices} --i7-indices ${i7_indices} --barcodes ${barcodeSets} ${read_files_by_lane[0]} ${read_files_by_lane[1]} ${read_files_by_lane[2]} ${read_files_by_lane[3]} > barcodeCount.stats
	}
	output{
		File barcode_count_statistics = "barcodeCount.stats"
	}
	runtime{
		runtime_minutes: 300
		requested_memory_mb_per_core: 4000
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
		java -Xmx7500m -jar ${adna_screen_jar} IndexAndBarcodeScreener --i5-indices ${i5_indices} --i7-indices ${i7_indices} --barcodes ${barcodeSets} --barcode-count ${barcode_count_statistics} ${read_files_by_lane[0]} ${read_files_by_lane[1]} ${read_files_by_lane[2]} ${read_files_by_lane[3]} ${label} > ${label}.stats
	}
	
	output{
		Array[File] fastq_to_align = glob("${label}*.fastq.gz")
		File statistics = "${label}.stats"
		File read_group = "read_group"
	}
	runtime{
		runtime_minutes: 300
		requested_memory_mb_per_core: 8000
	}
}

task collect_read_group_info{
	Array[File] read_groups_by_lane
	
	command{
		for file in ${sep=' ' read_groups_by_lane}  ; do 
			cat $file >> read_groups
		done
	}
	output{
		File read_groups = "read_groups"
	}
	runtime{
		runtime_minutes: 5
		requested_memory_mb_per_core: 100
	}
}

task aggregate_statistics{
	File adna_screen_jar
	Array [File] statistics_by_group
	
	command{
		java -Xmx500m -jar ${adna_screen_jar} AggregateStatistics ${sep=' ' statistics_by_group} > aggregated_statistics
	}
	output{
		File statistics = "aggregated_statistics"
	}
	runtime{
		runtime_minutes: 20
		requested_memory_mb_per_core: 1000
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
		runtime_minutes: 480
		requested_memory_mb_per_core: 6000
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
		requested_memory_mb_per_core: 1000
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
		runtime_minutes: 480
		requested_memory_mb_per_core: 8000
	}
}

# filter out unaligned reads, then sort
task filter_aligned_only{
	File picard_jar
	Array[File] bams
	Int processes = 5
	
	# picard complains "MAPQ should be 0 for unmapped read." while trying to filter unmapped reads
	#java -jar ${picard_jar} SortSam I=${bam} O=sorted_queryname.bam SORT_ORDER=queryname
	#java -jar ${picard_jar} FilterSamReads I=sorted_queryname.bam O=${filename} FILTER=includeAligned
	command{
		set -e
		mkdir filtered_sorted
		python <<CODE
		from multiprocessing import Pool
		from os.path import basename
		import subprocess
		
		def filter_bam_aligned_only(bam):
			output_filename = basename(bam)
			subprocess.check_output("samtools view -h -b -F 4 -o %s %s" % (output_filename, bam), shell=True)
			subprocess.check_output("java -Xmx3500m -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (output_filename, "filtered_sorted/" + output_filename), shell=True)
		
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		pool = Pool(processes=${processes})
		[pool.apply_async(filter_bam_aligned_only, args=(bam,)) for bam in bams]
		pool.close()
		pool.join()
		CODE
	}
	output{
		Array[File] filtered = glob("filtered_sorted/*.bam")
	}
	runtime{
		cpus: processes
		runtime_minutes: 480
		requested_memory_mb_per_core: 4000
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

task copy_and_rename{
	File source_file
	String output_path
	String output_filename_no_path
	
	command{
		mkdir -p ${output_path};
		cp -l $file "${output_path}/${output_filename_no_path}" || cp $file "${output_path}/${output_filename_no_path}"
	}
	runtime{
		runtime_minutes: 240
		requested_memory_mb_per_core: 1000
	}
}

task kmer_analysis{
	File python_kmer_analysis
	File barcodes_q_only
	File counts_by_index_barcode_key
	File index_barcode_keys
	
	String dataset_label
	String date

	command{
		python ${python_kmer_analysis} ${barcodes_q_only} ${counts_by_index_barcode_key} ${index_barcode_keys} > ${date}_${dataset_label}.kmer
	}
	output{
		File analysis = "${date}_${dataset_label}.kmer"
	}
	runtime{
		runtime_minutes: 60
		requested_memory_mb_per_core: 1000
	}
}
