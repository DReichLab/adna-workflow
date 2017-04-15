workflow ancientDNA_screen{
	String blc_input_directory
	File i5_indices
	File i7_indices
	File barcodeSets
	String alignment_reference
	
	File adna_screen_jar
	File picard_jar
	
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_fai
	File reference_pac
	File reference_rbwt
	File reference_rpac
	File reference_sa
	File reference_rsa

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
	call collect_filenames{ input:
		filename_arrays = merge_and_trim_lane.fastq_to_align
	}
	scatter(fastq_to_align in collect_filenames.filenames){
		call align { input:
			fastq_to_align = fastq_to_align,
			reference = reference,
			reference_amb = reference_amb,
			reference_ann = reference_ann,
			reference_bwt = reference_bwt,
			reference_fai = reference_fai,
			reference_pac = reference_pac,
			reference_rbwt = reference_rbwt,
			reference_rpac = reference_rpac,
			reference_rsa = reference_rsa,
			reference_sa = reference_sa
		}
		call convert_to_sam { input: 
			sai_to_sam=align.sai, 
			fastq_to_align=fastq_to_align,
			reference = reference,
			reference_amb = reference_amb,
			reference_ann = reference_ann,
			reference_bwt = reference_bwt,
			reference_fai = reference_fai,
			reference_pac = reference_pac,
			reference_rbwt = reference_rbwt,
			reference_rpac = reference_rpac,
			reference_rsa = reference_rsa,
			reference_sa = reference_sa
		}
	}
	call aggregate_lane_statistics{ input :
		adna_screen_jar=adna_screen_jar,
		statistics_by_lane=merge_and_trim_lane.statistics
	}
	call demultiplex {input:
		adna_screen_jar = adna_screen_jar,
		summary_statistics = aggregate_lane_statistics.statistics,
		aligned_sam_files = convert_to_sam.sam
	}
	scatter(sam in demultiplex.demultiplexed_sam){
		call sort{ input: 
			picard_jar = picard_jar,
			unsorted = sam
		}
		call deduplicate{ input:
			picard_jar = picard_jar,
			sorted = sort.sorted
		}
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
	}
}

task discover_lane_name_from_filename{
	String filename
	File python_lane_name
	
	command{
		python3 ${python_lane_name} ${filename} > lane_name
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
		java -jar ${adna_screen_jar} ${i5_indices} ${i7_indices} ${barcodeSets} ${read_files_by_lane[0]} ${read_files_by_lane[1]} ${read_files_by_lane[2]} ${read_files_by_lane[3]} ${label} > ${label}.stats
	}
	
	output{
		Array[File] fastq_to_align = glob("${label}*.fastq.gz")
		File statistics = "${label}.stats"
	}
	runtime{
			cpus: 1
			runtime_minutes: 300
			requested_memory_mb_per_core: 16384
			queue: "short"
	}
}

task aggregate_lane_statistics{
	File adna_screen_jar
	Array [File] statistics_by_lane
	
	command{
		java -cp ${adna_screen_jar} adnascreen.AggregateStatistics ${sep=' ' statistics_by_lane} > aggregated_statistics
	}
	output{
		File statistics = "aggregated_statistics"
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
	File reference_fai
	File reference_pac
	File reference_rbwt
	File reference_rpac
	File reference_rsa
	File reference_sa
	
	command {
		bwa aln -t ${threads} -o ${max_open_gaps} -n ${missing_alignments_fraction} -l ${seed_length} ${reference} ${fastq_to_align} > ${fastq_to_align}.sai
	}
	output{
		File sai = "${fastq_to_align}.sai"
	}
	runtime{
			cpus: "${threads}"
			runtime_minutes: 300
			requested_memory_mb_per_core: 8192
			queue: "short"
	}
}

task convert_to_sam{
	File fastq_to_align
	File sai_to_sam
	
	File reference
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_fai
	File reference_pac
	File reference_rbwt
	File reference_rpac
	File reference_rsa
	File reference_sa
	
	command{
		bwa samse ${reference} ${sai_to_sam} ${fastq_to_align} > ${sai_to_sam}.sam
	}
	output{
		File sam = "${sai_to_sam}.sam"
	}
}

# use String instead of filename to avoid file copying overhead
task collect_filenames{
	Array[Array[String]] filename_arrays
	File python_flatten
	
	command{
		echo "${sep='\n' filename_arrays}" > raw_array
		python3 ${python_flatten} < raw_array > file_of_filenames
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
	File summary_statistics
	Array[File] aligned_sam_files
	
	command{
		java -cp ${adna_screen_jar} adnascreen.DemultiplexSAM ${summary_statistics} ${sep=' ' aligned_sam_files}
	}
	output{
		Array[File] demultiplexed_sam = glob("*.sam")
	}
}

task sort{
	File picard_jar
	File unsorted
	
	command{
		java -jar ${picard_jar} SortSam I=${unsorted} O=${unsorted}.sorted SORT_ORDER=coordinate
	}
	output{
		File sorted = "${unsorted}.sorted"
	}
}

task deduplicate{
	File picard_jar
	File sorted

	command{
		java -jar ${picard_jar} MarkDuplicates I=${sorted} O=${sorted}.dedup M=${sorted}.dedup_stats
	}
	output{
		File deduplicated = "${sorted}.dedup"
		File statistics = "${sorted}.dedup_stats"
	}
}
