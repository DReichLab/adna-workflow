# tasks reworked from analysis.wdl to remove clipping
# clipping should be done beforehand, especially because it becomes complicated with mixed UDG cases

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
		cpus: 1
		runtime_minutes: 120
		requested_memory_mb_per_core: 3000
	}
}

task clip_deamination{
	File adna_screen_jar
	File picard_jar
	Array[File] bams
	Int deamination_bases_to_clip_half
	Int deamination_bases_to_clip_minus
	Int deamination_bases_to_clip_plus
	Array[String] udg_minus_libraries
	Array[String] udg_plus_libraries
	File python_read_groups_from_bam
	
	Int processes = 8
	
	command{
		set -e
		mkdir clipped_sorted
		
		python3 <<CODE
		from multiprocessing import Pool
		from os.path import basename, splitext
		import subprocess
		
		def libraries_in_bam(bam):
			output = subprocess.run(['python3', "${python_read_groups_from_bam}", '-l', bam], stdout=subprocess.PIPE, check=True)
			libraries = output.stdout.decode('utf-8').split()
			return libraries
		
		def clip_ends(bam, half_bases, minus_libraries, minus_bases, plus_libraries, plus_bases):
			sample_id_filename = basename(bam)
			sample_id_filename_no_extension, extension = splitext(sample_id_filename)
			
			clipped_bam = sample_id_filename_no_extension + ".clipped.bam"
			clipped_sorted_bam = "clipped_sorted/" + sample_id_filename_no_extension + ".bam"
			
			present_libraries = set(libraries_in_bam(bam))
			
			minus_args = []
			for minus_library in minus_libraries:
				if minus_library in present_libraries:
					minus_args.append('-s')
					minus_args.append(minus_library)
			plus_args = []
			for plus_library in plus_libraries:
				if plus_library in present_libraries:
					plus_args.append('-t')
					plus_args.append(plus_library)
			
			subprocess.run(["java", "-Xmx3500m", "-jar", "${adna_screen_jar}", "softclip", "-b", "-n", "%d" % (half_bases,), "-i", bam, "-o", clipped_bam, "-x", "%d" % (minus_bases,), "-y", "%d" % (plus_bases,)] + minus_args + plus_args, check=True)
			subprocess.check_output("java -Xmx3500m -jar ${picard_jar} SortSam I=%s O=%s SORT_ORDER=coordinate" % (clipped_bam, clipped_sorted_bam), shell=True)
			return clipped_sorted_bam
			
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		udg_minus_string = "${sep=',' udg_minus_libraries}"
		udg_minus_libs = udg_minus_string.split(',')
		udg_plus_string = "${sep=',' udg_plus_libraries}"
		udg_plus_libs = udg_plus_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(clip_ends, args=(bam, ${deamination_bases_to_clip_half}, udg_minus_libs, ${deamination_bases_to_clip_minus}, udg_plus_libs, ${deamination_bases_to_clip_plus} )) for bam in bams]
		pool.close()
		pool.join()
		[result.get() for result in results]
		CODE
	}
	runtime{
		cpus: processes
		runtime_minutes: 600
		requested_memory_mb_per_core: 3600
	}
	output{
		Array[File] clipped_bams = glob("clipped_sorted/*.bam")
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
	String label
	
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
			
			subprocess.check_output("samtools index %s" % (bam,), shell=True)
			subprocess.check_output("samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_autosome} %s > %s.autosome" % (bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_x}        %s > %s.x" % (bam, sample_id_filename), shell=True)
			subprocess.check_output("samtools view -c -q ${minimum_mapping_quality} -L ${coordinates_y}        %s > %s.y" % (bam, sample_id_filename), shell=True)
			subprocess.check_output("python3 ${python_snp_target_bed} ${label} %s.autosome %s.x %s.y > %s.snp_target_stats" % (sample_id_filename, sample_id_filename, sample_id_filename, sample_id_filename), shell=True)
			
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
		runtime_minutes: 480
		requested_memory_mb_per_core: 4000
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
	String region = "MT"
	Array[File] bams
	File haplogrep_jar
	Int phylotree_version
	File adna_screen_jar
	File picard_jar
	
	Int deamination_bases_to_clip_half
	Int deamination_bases_to_clip_minus
	Int deamination_bases_to_clip_plus
	Array[String] udg_minus_libraries
	Array[String] udg_plus_libraries
	File python_read_groups_from_bam
	
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
		
		def libraries_in_bam(bam):
			output = subprocess.run(['python3', "${python_read_groups_from_bam}", '-l', bam], stdout=subprocess.PIPE, check=True)
			libraries = output.stdout.decode('utf-8').split()
			return libraries
		
		def haplogrep_run(bam, half_bases, minus_libraries, minus_bases, plus_libraries, plus_bases):
			sample_id_filename = basename(bam)
			sample_id, extension = splitext(sample_id_filename)
			
			present_libraries = set(libraries_in_bam(bam))
			
			minus_args = []
			for minus_library in minus_libraries:
				if minus_library in present_libraries:
					minus_args.append('-s')
					minus_args.append(minus_library)
			plus_args = []
			for plus_library in plus_libraries:
				if plus_library in present_libraries:
					plus_args.append('-t')
					plus_args.append(plus_library)
			
			subprocess.check_output("java -Xmx2600m -jar ${picard_jar} SamToFastq I=%s FASTQ=%s.fastq" % (bam, sample_id), shell=True)
			subprocess.check_output("bwa aln -t 2 -o ${max_open_gaps} -n ${missing_alignments_fraction} -l ${seed_length} ${reference} %s.fastq > %s.sai" % (sample_id, sample_id), shell=True)
			subprocess.check_output("bwa samse ${reference} %s.sai %s.fastq | samtools view -bS - > %s.realigned.bam" % (sample_id, sample_id, sample_id), shell=True)
			clipped_bam = "%s.clipped_unsorted_realigned.bam" % (sample_id,)
			subprocess.run(["java", "-Xmx2600m", "-jar", "${adna_screen_jar}", "softclip", "-b", "-n", "%d" % (half_bases,), "-i", bam, "-o", clipped_bam, "-x", "%d" % (minus_bases,), "-y", "%d" % (plus_bases,)] + minus_args + plus_args, check=True)
			subprocess.check_output("java -Xmx2600m -jar ${picard_jar} SortSam I=%s.clipped_unsorted_realigned.bam O=%s.bam SORT_ORDER=coordinate" % (sample_id, sample_id), shell=True)
			subprocess.check_output("samtools index %s.bam" % (sample_id,), shell=True)
			subprocess.check_output("samtools mpileup -q ${minimum_mapping_quality} -Q ${minimum_base_quality} -r ${region} -u -f ${reference} %s.bam | bcftools call -c -v --ploidy 1 > %s.vcf" % (sample_id, sample_id), shell=True)
			subprocess.check_output("java -Xmx2600m -jar ${haplogrep_jar} --format vcf --phylotree ${phylotree_version} --in %s.vcf --out %s.haplogroup" % (sample_id, sample_id), shell=True)
		
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		udg_minus_string = "${sep=',' udg_minus_libraries}"
		udg_minus_libs = udg_minus_string.split(',')
		udg_plus_string = "${sep=',' udg_plus_libraries}"
		udg_plus_libs = udg_plus_string.split(',')
		
		pool = Pool(processes=${processes})
		results = [pool.apply_async(haplogrep_run, args=(bam, ${deamination_bases_to_clip_half}, udg_minus_libs, ${deamination_bases_to_clip_minus}, udg_plus_libs, ${deamination_bases_to_clip_plus})) for bam in bams]
		pool.close()
		pool.join()
		[result.get() for result in results]
		CODE
	}
	output{
		Array[File] haplogroup_report = glob("*.haplogroup")
	}
	runtime{
		cpus: processes
		runtime_minutes: 360
		requested_memory_mb_per_core: 3000
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
	Int deamination_bases_to_clip = 0
	
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
	File statistics
	
	File adna_screen_jar
	File python_depth_histogram
	File python_preseq_process
	
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
			subprocess.check_output("java -Xmx4500m -jar ${adna_screen_jar} FilterSAM -b -i %s -o %s -m ${minimum_mapping_quality} -q ${minimum_base_quality} -p ${targets_bed}" % (bam, filtered_filename), shell=True)
			
			# demultiplexing reads
			raw_count_str = subprocess.check_output("java -Xmx4500m -jar ${adna_screen_jar} AggregateStatistics -k %s -l raw ${statistics}" % (sample_id_key_not_filename), shell=True)
			raw_count = int(raw_count_str)
			
			# build histogram
			unique_reads_histogram_filename = sample_id + ".unique_reads_histogram"
			subprocess.check_output("java -Xmx4500m -jar ${adna_screen_jar} DuplicatesHistogram -i %s > %s" % (filtered_filename, unique_reads_histogram_filename), shell=True)
			unique_read_count, total_count = count_unique_reads(unique_reads_histogram_filename)
			
			read_ratio = (raw_count / total_count) if total_count > 0 else float('inf')
			step = int(total_count / 4)
			extrapolation_max = int((total_count * 5) if read_ratio > 0 else 0)
			preseq_table_filename = sample_id + ".preseq_table"
			if (unique_read_count > 0) and ((total_count  / unique_read_count) < 100):
				subprocess.run("preseq lc_extrap -H %s -s %d -e %d > %s" % (unique_reads_histogram_filename, step, extrapolation_max, preseq_table_filename), shell=True)
			else: # avoid running preseq for low complexity samples
				subprocess.run("touch %s" % (preseq_table_filename), shell=True)
			
			
			targets_histogram_filename = sample_id + ".targets_histogram"
			subprocess.check_output("samtools depth -b ${targets_bed} -q ${minimum_base_quality} -Q ${minimum_mapping_quality} %s | python3 ${python_depth_histogram} > %s" % (filtered_filename, targets_histogram_filename), shell=True)
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
		cpus: processes
		runtime_minutes: 200
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
			
			subprocess.run("samtools index %s" % (bam,), shell=True, check=True)
			subprocess.run("${angsd} -i %s -r X:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ ${minimum_mapping_quality} -minQ ${minimum_base_quality} -out %s" % (bam, sample_id), shell=True, check=True)
			
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
		cpus: processes
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
		ssh -t mym11@orchestra-legacy.med.harvard.edu ssh rc-app-shared01.orchestra /opt/python-3.4.2/bin/python ${django_manage_for_command} preliminary_report_ready --date_string ${date_string} --name ${name}
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
