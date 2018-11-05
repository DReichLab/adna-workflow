import "demultiplex.wdl" as demultiplex_align_bams
import "analysis.wdl" as analysis
import "release_and_pulldown.wdl" as pulldown

workflow sample_merge_and_pulldown_with_analysis{
	File sample_library_list
	String release_directory
	String genome_reference_string
	String mt_reference_string

	File adna_screen_jar
	File picard_jar
	File htsbox
	File pmdtools
	File haplogrep_jar
	
	File mt_reference_rsrs_in
	File mt_reference_rcrs_in

	File python_damage_two_bases
	File python_angsd_results
	File python_target
	File python_coverage

	Float missing_alignments_fraction
	Int max_open_gaps
	Int seed_length
	Int minimum_mapping_quality
	Int minimum_base_quality
	Int deamination_bases_to_clip
	
	call demultiplex_align_bams.prepare_reference as prepare_reference_rsrs{ input:
		reference = mt_reference_rsrs_in
	}
	call demultiplex_align_bams.prepare_reference as prepare_reference_rcrs{ input:
		reference = mt_reference_rcrs_in
	}

	call prepare_bam_list{ input:
		library_id_file = sample_library_list,
		genome_reference_string = genome_reference_string,
		mt_reference_string = mt_reference_string
	}
	call merge_bams as merge_bams_nuclear{ input : 
		bam_lists_per_individual = prepare_bam_list.nuclear_list,
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		reference = genome_reference_string
	}
	call release_samples as release_samples_nuclear { input:
		release_directory = release_directory,
		bams = merge_bams_nuclear.bams,
		sample_library_list = sample_library_list
	}
	call merge_bams as merge_bams_mt{ input:
		bam_lists_per_individual = prepare_bam_list.mt_list,
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		reference = genome_reference_string
	}
	call release_samples as release_samples_mt{ input:
		release_directory = release_directory,
		bams = merge_bams_mt.bams,
		sample_library_list = sample_library_list
	}
	call analysis.damage_loop as damage_nuclear{ input :
		pmdtools = pmdtools,
		python_damage_two_bases = python_damage_two_bases,
		bams = merge_bams_nuclear.bams,
		damage_label = "damage_nuclear",
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		processes = 1
	}
	call analysis.damage_loop as damage_mt{ input :
		pmdtools = pmdtools,
		python_damage_two_bases = python_damage_two_bases,
		bams = merge_bams_mt.bams,
		damage_label = "damage_mt",
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		processes = 1
	}
	call analysis.angsd_contamination{ input:
		bams = merge_bams_nuclear.bams,
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		python_angsd_results = python_angsd_results,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip,
		processes = 1
	}
	call analysis.haplogrep as haplogrep_rcrs{ input:
		missing_alignments_fraction = missing_alignments_fraction,
		max_open_gaps = max_open_gaps,
		seed_length = seed_length,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip,
		bams = merge_bams_mt.bams,
		reference = prepare_reference_rcrs.reference_fa,
		reference_amb = prepare_reference_rcrs.reference_amb,
		reference_ann = prepare_reference_rcrs.reference_ann,
		reference_bwt = prepare_reference_rcrs.reference_bwt,
		reference_pac = prepare_reference_rcrs.reference_pac,
		reference_sa = prepare_reference_rcrs.reference_sa,
		adna_screen_jar = adna_screen_jar,
		picard_jar = picard_jar,
		haplogrep_jar = haplogrep_jar,
		processes = 1
	}
	call analysis.summarize_haplogroups{ input:
		haplogrep_output = haplogrep_rcrs.haplogroup_report
	}
	call analysis.chromosome_target as rsrs_chromosome_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = merge_bams_mt.bams,
		targets="\"{'MT_post':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}
	call analysis.chromosome_coverage as rsrs_coverage{ input:
		bam_stats = rsrs_chromosome_target_post.target_stats,
		python_coverage = python_coverage,
		reference_length = 16569,
		coverage_field = "MT_post-coverageLength"
	}
	scatter(bam in merge_bams_mt.bams){
		call analysis.contammix{ input:
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
	call analysis.concatenate as concatenate_contammix{ input:
		to_concatenate = contammix.contamination_estimate
	}
	call analysis.snp_target_bed as count_1240k_post { input:
#		coordinates_autosome = coordinates_1240k_autosome,
#		coordinates_x = coordinates_1240k_x,
#		coordinates_y = coordinates_1240k_y,
		bams = merge_bams_nuclear.bams,
		minimum_mapping_quality = minimum_mapping_quality,
		minimum_base_quality = minimum_base_quality,
		deamination_bases_to_clip = deamination_bases_to_clip,
		label = "1240k_post",
		picard_jar = picard_jar,
		adna_screen_jar = adna_screen_jar,
#		python_snp_target_bed = python_snp_target_bed,
		processes = 1
	}
	call analysis.concatenate as concatenate_count_1240k_post{ input:
		to_concatenate = count_1240k_post.snp_target_stats
	}
	call pulldown_merged_samples{ input:
		bams = release_samples_nuclear.released_bams,
		sex_by_instance_id = concatenate_count_1240k_post.concatenated,
		sample_bam_list = sample_library_list
	}
	call analysis_results{ input:
		keyed_value_results = [
			damage_nuclear.damage_all_samples_two_bases,
			damage_mt.damage_all_samples_two_bases,
			angsd_contamination.contamination,
			concatenate_count_1240k_post.concatenated,
			rsrs_coverage.coverages,
			concatenate_contammix.concatenated,
			summarize_haplogroups.haplogroups
		]
	}
}

# take a list of instance ids with library ids, build corresponding list of bam paths for nuclear and mt components
task prepare_bam_list{
	# each line has an individual/instance ID, UDG treatment, then component library ids
	# does each line also need version number?
	File library_id_file
	File python_bam_finder
	File python_library_id
	String genome_reference_string
	String mt_reference_string
	
	command{
		python3 <<CODE
		import subprocess
		
		def bam_path_from_library_id(library_id, reference):
			result = subprocess.run(['python3', '${python_bam_finder}', '--reference', reference, library_id], check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
			return result
		
		with open("${library_id_file}", 'r') as f, open('nuclear_list', 'w') as nuclear_filelist, open('mt_list', 'w') as mt_filelist:
			shop_bam_root = 'aln.sort.mapped.rmdupse_adna_v2.md'
			for line in f:
				fields = line.split('\t')
				individual_id = fields[0]
				instance_id = fields[1]
				udg = fields[2] # not needed for merge
				library_ids = fields[3:]
				
				nuclear_fields = [instance_id]
				mt_fields = [instance_id]
				
				for library_id in library_ids:
					# whole genome
					nuclear = bam_path_from_library_id(library_id, "${genome_reference_string}")
					mt = bam_path_from_library_id(library_id, "${mt_reference_string}")
					
					# all libraries should have a nuclear and MT component
					if nuclear == '' or mt == '':
						raise ValueError('missing bam for %s' % library_id)
					else:
						nuclear_fields.append(library_id)
						nuclear_fields.append(nuclear)
						mt_fields.append(library_id)
						mt_fields.append(mt)
				print('\t'.join(nuclear_fields), file=nuclear_filelist)
				print('\t'.join(mt_fields), file=mt_filelist)
		CODE
	}
	output{
		File nuclear_list = "nuclear_list"
		File mt_list = "mt_list"
	}
	runtime{
		runtime_minutes: 10
		requested_memory_mb_per_core: 100
	}
}

# perform a "by-sample" (individual/instance) merge
task merge_bams{
	# each line has an instance ID, then component bam paths
	File bam_lists_per_individual
	File adna_screen_jar
	File picard_jar
	String reference
	
	Int processes = 1
	
	command{
		python3 <<CODE
		from multiprocessing import Pool
		from os.path import basename
		import subprocess
		import os
		
		def merge_bam(instance_id, library_ids, bam_paths):
			instance_id_filename = "%s.${reference}.bam" % (instance_id)
			# make a directory for this instance ID
			os.makedirs(instance_id)
			# write instance ID into read groups
			bams_with_altered_read_groups = []
			for library_id, bam in zip(library_ids, bam_paths):
				bam_with_altered_read_groups = instance_id + '/' + basename(bam)
				subprocess.run(["java", "-Xmx2500m", "-jar", "${adna_screen_jar}", "ReadGroupRewrite", "-i", bam, "-o", bam_with_altered_read_groups, "-s", instance_id, "-l", library_id], check=True)
				bams_with_altered_read_groups.append(bam_with_altered_read_groups)
			# merge
			merge_file_list = 'I=' + ' I='.join(bams_with_altered_read_groups)
			command = "java -Xmx2500m -jar ${picard_jar} MergeSamFiles %s O=%s SORT_ORDER=coordinate" % (merge_file_list, instance_id_filename)
			#print('combine bam lists ' + command)
			subprocess.check_output(command, shell=True)
		
		pool = Pool(processes=${processes})
		with open("${bam_lists_per_individual}") as f:
			for line in f:
				fields = line.split()
				instance_id = fields[0]
				library_ids = fields[1::2]
				bam_paths = fields[2::2]
				pool.apply_async(merge_bam, args=(instance_id, library_ids, bam_paths))
		pool.close()
		pool.join()
		CODE
	}

	output{
		Array[File] bams = glob("*.bam")
	}
	runtime{
		cpus: processes
		runtime_minutes: 200
		requested_memory_mb_per_core: 3000
	}
}

task analysis_results{
	File python_combine_dictionary_results
	Array[File] keyed_value_results
	
	command{
		python3 ${python_combine_dictionary_results} ${sep=' ' keyed_value_results} > results
	}
	output{
		File results = "results"
	}
	runtime{
		runtime_minutes: 5
		requested_memory_mb_per_core: 100
	}
}

task release_samples{
	String release_directory
	Array[File] bams
	File sample_library_list
	
	command{
		python3 <<CODE
		import os
		import shutil
		from pathlib import Path
		
		instance_to_individual = dict()
		with open("${sample_library_list}") as f:
			for line in f:
				fields = line.split('\t')
				individual_id = fields[0]
				instance_id = fields[1]
				instance_to_individual[instance_id] = individual_id
		
		bams_string = "${sep=',' bams}"
		bams = bams_string.split(',')
		
		with open('bam_list', 'w') as bam_list:
			for bam in bams:
				source_file = Path(bam)
				# create a directory for the individual if it does not exist yet
				instance_id = os.path.splitext(os.path.basename(bam))[0]
				individual_id = instance_to_individual[instance_id]
				bam_directory = Path(release_directory) / individual_id
				bam_directory.mkdir(exist_ok=True)
				# copy file
				bam_destination = bam_directory / instance_id
				if bam_destination.exists():
					sys.stderr.write('%s already exists' % (source_file))
				else:
					create = shutil.copy(source_file, bam_destination)
					os.chmod(created, 0o440)
				print(str(bam_destination), file=bam_list)
		CODE
	}
	output{
		Array[String] released_bams = read_lines('bam_list')
	}
	runtime{
		runtime_minutes: 60
		requested_memory_mb_per_core: 100
	}
}

task pulldown_merged_samples{
	File python_pulldown_sample
	File pulldown_executable
	String label
	String release_directory
	
	File sex_by_instance_id
	Array[File] bams
	
	File sample_bam_list
	
	command{
		python3 ${python_pulldown_sample} --pulldown_executable ${pulldown_executable} --pulldown_label ${label} --release_directory ${release_directory} ${sample_bam_list} ${sex_by_instance_id} ${sep=' ' bams}
	}
}
