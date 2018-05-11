workflow release_and_pulldown{
	File adna_jar
	File picard_jar
	File python_release_libraries
	File python_pulldown
	File python_read_groups_from_bam
	File pulldown_binary
	
	File report
	File nuclear_bamlist
	File mt_bamlist
	
	String release_parent_directory
	String label

	call build_release_libraries as build_mt{input:
		bamlist = mt_bamlist,
		adna_jar = adna_jar,
		picard_jar = picard_jar,
		python_release_libraries = python_release_libraries,
		python_read_groups_from_bam = python_read_groups_from_bam,
		release_parent_directory = release_parent_directory
	}
	call build_release_libraries as build_nuclear{input:
		bamlist = nuclear_bamlist,
		adna_jar = adna_jar,
		picard_jar = picard_jar,
		python_release_libraries = python_release_libraries,
		python_read_groups_from_bam = python_read_groups_from_bam,
		release_parent_directory= release_parent_directory
	}
	call pulldown{ input:
		bamlist = nuclear_bamlist,
		report = report,
		python_pulldown = python_pulldown,
		python_release_libraries = python_release_libraries,
		python_read_groups_from_bam = python_read_groups_from_bam,
		release_parent_directory = release_parent_directory,
		label = label,
		unused = build_nuclear.unused
	}
}

task build_release_libraries{
	File bamlist
	
	File adna_jar
	File picard_jar
	File python_release_libraries
	File python_read_groups_from_bam
	
	String release_parent_directory
	
	Int processes = 8
	
	command{
		python3 ${python_release_libraries} --num_threads ${processes} --release_directory ${release_parent_directory} ${bamlist} ${adna_jar} ${picard_jar} > out
	}
	output{
		File unused = "out"
	}
	runtime{
		cpus: processes
		requested_memory_mb_per_core: 6000
	}
}

task pulldown{
	File bamlist
	File report
	#File pulldown_binary
	File python_pulldown
	File python_release_libraries
	File python_read_groups_from_bam
	
	String release_parent_directory
	String label
	
	File unused
	
	command{
		python3 ${python_pulldown} --pulldown_label ${label} --release_directory ${release_parent_directory} ${bamlist} ${report}
	}
	output{
		#normal_ind = "${label}.normal.ind"
		#normal_parameters = "${label}.normal.parameters"
		#damage_restricted_ind = "${label}.damage_restricted.ind"
		#damage_restricted_parameters = "${label}.damage_restricted.parameters"
		File dblist = "${label}.dblist"
	}
	runtime{
		requested_memory_mb_per_core: 6000
	}
}
