import "demultiplex.wdl" as demultiplex_align_bams

workflow prepare_all_references{
	File reference_in
	File mt_reference_rsrs_in
	File mt_reference_rcrs_in

	call demultiplex_align_bams.prepare_reference as prepare_reference_nuclear{ input:
		reference = reference_in
	}
	call demultiplex_align_bams.prepare_reference as prepare_reference_rsrs{ input:
		reference = mt_reference_rsrs_in
	}
	call demultiplex_align_bams.prepare_reference as prepare_reference_rcrs{ input:
		reference = mt_reference_rcrs_in
	}
}
