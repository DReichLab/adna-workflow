import "demultiplex.wdl" as demultiplex
import "analysis.wdl" as analysis

workflow mt_haplogroup{
	call demultiplex.prepare_reference as prepare_reference_rcrs{
	}

	call analysis.haplogrep{ input:
		reference = prepare_reference_rcrs.reference_fa,
		reference_amb = prepare_reference_rcrs.reference_amb,
		reference_ann = prepare_reference_rcrs.reference_ann,
		reference_bwt = prepare_reference_rcrs.reference_bwt,
		reference_pac = prepare_reference_rcrs.reference_pac,
		reference_sa = prepare_reference_rcrs.reference_sa
	}
}
