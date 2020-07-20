import "demultiplex.wdl" as demultiplex
import "analysis.wdl" as analysis

workflow mt_contammix{
	Array[File] bams
	Int minimum_mapping_quality
	Int minimum_base_quality

	File adna_screen_jar
	File python_target
	File python_coverage

	call demultiplex.prepare_reference as prepare_reference_rsrs{
	}

	call analysis.chromosome_target as rsrs_chromosome_target_post{ input:
		python_target = python_target,
		adna_screen_jar = adna_screen_jar,
		bams = bams,
		targets="\"{'MT_post':'MT'}\"",
		minimum_mapping_quality = minimum_mapping_quality
	}

	call analysis.chromosome_coverage as rsrs_coverage{ input:
		bam_stats = rsrs_chromosome_target_post.target_stats,
		python_coverage = python_coverage,
		reference_length = 16569,
		coverage_field = "MT_post-coverageLength"
	}

	scatter(bam in bams){
		call analysis.contammix{ input:
			minimum_mapping_quality = minimum_mapping_quality,
			minimum_base_quality = minimum_base_quality,

			bam = bam,

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
}
