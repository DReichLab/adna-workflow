import unittest
from preseq_process import read_preseq_file, total_and_unique_target_hits, EmpiricalTargetEstimator, preseq_analysis

class TestPreseq(unittest.TestCase):
	preseq_table_filename = 'test/AGTTGGT_CGACCTG_ACGGTCT-CGTTAGA-GTAACTC-TACCGAG_CTAGACA-GACTCGC-TCGAGTG-AGTCTAT.preseq_table'
	histogram_filename = 'test/AGTTGGT_CGACCTG_ACGGTCT-CGTTAGA-GTAACTC-TACCGAG_CTAGACA-GACTCGC-TCGAGTG-AGTCTAT.targets_histogram'
	
	def test_read_table(self):
		
		reads_hitting_any_target, unique_reads = read_preseq_file(self.preseq_table_filename)
		
		self.assertEqual(0, reads_hitting_any_target[0])
		self.assertEqual(506042.0, reads_hitting_any_target[1])
		self.assertEqual(5060420.0, reads_hitting_any_target[10])
		
		self.assertEqual(0, unique_reads[0])
		self.assertEqual(500627.5, unique_reads[1])
		self.assertEqual(4589071.5, unique_reads[10])
		
	def test_total_and_unique_target_hits(self):
		total_hits, unique_targets = total_and_unique_target_hits(self.histogram_filename)
		self.assertEqual(2226266, total_hits)
		self.assertEqual(574843, unique_targets)

	def test_empirical_estimator_default(self):
		estimator = EmpiricalTargetEstimator(1, 1, -1)
		expected = 1.0 / 3.0 * 1150639
		estimate = estimator.unique_reads_to_empirical_targets(575319.5)
		self.assertAlmostEqual(expected, estimate)
		
	def test_empirical_estimator_custom(self):
		estimator = EmpiricalTargetEstimator(1.1, 0.9, -1)
		expected = 396772.068965517
		estimate = estimator.unique_reads_to_empirical_targets(575319.5)
		self.assertAlmostEqual(expected, estimate)
		
	def test_preseq_analysis(self):
		number_raw_reads = 11036945
		minimum_raw_reads = 3e6
		expected_targets_per_raw_read_threshold = 0.01
		reads_hitting_any_target, unique_reads = read_preseq_file(self.preseq_table_filename)
		empiricalTargetEstimator = EmpiricalTargetEstimator(1.1, 0.9, -1)
		total_hits, unique_targets = total_and_unique_target_hits(self.histogram_filename)
		
		values = preseq_analysis(reads_hitting_any_target, unique_reads, number_raw_reads, total_hits, unique_targets, minimum_raw_reads, expected_targets_per_raw_read_threshold, empiricalTargetEstimator)
		
		expected_raw_reads = 17561290.5429225
		self.assertAlmostEqual(number_raw_reads, values['number_raw_reads'], places=0)
		self.assertAlmostEqual(expected_raw_reads, values['preseq_total_reads_required'], places=0)
		self.assertAlmostEqual(expected_raw_reads - number_raw_reads, values['preseq_additional_reads_required'], places=0)

if __name__ == '__main__':
	unittest.main()
