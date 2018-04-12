import unittest
import tempfile
import os.path
from preseq_process import read_preseq_file, total_and_unique_target_hits, EmpiricalTargetEstimator, preseq_analysis, find_xy_for_slope

class TestPreseq(unittest.TestCase):
	# This is from Ellora 20180312, processed with step size = reads / 4
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
		
		self.assertAlmostEqual(0, values['preseq_raw_reads_inverse_e'], places=0)
		
		expected_raw_reads_tenth = 3991851.51647285
		self.assertAlmostEqual(expected_raw_reads_tenth, values['preseq_raw_reads_tenth'], places=0)
		
		expected_raw_reads_threshold = 18514435
		self.assertAlmostEqual(number_raw_reads, values['number_raw_reads'], places=0)
		self.assertAlmostEqual(expected_raw_reads_threshold, values['preseq_total_reads_required'], places=0)
		self.assertAlmostEqual(expected_raw_reads_threshold - number_raw_reads, values['preseq_additional_reads_required'], places=0)
		
		expected_targets_at_threshold = 823460.771239422
		self.assertAlmostEqual(expected_targets_at_threshold, values['preseq_expected_unique_targets_at_threshold'], places=0)
		
	def test_fail_to_open_preseq_file(self):
		reads_hitting_any_target, unique_reads = read_preseq_file('does_not_exist')
		self.assertEqual(0, len(reads_hitting_any_target))
		self.assertEqual(0, len(unique_reads))
		
	def test_empty_preseq_file(self):
		with tempfile.TemporaryDirectory() as temp_directory:
			empty_filename = os.path.join(temp_directory, 'empty_file_example')
			with open(empty_filename, 'w') as empty:
				 pass
			self.assertTrue(os.path.isfile(empty_filename))
			
			reads_hitting_any_target, unique_reads = read_preseq_file(empty_filename)
			self.assertIsNotNone(reads_hitting_any_target)
			self.assertEqual(0, len(reads_hitting_any_target))
			self.assertIsNotNone(unique_reads)
			self.assertEqual(0, len(unique_reads))
			
	def test_fail_to_open_target_histogram_file(self):
		total_hits, unique_targets = total_and_unique_target_hits('does_not_exist')
		self.assertEqual(0, total_hits)
		self.assertEqual(0, unique_targets)
			
	def test_empty_target_histogram_file(self):
		with tempfile.TemporaryDirectory() as temp_directory:
			empty_filename = os.path.join(temp_directory, 'empty_file_example')
			with open(empty_filename, 'w') as empty:
				 pass
			self.assertTrue(os.path.isfile(empty_filename))
			
			total_hits, unique_targets = total_and_unique_target_hits(empty_filename)
			self.assertEqual(0, total_hits)
			self.assertEqual(0, unique_targets)
			
	def test_find_xy_for_slope_empty(self):
		X = []
		Y = []
		
		x, y = find_xy_for_slope(X, Y, 1.0)
		self.assertAlmostEqual(0, x)
		self.assertAlmostEqual(0, y)
	
	def test_find_xy_for_slope_too_high(self):
		X = [0, 100, 200, 300, 400]
		Y = [0, 40, 60, 70, 75]
		
		x, y = find_xy_for_slope(X, Y, 1.0)
		self.assertAlmostEqual(0, x)
		self.assertAlmostEqual(0, y)
	
	def test_find_xy_for_slope_exact(self):
		X = [0, 100, 200, 300, 400]
		Y = [0, 40, 60, 70, 75]
		
		x, y = find_xy_for_slope(X, Y, 0.4)
		self.assertAlmostEqual(100, x)
		self.assertAlmostEqual(40, y)
		
		x, y = find_xy_for_slope(X, Y, 0.2)
		self.assertAlmostEqual(200, x)
		self.assertAlmostEqual(60, y)
		
		x, y = find_xy_for_slope(X, Y, 0.1)
		self.assertAlmostEqual(300, x)
		self.assertAlmostEqual(70, y)
		
		x, y = find_xy_for_slope(X, Y, 0.05)
		self.assertAlmostEqual(400, x)
		self.assertAlmostEqual(75, y)
	
	def test_find_xy_for_slope_interpolate(self):
		X = [0, 100, 200, 300, 400]
		Y = [0, 40, 60, 70, 75]
		
		x, y = find_xy_for_slope(X, Y, 0.3)
		self.assertAlmostEqual(150, x)
		self.assertAlmostEqual(55, y)
		
		x, y = find_xy_for_slope(X, Y, 0.15)
		self.assertAlmostEqual(250, x)
		self.assertAlmostEqual(67.5, y)
		
		x, y = find_xy_for_slope(X, Y, 0.35)
		self.assertAlmostEqual(125, x)
		self.assertAlmostEqual(48.75, y)
		
		x, y = find_xy_for_slope(X, Y, 0.25)
		self.assertAlmostEqual(175, x)
		self.assertAlmostEqual(58.75, y)
		
		x, y = find_xy_for_slope(X, Y, 0.075)
		self.assertAlmostEqual(350, x)
		self.assertAlmostEqual(73.75, y)
		
if __name__ == '__main__':
	unittest.main()
