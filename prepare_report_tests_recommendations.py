import unittest
import copy

from prepare_report import recommendation

LOW_DATA = 'low data'
FAIL = 'fail'
WARNING = 'warning'
PASS = 'pass'

class TestPrepareReport(unittest.TestCase):
	passing_sample = {
		'merged': 1000000,
		'damage_rsrs_ct1': 0.1,
		'damage_rsrs_ga1': 0.1,
		'contamination_contammix': 0.98,
		'MT_post-coverageLength': 200000,
		'spike3k_complexity': 0.2,
		}
	
	def test_sample_pass(self):
		self.assertEqual(PASS, recommendation(self.passing_sample))
		
	def test_low_data(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['merged'] = 200
		self.assertEqual(LOW_DATA, recommendation(test_sample))
		
	def test_fail_damage(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['damage_rsrs_ct1'] = 0.001
		test_sample['damage_rsrs_ga1'] = 0.001
		self.assertEqual(FAIL, recommendation(test_sample))
		
	def test_warn_damage1(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['damage_rsrs_ct1'] = 0.02
		test_sample['damage_rsrs_ga1'] = 0.02
		self.assertEqual(WARNING, recommendation(test_sample))
		
	def test_warn_damage2(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['damage_rsrs_ct1'] = 0.01
		test_sample['damage_rsrs_ga1'] = 0.04
		self.assertEqual(WARNING, recommendation(test_sample))
		
	def test_fail_contamination(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['contamination_contammix'] = 0.55
		self.assertEqual(FAIL, recommendation(test_sample))
		
	def test_warn_contamination(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['contamination_contammix'] = 0.9
		self.assertEqual(WARNING, recommendation(test_sample))
		
	def test_fail_coverage(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['MT_post-coverageLength'] = 2000
		print('fail coverage')
		self.assertEqual(FAIL, recommendation(test_sample))
		
	def test_warn_coverage(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['MT_post-coverageLength'] = 17000
		self.assertEqual(WARNING, recommendation(test_sample))
		
	def test_fail_spike_complexity(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['spike3k_complexity'] = 0.001
		self.assertEqual(FAIL, recommendation(test_sample))
	
	def test_warn_spike_complexity(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['spike3k_complexity'] = 0.01
		self.assertEqual(WARNING, recommendation(test_sample))
		
	def test_fail_damage_complexity(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['damage_rsrs_ct1'] = 0.01
		test_sample['damage_rsrs_ga1'] = 0.01
		test_sample['spike3k_complexity'] = 0.001
		self.assertEqual(FAIL, recommendation(test_sample))
	
	def test_warn_damage_fail_coverage(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['damage_rsrs_ct1'] = 0.02
		test_sample['damage_rsrs_ga1'] = 0.02
		test_sample['MT_post-coverageLength'] = 2000
		self.assertEqual(FAIL, recommendation(test_sample))
		
	def test_fail_contamination_warn_spike(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['contamination_contammix'] = 0.52
		test_sample['spike3k_complexity'] = 0.01
		self.assertEqual(FAIL, recommendation(test_sample))
		
	def test_pass_rescue(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['spike3k_complexity'] = 10
		self.assertEqual(PASS, recommendation(test_sample))
		
	def test_fail_contamination_rescue(self):
		test_sample = copy.deepcopy(self.passing_sample)
		test_sample['contamination_contammix'] = 0.59
		test_sample['spike3k_complexity'] = 5
		self.assertEqual(WARNING, recommendation(test_sample))
	
if __name__ == '__main__':
	unittest.main()
 
