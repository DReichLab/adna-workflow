import unittest
import tempfile
import os.path

from snp_target_coverage_bed import sex_determination

class TestSexDetermination(unittest.TestCase):
	
	def test_male(self):
		x = 100
		y = 100
		sex = sex_determination(x, y)
		self.assertEqual('M', sex)
	
	def test_male_threshold_met(self):
		x = 70
		y = 30
		sex = sex_determination(x, y)
		self.assertEqual('M', sex)
		
	def test_male_threshold_unmet(self):
		x = 71
		y = 30
		sex = sex_determination(x, y)
		self.assertEqual('U', sex)
	
	def test_female(self):
		x = 100
		y = 10
		sex = sex_determination(x, y)
		self.assertEqual('F', sex)
		
	def test_female_threshold_met(self):
		x = 90
		y = 10
		sex = sex_determination(x, y)
		self.assertEqual('F', sex)
		
	def test_female_threshold_unmet(self):
		x = 90
		y = 11
		sex = sex_determination(x, y)
		self.assertEqual('U', sex)
		
	def test_low_count(self):
		x = 10
		y = 10
		sex = sex_determination(x, y)
		self.assertEqual('U', sex)
		
	def test_ambiguous(self):
		x = 200
		y = 50
		sex = sex_determination(x, y)
		self.assertEqual('U', sex)
		
if __name__ == '__main__':
	unittest.main()
