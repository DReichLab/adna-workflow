import unittest
from angsd_results import parse_angsd_results

class TestAngsd(unittest.TestCase):
	def test_file_example(self):
		results = parse_angsd_results("angsd_contamination_output")
		self.assertEqual(results.get("nsites"), 2092)
		self.assertAlmostEqual(results.get("MoM"), 0.006946)
		self.assertAlmostEqual(results.get("SE(MoM)"), 2.138220e-03)
		self.assertAlmostEqual(results.get("ML"), 0.006724)
		self.assertAlmostEqual(results.get("SE(ML)"), 7.099542e-15)

if __name__ == '__main__':
	unittest.main()
