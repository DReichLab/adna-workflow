import unittest
import tempfile
import os

import merge_pulldown

class TestPreseq(unittest.TestCase):
	def test_duplicate_ind_one_file(self):
		input_ind_filenames = ['test/pulldown/duplicate.ind']
		with tempfile.TemporaryDirectory() as temp_directory:
			output_filename = os.path.join(temp_directory, 'duplicate_ind_one_file.ind')
			self.assertRaisesRegex(ValueError, "individual duplicate appears more than once", merge_pulldown.write_ind_file, output_filename, input_ind_filenames)
			
	def test_duplicate_ind_two_files(self):
		input_ind_filenames = ['test/pulldown/x.ind', 'test/pulldown/x.ind']
		with tempfile.TemporaryDirectory() as temp_directory:
			output_filename = os.path.join(temp_directory, 'duplicate_ind_two_files.ind')
			self.assertRaisesRegex(ValueError, "individual x appears more than once", merge_pulldown.write_ind_file, output_filename, input_ind_filenames)
			
	def test_bad_sex(self):
		input_ind_filenames = ['test/pulldown/bad_sex.ind']
		with tempfile.TemporaryDirectory() as temp_directory:
			output_filename = os.path.join(temp_directory, 'bad_sex.ind')
			self.assertRaisesRegex(ValueError, "invalid sex", merge_pulldown.write_ind_file, output_filename, input_ind_filenames)
			
	def test_bad_snp1(self):
		self.assertRaisesRegex(ValueError, "does not have expected 6 fields:", merge_pulldown.check_snp_file, 'test/pulldown/duplicate.ind')
		
	def test_bad_snp2(self):
		self.assertRaisesRegex(ValueError, "does not have expected 6 fields:", merge_pulldown.check_snp_file, 'test/pulldown/x.geno')
		
	def test_bad_snp_allele(self):
		self.assertRaisesRegex(ValueError, "bad alleles", merge_pulldown.check_snp_file, 'test/pulldown/bad_allele.snp')
		
	def test_bad_snp_position(self):
		self.assertRaisesRegex(ValueError, "bad position", merge_pulldown.check_snp_file, 'test/pulldown/bad_position.snp')
		
	def test_four_snp(self):
		snp_check = merge_pulldown.check_snp_file('test/pulldown/four.snp')
		self.assertTrue(snp_check)
		
	def test_snp_comparison(self):
		snp_filename = 'test/pulldown/four.snp'
		snp_filenames = [snp_filename, snp_filename]
		check = merge_pulldown.compare_snp_files(snp_filenames)
		self.assertTrue(check)
		
	def test_combine(self):
		ind_counts = [1, 2]
		num_snps = 4
		with tempfile.TemporaryDirectory() as temp_directory:
			output_filename = os.path.join(temp_directory, 'combined.geno')
			merge_pulldown.combine_geno_files(output_filename, ['test/pulldown/x.geno', 'test/pulldown/yz.geno'], ind_counts, num_snps)
			
			with open(output_filename, 'r') as combined, open('test/pulldown/xyz.geno', 'r') as expected:
				combined_lines = combined.readlines()
				expected_lines = expected.readlines()
				
				for i in range(num_snps):
					self.assertEqual(combined_lines[i].strip(), expected_lines[i].strip())
					
	def test_bad_geno(self):
		ind_counts = [1, 2]
		num_snps = 4
		with tempfile.TemporaryDirectory() as temp_directory:
			output_filename = os.path.join(temp_directory, 'combined.geno')
			self.assertRaisesRegex(ValueError, "invalid genotype value\(s\)", merge_pulldown.combine_geno_files, output_filename, ['test/pulldown/bad.geno', 'test/pulldown/yz.geno'], ind_counts, num_snps)
			
	def test_geno_too_many(self):
		ind_counts = [1, 2]
		num_snps = 4
		with tempfile.TemporaryDirectory() as temp_directory:
			output_filename = os.path.join(temp_directory, 'combined.geno')
			self.assertRaisesRegex(ValueError, "extra values in", merge_pulldown.combine_geno_files, output_filename, ['test/pulldown/too_many.geno', 'test/pulldown/yz.geno'], ind_counts, num_snps)

if __name__ == '__main__':
	unittest.main()
