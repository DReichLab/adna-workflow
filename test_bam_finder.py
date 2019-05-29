import unittest
from bam_finder import ShopVersion, getBamPath, library_default_dir, default_bam_root, MT_default_dir, bamInThisDirectory, getShopBamPath, find_pipeline_bam

class TestShopVersion(unittest.TestCase):
	
	def testValid(self):
		version_string = 'v0030.2__2018_05_02'
		self.assertTrue(ShopVersion.isValidVersionString(version_string))
		version = ShopVersion(version_string)
		self.assertEqual(30, version.major)
		self.assertEqual(2, version.minor)
		
	def testValidMT(self):
		version_string = 'MT.v0002.0__2018_03_08'
		self.assertTrue(ShopVersion.isValidVersionString(version_string))
		version = ShopVersion(version_string)
		self.assertEqual(2, version.major)
		self.assertEqual(0, version.minor)
		
	def testInvalid(self):
		self.assertFalse(ShopVersion.isValidVersionString('merged'))

	# this is the shop MT
	def testShopBamInDirectory(self):
		path = '/n/data1/hms/genetics/reich/1000Genomes/amh_samples/ancientMergeSets__MT/B-per_library_versions/S8861.E1.L1/MT.v0030.1__2018_08_08/merged/'
		expected = path + 'aln.sort.mapped.rmdupse_adna_v2.md.bam'
		found = bamInThisDirectory(path, default_bam_root)
		self.assertEqual(found, expected)

	# also Shop MT
	def testShopBamMT(self):
		expected = '/n/data1/hms/genetics/reich/1000Genomes/amh_samples/ancientMergeSets__MT/B-per_library_versions/S8861.E1.L1/MT.v0030.1__2018_08_08/merged/aln.sort.mapped.rmdupse_adna_v2.md.bam'
		found = getShopBamPath('S8861.E1.L1', MT_default_dir, default_bam_root)
		self.assertEqual(found, expected)

	def testShopNuclear(self):
		expected = '/n/data1/hms/genetics/reich/1000Genomes/amh_samples/ancientMergeSets__CAPTURE/B-per_library_versions/S8861.E1.L1/v0029.8__2018_03_01/merged/aln.sort.mapped.rmdupse_adna_v2.md.bam'
		found = getBamPath('S8861.E1.L1', library_default_dir, default_bam_root, 'hg19')
		self.assertEqual(expected, found)

	def testShopMT(self):
		expected = '/n/data1/hms/genetics/reich/1000Genomes/amh_samples/ancientMergeSets__MT/B-per_library_versions/S8861.E1.L1/MT.v0030.1__2018_08_08/merged/aln.sort.mapped.rmdupse_adna_v2.md.bam'
		found = getBamPath('S8861.E1.L1', MT_default_dir, default_bam_root, 'rsrs')
		self.assertEqual(expected, found)

	def testPipelineNuclear(self):
		expected = '/n/groups/reich/matt/pipeline/released_libraries/S8861.E1.L2/S8861.E1.L2.1240k_plus.hg19.v1.bam'
		found = str(getBamPath('S8861.E1.L2', '', '', 'hg19'))
		self.assertEqual(expected, found)

	def testPipelineMT(self):
		expected = '/n/groups/reich/matt/pipeline/released_libraries/S8861.E1.L2/S8861.E1.L2.1240k_plus.rsrs.v1.bam'
		found = str(getBamPath('S8861.E1.L2', '', '', 'rsrs'))
		self.assertEqual(expected, found)

	def testAlternateBigYoruba(self):
		expected = '/n/data1/hms/genetics/reich/1000Genomes/amh_samples/ancientMergeSets__CAPTURE/B-per_library_versions/S6441.E1.L2/v0025.1__2017_08_26/merged/aln.sort.mapped.rmdupse_adna_v2.md.bam'
		#former = '/n/data1/hms/genetics/reich/1000Genomes/amh_samples/ancientMergeSets__CAPTURE/B-per_library_versions/S6441.E1.L2/v0029.5b__2018_02_27__bigyoruba/merged/aln.sort.mapped.rmdupse_adna_v2.md.bam'
		# we do not want the BigYoruba version
		found = str(getShopBamPath('S6441.E1.L2', library_default_dir, default_bam_root))
		self.assertEqual(expected, found)
	
	# versions
	def testOnlyBAM_default(self):
		expected = 'S10123.E1.L1/S10123.E1.L1.1240k_plus.hg19.v1.bam'
		found = find_pipeline_bam('S10123.E1.L1', 'hg19', '1240k', pipeline_parent_bam_dir='bam_finder_test')
		self.assertTrue(str(found).endswith(expected))
	
	def testOnlyBAM(self):
		expected = 'S10123.E1.L1/S10123.E1.L1.1240k_plus.hg19.v1.bam'
		found = find_pipeline_bam('S10123.E1.L1', 'hg19', '1240k', version_policy='only', pipeline_parent_bam_dir='bam_finder_test')
		self.assertTrue(str(found).endswith(expected))
			
	def testOnlyBAMOfMany(self):
		found = find_pipeline_bam('S10124.E1.L1', 'hg19', '1240k', version_policy='only', pipeline_parent_bam_dir='bam_finder_test')
		self.assertEqual('', found)
		
	def testLatestBAMOne(self):
		expected = 'S10123.E1.L1/S10123.E1.L1.1240k_plus.hg19.v1.bam'
		found = find_pipeline_bam('S10123.E1.L1', 'hg19', '1240k', version_policy='latest', pipeline_parent_bam_dir='bam_finder_test')
		self.assertTrue(str(found).endswith(expected))
		
	def testLatestBAMMany(self):
		expected = 'S10124.E1.L1/S10124.E1.L1.1240k_plus.hg19.v2.bam'
		found = find_pipeline_bam('S10124.E1.L1', 'hg19', '1240k', version_policy='latest', pipeline_parent_bam_dir='bam_finder_test')
		self.assertTrue(str(found).endswith(expected))
		
if __name__ == '__main__':
	unittest.main()
