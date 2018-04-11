import unittest

from prepare_report import findSampleSheetEntry

class TestPrepareReport(unittest.TestCase):
#	CCTGCGA 1 I5
#	TCGCAGG 1 I7
#	GGTATCG:TTACAGT:AACGCTA:CCGTGAC Q2      13      14      15      16
#	TACGTTC:ACGTAAG:CGTACCT:GTACGGA Q3      25      26      27      28
	def test_simpleLookup(self):
		full = 'CCTGCGA_TCGCAGG_GGTATCG:TTACAGT:AACGCTA:CCGTGAC_TACGTTC:ACGTAAG:CGTACCT:GTACGGA'
		expectedLibraryID = 'testLibrary'
		expectedPlateID = 'testPlate'
		expectedExperiment = 'test1240kPlus'
		keyMapping = {full: [expectedLibraryID, expectedPlateID, expectedExperiment]}
		
		sampleSheetID, libraryID, plateID, experiment = findSampleSheetEntry(full, keyMapping)
		self.assertEqual(full, sampleSheetID)
		self.assertEqual(expectedLibraryID, libraryID)
		self.assertEqual(expectedPlateID, plateID)
		self.assertEqual(expectedExperiment, experiment)
		
	def test_subsetLookup(self):
		full = 'CCTGCGA_TCGCAGG_GGTATCG:TTACAGT:AACGCTA:CCGTGAC_TACGTTC:ACGTAAG:CGTACCT:GTACGGA'
		singleBarcodes = 'CCTGCGA_TCGCAGG_CCGTGAC_ACGTAAG'
		expectedLibraryID = 'testLibrary'
		expectedPlateID = 'testPlate'
		expectedExperiment = 'test1240kPlus'
		keyMapping = dict()
		keyMapping = {singleBarcodes: [expectedLibraryID, expectedPlateID, expectedExperiment]}
		
		sampleSheetID, libraryID, plateID, experiment = findSampleSheetEntry(full, keyMapping)
		self.assertEqual(singleBarcodes, sampleSheetID)
		self.assertEqual(expectedLibraryID, libraryID)
		self.assertEqual(expectedPlateID, plateID)
		self.assertEqual(expectedExperiment, experiment)
		
	def test_noMatch(self):
		full = 'CCTGCGA_TCGCAGG_GGTATCG:TTACAGT:AACGCTA:CCGTGAC_TACGTTC:ACGTAAG:CGTACCT:GTACGGA'
		expectedLibraryID = 'testLibrary'
		expectedPlateID = 'testPlate'
		expectedExperiment = 'test1240kPlus'
		keyMapping = dict()
		
		sampleSheetID, libraryID, plateID, experiment = findSampleSheetEntry(full, keyMapping)
		self.assertEqual('', sampleSheetID)
		self.assertEqual('', libraryID)
		self.assertEqual('', plateID)
		self.assertEqual('', experiment)
	
	def test_multiple_matches(self):
		full = 'CCTGCGA_TCGCAGG_GGTATCG:TTACAGT:AACGCTA:CCGTGAC_TACGTTC:ACGTAAG:CGTACCT:GTACGGA'
		singleBarcodes1 = 'CCTGCGA_TCGCAGG_CCGTGAC_ACGTAAG'
		singleBarcodes2 = 'CCTGCGA_TCGCAGG_CCGTGAC_GTACGGA'
		
		keyMapping = dict()
		keyMapping[singleBarcodes1] = ['test1a', 'test1b', 'test1c']
		keyMapping[singleBarcodes2] = ['test2a', 'test2b', 'test2c']
		
		MULTIPLE = 'MULTIPLE'
		expectedLibraryID = MULTIPLE
		expectedSampleSheetID = MULTIPLE
		expectedPlateID = MULTIPLE
		expectedExperiment = MULTIPLE
		
		sampleSheetID, libraryID,plateID, experiment = findSampleSheetEntry(full, keyMapping)
		self.assertEqual(expectedLibraryID, libraryID)
		self.assertEqual(expectedSampleSheetID, sampleSheetID)
		self.assertEqual(expectedPlateID, plateID)
		self.assertEqual(expectedExperiment, experiment)
	
if __name__ == '__main__':
	unittest.main()
