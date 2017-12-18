import unittest

from prepare_report import findSampleSheetEntry

class TestPrepareReport(unittest.TestCase):
#	CCTGCGA 1 I5
#	TCGCAGG 1 I7
#	GGTATCG:TTACAGT:AACGCTA:CCGTGAC Q2      13      14      15      16
#	TACGTTC:ACGTAAG:CGTACCT:GTACGGA Q3      25      26      27      28
	def test_simpleLookup(self):
		full = 'CCTGCGA_TCGCAGG_GGTATCG:TTACAGT:AACGCTA:CCGTGAC_TACGTTC:ACGTAAG:CGTACCT:GTACGGA'
		expectedLibraryID = 'test'
		keyMapping = {full: expectedLibraryID}
		
		sampleSheetID, libraryID = findSampleSheetEntry(full, keyMapping)
		self.assertEqual(full, sampleSheetID)
		self.assertEqual(expectedLibraryID, libraryID)
		
	def test_subsetLookup(self):
		full = 'CCTGCGA_TCGCAGG_GGTATCG:TTACAGT:AACGCTA:CCGTGAC_TACGTTC:ACGTAAG:CGTACCT:GTACGGA'
		singleBarcodes = 'CCTGCGA_TCGCAGG_CCGTGAC_ACGTAAG'
		expectedLibraryID = 'test'
		keyMapping = dict()
		keyMapping[singleBarcodes] = expectedLibraryID
		
		sampleSheetID, libraryID = findSampleSheetEntry(full, keyMapping)
		self.assertEqual(singleBarcodes, sampleSheetID)
		self.assertEqual(expectedLibraryID, libraryID)
		
	def test_noMatch(self):
		full = 'CCTGCGA_TCGCAGG_GGTATCG:TTACAGT:AACGCTA:CCGTGAC_TACGTTC:ACGTAAG:CGTACCT:GTACGGA'
		expectedLibraryID = 'test'
		keyMapping = dict()
		
		sampleSheetID, libraryID = findSampleSheetEntry(full, keyMapping)
		self.assertEqual('', sampleSheetID)
		self.assertEqual('', libraryID)
	
	def test_multiple_matches(self):
		full = 'CCTGCGA_TCGCAGG_GGTATCG:TTACAGT:AACGCTA:CCGTGAC_TACGTTC:ACGTAAG:CGTACCT:GTACGGA'
		singleBarcodes1 = 'CCTGCGA_TCGCAGG_CCGTGAC_ACGTAAG'
		singleBarcodes2 = 'CCTGCGA_TCGCAGG_CCGTGAC_GTACGGA'
		
		keyMapping = dict()
		keyMapping[singleBarcodes1] = 'test1'
		keyMapping[singleBarcodes2] = 'test2'
		
		expectedLibraryID = 'MULTIPLE'
		expectedSampleSheetID = 'MULTIPLE'
		
		sampleSheetID, libraryID = findSampleSheetEntry(full, keyMapping)
		self.assertEqual(expectedLibraryID, libraryID)
		self.assertEqual(expectedSampleSheetID, sampleSheetID)
	
if __name__ == '__main__':
	unittest.main()
