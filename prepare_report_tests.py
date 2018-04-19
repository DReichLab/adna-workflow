import unittest

from prepare_report import findSampleSheetEntry, parse_index_barcode_key_into_labels, sequence_labels, sequence_to_label

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
		
	def test_barcode_labels(self):
		barcode_filename = 'test/barcodes_q_only'
		barcode_mapping = sequence_labels(barcode_filename)
		
		# ATCGATT:CAGTCAA:GCTAGCC:TGACTGG	Q1	1	2	3	4
		self.assertEqual('Q1', barcode_mapping['ATCGATT:CAGTCAA:GCTAGCC:TGACTGG'])
		self.assertEqual('1', barcode_mapping['ATCGATT'])
		self.assertEqual('2', barcode_mapping['CAGTCAA'])
		self.assertEqual('3', barcode_mapping['GCTAGCC'])
		self.assertEqual('4', barcode_mapping['TGACTGG'])
		
		# GGTCATT:TTAGCAA:AACTGCC:CCGATGG	Q30	157	158	159	160
		self.assertEqual('Q30', barcode_mapping['GGTCATT:TTAGCAA:AACTGCC:CCGATGG'])
		self.assertEqual('157', barcode_mapping['GGTCATT'])
		self.assertEqual('158', barcode_mapping['TTAGCAA'])
		self.assertEqual('159', barcode_mapping['AACTGCC'])
		self.assertEqual('160', barcode_mapping['CCGATGG'])
		
	def test_i5_labels(self):
		i5_filename = 'test/i5_Reich20170725'
		i5_mapping = sequence_labels(i5_filename)
		#GGACGCA	10
		self.assertEqual('10', i5_mapping['GGACGCA'])
		self.assertEqual('10', sequence_to_label('GGACGCA', i5_mapping))
		#CGCCGTC	48
		self.assertEqual('48', i5_mapping['CGCCGTC'])
		self.assertEqual('48', sequence_to_label('CGCCGTC', i5_mapping))
		
		not_in_set = 'AAAAAAA'
		self.assertEqual(not_in_set, sequence_to_label(not_in_set, i5_mapping))
	
	def test_i7_labels(self):
		i7_filename = 'test/i7_Reich20170725'
		i7_mapping = sequence_labels(i7_filename)
		# TCGCAGG	1
		self.assertEqual('1', i7_mapping['TCGCAGG'])
		self.assertEqual('1', sequence_to_label('TCGCAGG', i7_mapping))
		# CGAGATC	81
		self.assertEqual('81', i7_mapping['CGAGATC'])
		self.assertEqual('81', sequence_to_label('CGAGATC', i7_mapping))
		# CCGTTGA	96
		self.assertEqual('96', i7_mapping['CCGTTGA'])
		self.assertEqual('96', sequence_to_label('CCGTTGA', i7_mapping))
		
	def test_index_barcode_labeling(self):
		barcode_filename = 'test/barcodes_q_only'
		barcode_mapping = sequence_labels(barcode_filename)
		i5_filename = 'test/i5_Reich20170725'
		i5_mapping = sequence_labels(i5_filename)
		i7_filename = 'test/i7_Reich20170725'
		i7_mapping = sequence_labels(i7_filename)
		
		key = 'GGACGCA_CGAGATC_ATCGATT:CAGTCAA:GCTAGCC:TGACTGG_GGTCATT:TTAGCAA:AACTGCC:CCGATGG'
		i5, i7, p5, p7 = parse_index_barcode_key_into_labels(key, i5_mapping, i7_mapping, barcode_mapping)
		self.assertEqual('10', i5)
		self.assertEqual('81', i7)
		self.assertEqual('Q1', p5)
		self.assertEqual('Q30', p7)
		
	def test_index_barcode_labeling_single_barcodes(self):
		barcode_filename = 'test/barcodes_q_only'
		barcode_mapping = sequence_labels(barcode_filename)
		i5_filename = 'test/i5_Reich20170725'
		i5_mapping = sequence_labels(i5_filename)
		i7_filename = 'test/i7_Reich20170725'
		i7_mapping = sequence_labels(i7_filename)
		
		key = 'GGACGCA_CGAGATC_ATCGATT_TTAGCAA'
		i5, i7, p5, p7 = parse_index_barcode_key_into_labels(key, i5_mapping, i7_mapping, barcode_mapping)
		self.assertEqual('10', i5)
		self.assertEqual('81', i7)
		self.assertEqual('1', p5)
		self.assertEqual('158', p7)
		
	def test_index_barcode_labeling_not_in_set(self):		
		key = 'AAAAAAA_CCCCCCC_GGGGGGG_TTTTTTT'
		i5, i7, p5, p7 = parse_index_barcode_key_into_labels(key, {}, {}, {})
		self.assertEqual('AAAAAAA', i5)
		self.assertEqual('CCCCCCC', i7)
		self.assertEqual('GGGGGGG', p5)
		self.assertEqual('TTTTTTT', p7)
		
	def test_index_barcode_labeling_empty(self):		
		key = ''
		i5, i7, p5, p7 = parse_index_barcode_key_into_labels(key, {}, {}, {})
		self.assertEqual('', i5)
		self.assertEqual('', i7)
		self.assertEqual('', p5)
		self.assertEqual('', p7)
		
	def test_index_barcode_labeling_none(self):
		i5, i7, p5, p7 = parse_index_barcode_key_into_labels(None, {}, {}, {})
		self.assertEqual('', i5)
		self.assertEqual('', i7)
		self.assertEqual('', p5)
		self.assertEqual('', p7)
		
	
if __name__ == '__main__':
	unittest.main()
