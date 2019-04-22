import argparse
from pathlib import Path
import re
from collections import OrderedDict

class LaneFiles:
	def __init__(self):
		self.I1 = None
		self.I2 = None
		self.R1 = None
		self.R2 = None
	
	# this is consumed by WDL read_tsv()
	# this order matches the inputs for merge_and_trim_lane and barcode_count 
	def __repr__(self):
		return '{}\t{}\t{}\t{}'.format(self.R1, self.R2, self.I1, self.I2)
		
	def is_complete(self):
		return self.I1 != None and self.I2 != None and self.R1 != None and self.R2 != None

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Group fastq files by lane for index-barcode analysis and merging", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("fastq", help="fastq files: index and regular reads", nargs='+')
	args = parser.parse_args()
	
	
	files_by_lane = OrderedDict()
	for fullpath in args.fastq:
		# fullpaths are expected to look like this Undetermined_S0_L001_I1_001.fastq
		m = re.match('Undetermined_S0_L(\d+)_([IR]\d)_001.fastq', Path(fullpath).name)
		lane_number = int(m.group(1))
		read_type = m.group(2)
		#print('{}\t{:d}\t{}'.format(fullpath, lane_number, read_type))
		
		if lane_number not in files_by_lane:
			files_by_lane[lane_number] = LaneFiles()
			
		if read_type == 'I1':
			files_by_lane[lane_number].I1 = fullpath
		elif read_type == 'I2':
			files_by_lane[lane_number].I2 = fullpath
		elif read_type == 'R1':
			files_by_lane[lane_number].R1 = fullpath
		elif read_type == 'R2':
			files_by_lane[lane_number].R2 = fullpath
		else:
			raise ValueError('Unhandled read type {}'.format(read_type))
		
	# print out data by lane numbers
	for lane_number, lane_data in files_by_lane.items():
		if lane_data.is_complete():
			print(str(lane_data))
		else: # all lanes should have complete data
			raise ValueError('Lane {:d} is incomplete'.format(lane_number))
