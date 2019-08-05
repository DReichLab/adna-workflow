#!/usr/bin/python3
import argparse
import os
import re
import sys
from pathlib import Path
from operator import attrgetter
from functools import reduce

from library_id import LibraryID

ONLY = 'only'
LATEST = 'latest'

# Shop's versioning strings look like v0030.2__2018_05_02
class ShopVersion():
	def __init__(self, version_directory):
		#m = re.match('(MT.)?v(\d+)\.(\d+)_', version_directory)
		m = re.search('(MT.)?v(\d+)\.(\d+)__(\d+)_(\d+)_(\d+)', version_directory)
		self.major = int(m.group(2))
		self.minor = int(m.group(3))
		self.directory = version_directory

		year = int(m.group(4))
		month = int(m.group(5))
		day = int(m.group(6))
		self.date_string = '{:04d}{:02d}{:02d}'.format(year, month, day)
		
	def __str__(self):
		return self.directory
	
	def isValidVersionString(string):
		try:
			version = ShopVersion(string)
			return True
		except:
			return False

# given an ID string, return a bam path
# This looks for a subdirectory of the parent_directory (either Shop's library or sample directories)
def getShopBamPath(requestedID, parent_directory, bam_root):
	pathname = Path(str(requestedID))
	parent_path = Path(parent_directory)
	
	path = parent_path / pathname
	if os.path.exists(path):
		# we expect to see a version directory
		with os.scandir(path) as top_directory:
			version_directories = [ShopVersion(x.name) for x in top_directory if x.is_dir() and ShopVersion.isValidVersionString(x.name) and (x.name.startswith('v') or x.name.startswith('MT.v'))]
			sorted_version_directories = sorted(version_directories, key=attrgetter('major', 'minor'))
			# find the most recent version
			for version_directory in reversed(sorted_version_directories):
				version_path = path / version_directory.directory
				bam = getBamInVersionDirectory(version_path, bam_root)
				if bam != '':
					return bam
	return ''

# If we already know the version directory, look for a bam within its merged directory
def getBamInVersionDirectory(path, bam_root):
	path = Path(path)
	# we expect to see a merged directory
	with os.scandir(path) as version_directory:
		merged_subdirectories = [x for x in version_directory if os.path.exists(x) and x.is_dir() and x.name == 'merged']
		if len(merged_subdirectories) == 1:
			path = path / merged_subdirectories[0].name
			return bamInThisDirectory(path, bam_root)
	return ''

def bamInThisDirectory(path, bam_root):
	try:
		with os.scandir(path) as final_directory:
			# find .bam file that starts with the specified root
			exact = [x for x in final_directory if isPossibleBAM(x) and x.name == (bam_root + '.bam')]
			if len(exact) == 1:
				return exact[0].path
			
			bams = [x for x in final_directory if isPossibleBAM(x) and x.name.startswith(bam_root)]
			if len(bams) == 1:
				return bams[0].path
			#elif len(bams) > 1:
			#	print('multiple possible bams in {}'.format(path), file=sys.stderr)
	except FileNotFoundError:
		pass
	return ''

def bamFromAnnoPath(anno_bam_path, bam_root):
	# fix paths from orchestra for O2
	if anno_bam_path.startswith('/groups/'):
		anno_bam_path = '/n' + anno_bam_path
	# check if this is a bam file
	if isPossibleBAM(anno_bam_path):
		return anno_bam_path
	else:
		# if this is a directory and all bams in it are symbolic links to another bam, use that bam
		try:
			with os.scandir(anno_bam_path) as directory:
				links = [Path(x) for x in directory if x.is_symlink() and x.name.endswith('.bam')]
				if len(links) > 0:
					
						candidateBAM = reduce(lambda x, y: x.resolve() if x.resolve() == y.resolve() else None, links)
						if candidateBAM != None and isPossibleBAM(candidateBAM):
							return candidateBAM
		except:
			pass # if any link does not resolve, give up and proceed to another option 
		# this is a directory; walk back up the path and look directly for a bam_candidate
		current_path = anno_bam_path
		while len(current_path) > 1 and bamInThisDirectory(current_path, bam_root) == '':
			current_path = os.path.dirname(current_path)
		if bamInThisDirectory(current_path, bam_root) != '':
			return bamInThisDirectory(current_path, bam_root)
		
		'''
		# this is a directory; walk back up the path to look for a version directory, and then a bam
		current_path = anno_bam_path
		while len(current_path) > 1 and not ShopVersion.isValidVersionString(os.path.basename(current_path)):
			current_path = os.path.dirname(current_path)
		# if we have found a version directory, look for a bam_candidate
		anno_bam_paths[instance_id] = ''
		if ShopVersion.isValidVersionString(os.path.basename(current_path)):
			anno_bam_paths[instance_id] = getBamInVersionDirectory(current_path, bam_root)
		'''
	return ''

def isPossibleBAM(bam_candidate):
	 return os.path.exists(bam_candidate) and os.path.isfile(bam_candidate) and Path(bam_candidate).suffix == '.bam'
 
# find a bam in pipeline results
def find_pipeline_bam(library_id, reference, experiment, version_policy=ONLY, pipeline_parent_bam_dir = '/n/groups/reich/matt/pipeline/released_libraries'):
	path = Path(pipeline_parent_bam_dir) / library_id
	if path.exists() and path.is_dir():
		with os.scandir(path) as directory:
			#print('{}\t{}\t{}\t{}'.format(library_id, reference, experiment, version_policy))
			candidates = [Path(x) for x in directory if x.name.startswith(library_id) and x.name.endswith('.bam') and reference in x.name and experiment in x.name]
			
			if len(candidates) == 1:
				return candidates[0]
			elif len(candidates) > 1:
				if version_policy == LATEST:
					latest = max(candidates, key=lambda p: int(p.name.split('.')[-2][1:])) # filename ends with 'v2.bam', parse out 2 from v2
					return latest
				elif version_policy == ONLY:
					raise ValueError('There are {:d} versions, not one'.format(len(candidates)))
	return ''
 
def readAnnoFile(anno_filename, requestedIDDict, bam_root):
	anno_bam_paths = {}
	read_groups_by_id = {}
	with open(anno_filename, errors='surrogateescape') as anno_file:
		headers = anno_file.readline().split('\t')
		num_headers = len(headers)
		# find header indices, and map them to correct IDs
		instance_id_index = headers.index('Instance ID')
		anno_bam_index = headers.index('pulldown: 3rd column of nickdb (bam)')
		read_groups_index = headers.index('pulldown: 5th column of nickdb (diploid source)')

		for line in anno_file:
			fields = line.split('\t')
			instance_id = fields[instance_id_index]
			anno_bam_path = fields[anno_bam_index]
			read_groups = fields[read_groups_index]
			
			if instance_id in requestedIDDict:
				found_anno_bam_path = bamFromAnnoPath(anno_bam_path, bam_root)
				if found_anno_bam_path != '':
					anno_bam_paths[instance_id] = found_anno_bam_path
					
			read_groups_by_id[instance_id] = read_groups
						
	return anno_bam_paths, read_groups_by_id

library_default_dir = '/n/data1/hms/genetics/reich/1000Genomes/amh_samples/ancientMergeSets__CAPTURE/B-per_library_versions'
sample_default_dir = '/n/data1/hms/genetics/reich/1000Genomes/amh_samples/ancientMergeSets__CAPTURE/C-per_sample_versions'
MT_default_dir = '/n/data1/hms/genetics/reich/1000Genomes/amh_samples/ancientMergeSets__MT/B-per_library_versions'
default_bam_root = 'aln.sort.mapped.rmdupse_adna_v2.md'

# TODO parent directory should have to match reference for proper lookup of Shop's bams
def getBamPath(requestedID, shop_parent_directory=library_default_dir, bam_root=default_bam_root, reference='hg19', experiment='1240k', version_policy=ONLY):
	shop_bam_path = getShopBamPath(requestedID, shop_parent_directory, bam_root)
	pipeline_bam_path = find_pipeline_bam(requestedID, reference, experiment, version_policy=version_policy)
	
	if shop_bam_path != '' and pipeline_bam_path != '':
		#raise ValueError('multiple bams for {}'.format(requestedID))
		print('multiple bams for {}'.format(requestedID), file=sys.stderr)
		return pipeline_bam_path
	elif shop_bam_path != '':
		return shop_bam_path
	elif pipeline_bam_path != '':
		return pipeline_bam_path
	else:
		return ''

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Try to find the bam associated with ID(s)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	parser.add_argument("-a", "--anno", help="Use anno file for bam hints and read groups")
	parser.add_argument("-f", "--filename", help="File containing IDs to lookup")
	parser.add_argument("-p", "--parent_directory", help="Parent directory to look for Shop's files. Use full path to get full path results.", default=library_default_dir)
	parser.add_argument("-b", "--bam_root", help="Bam filename root to look for", default=default_bam_root)
	parser.add_argument("-l", "--library_filter", help="Restrict file inputs to Reich Lab library ID format, for example 'S0123.E1.L2'", action='store_true')
	parser.add_argument("-r", "--reference", help="For example: hg19, rsrs", default='hg19')
	parser.add_argument("-e", "--experiment", help="Examples: 1240k, BigYoruba", default='1240k')
	parser.add_argument("--version_policy", choices=[ONLY, LATEST], default=LATEST, help='Policy for pipeline bams. Only will raise exception if there is more than one version.')
	parser.add_argument("requested_ids", help="Individual IDs to process from command line", nargs='*')
	args = parser.parse_args()

	requestedIDs = []

	if args.filename:
		#print ('opening {}'.format(args.filename))
		with open(args.filename) as f:
			for line in f:
				if args.library_filter: # restrict to library format
					try:
						libraryID = LibraryID(line)
						requestedIDs.append(libraryID)
					except:
						pass
				else:
					requestedIDs.append(line.strip())
	
	requestedIDs.extend(args.requested_ids)
	
	parent_directory = args.parent_directory
	if args.reference == 'rsrs':
		parent_directory = MT_default_dir
	
	requestedIDDict = {x : x for x in requestedIDs}
	
	bam_paths = {}
	for requestedID in requestedIDs:
		# remove trailing _all, _d, _published
		requestedID_shortened = requestedID.split('_')[0]
		bam_paths[requestedID] = getBamPath(requestedID_shortened, parent_directory, args.bam_root, args.reference, args.experiment, version_policy=args.version_policy)
	
	anno_bam_paths = {}
	read_groups_by_id = {}
	if args.anno:
		anno_bam_paths, read_groups_by_id = readAnnoFile(args.anno, requestedIDDict, args.bam_root)
	
	for requestedID in requestedIDs:
		if args.anno:
			print('{0}\t{0}\t{1}\t{2}'.format(requestedID, anno_bam_paths.get(requestedID, ''), read_groups_by_id.get(requestedID, '')))
			anno_bam_str = anno_bam_paths.get(requestedID, '')
			anno_bam_path = Path(anno_bam_str)
			bam_str = bam_paths.get(requestedID, '')
			bam_path = Path(bam_str)
			if requestedID not in anno_bam_paths and requestedID not in read_groups_by_id:
				print('{} does not match an Instance ID in anno file'.format(requestedID), file=sys.stderr)
			elif anno_bam_str == ''  or not os.path.exists(anno_bam_path):
				print('anno bam not found for {}'.format(requestedID), file=sys.stderr)
			if ((os.path.exists(anno_bam_path) and anno_bam_str != '' and os.path.exists(bam_path) and bam_str != '' and not os.path.samefile(anno_bam_path, bam_path))
				or (not os.path.exists(anno_bam_path)) and os.path.exists(bam_path)):
				print('mismatch between anno bam and latest bam available: {}\t[{}]\t[{}]'.format(requestedID, anno_bam_path, bam_path), file=sys.stderr)
		else:
			#print('{0}\t{0}\t{1}\t{2}'.format(requestedID, bam_paths.get(requestedID, ''), read_groups_by_id.get(requestedID, '')))
			print('{}'.format(bam_paths.get(requestedID))) 
