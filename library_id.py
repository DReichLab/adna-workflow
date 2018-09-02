import re

class LibraryID:
	# parse a library id string into components
	# examples:
	#	S1.E1.L1
	#	S2.L1
	#	S3.Y1.E1.L1
	def __init__(self, library_id_string):
		self.sample_suffix = ''
		self.lysis = None
		self.extract = None
		
		fields = library_id_string.split('.')
		# parse sample
		m = re.match('S(\d+)([a-zA-Z]*)', fields[0])
		if m != None:
			self.sample = int(m.group(1))
			self.sample_suffix = m.group(2)
		else:
			raise ValueError('Unable to parse sample: ' + fields[0])
		# parse library
		if fields[-1].startswith('L'):
			self.library = int(fields[-1][1:])
		else:
			raise ValueError('Library not found last: ' + fields[-1])
		
		for field in fields[1:-1]:
			if field.startswith('Y'):
				self.lysis = int(field[1:])
			elif field.startswith('E'):
				self.extract = int(field[1:])
			else:
				raise ValueError('Unexpected in library ID: ' + field)
			
	def __str__(self):
		fields = []
		# samples were initially padded to 4 digits
		if self.sample < 10000:
			fields.append("S{:04d}{}".format(self.sample, self.sample_suffix))
		else: # later samples are not padded more
			fields.append("S{:d}{}".format(self.sample, self.sample_suffix))
			
		if self.lysis:
			fields.append("Y{:d}".format(self.lysis))
		if self.extract:
			fields.append("E{:d}".format(self.extract))
		fields.append("L{:d}".format(self.library))
		
		return '.'.join(fields)
