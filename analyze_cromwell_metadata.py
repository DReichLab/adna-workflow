# print out the elapsed time for tasks
# input is the cromwell metadata file
import json
import sys
from datetime import datetime

filename = sys.argv[1]

metadata = json.load(open(filename))
calls = metadata['calls']

# when seconds are whole, we need a separate format to read the datetime
# TODO handle timezones more robustly
datetimeFormat = "%Y-%m-%dT%H:%M:%S.%f-05:00"
datetimeFormat_alternate = "%Y-%m-%dT%H:%M:%S-05:00"

for task_name in calls:
	#print (task_name)
	task_attributes_by_shard = calls[task_name]
	for shard_num in range(0, len(task_attributes_by_shard)):
		task_attributes = task_attributes_by_shard[shard_num]
		
		# read the start and end times, falling back to alternate datetime format
		try:
			startTime = datetime.strptime(task_attributes['start'], datetimeFormat)
		except ValueError:
			startTime = datetime.strptime(task_attributes['start'], datetimeFormat_alternate)
		try:	
			endTime = datetime.strptime(task_attributes['end'], datetimeFormat)
		except ValueError:
			endTime = datetime.strptime(task_attributes['end'], datetimeFormat_alternate)
			
		elapsedSeconds = (endTime - startTime).total_seconds()
			
		executionEvents = task_attributes['executionEvents'][0]
		for event in executionEvents:
			elapsedRuntime = 0
			if executionEvents['description'] == 'RunningJob':
				try:
					runStartTime = datetime.strptime(executionEvents['startTime'], datetimeFormat)
				except ValueError:
					runStartTime = datetime.strptime(executionEvents['startTime'], datetimeFormat_alternate)
				try:	
					runEndTime = datetime.strptime(executionEvents['endTime'], datetimeFormat)
				except ValueError:
					runEndTime = datetime.strptime(executionEvents['endTime'], datetimeFormat_alternate)
				elapsedRuntime = (endTime - runStartTime).total_seconds()
				break
		
		# print task name and leave out the workflow name
		print ("{}\t{}\t{}".format(task_name.split('.')[1], elapsedSeconds, elapsedRuntime))
