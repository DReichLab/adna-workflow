from multiprocessing import Pool
from os.path import basename
import subprocess
import os
import argparse

def merge_bam(instance_id, library_ids, experiments, bam_paths, reference, picard_jar):
	instance_id_filename = "%s.%s.bam" % (instance_id, reference)
	# make a directory for this instance ID
	os.makedirs(instance_id)
	with open(instance_id + '/stdout_merge', 'w') as stdout_merge, \
		open(instance_id + '/stderr_merge', 'w') as stderr_merge:
		# write instance ID into read groups
		bams_with_altered_read_groups = []
		for library_id, experiment, bam in zip(library_ids, experiments, bam_paths):
			read_group_id = library_id + '.' + experiment
			bam_with_altered_read_groups = instance_id + '/' + read_group_id + '.' + basename(bam)
			subprocess.run(["java", "-Xmx2700m", "-jar", picard_jar, "AddOrReplaceReadGroups", 
				"I=%s" % (bam,), 
				"O=%s" % (bam_with_altered_read_groups,), 
				"RGID=%s" % (read_group_id,), 
				"RGLB=%s" % (library_id,),
				"RGPL=illumina",
				"RGPU=%s" % (library_id,),
				"RGSM=%s" % instance_id], 
				check=True, stdout=stdout_merge, stderr=stderr_merge)
			bams_with_altered_read_groups.append(bam_with_altered_read_groups)
		# merge
		merge_file_list = 'I=' + ' I='.join(bams_with_altered_read_groups)
		# TODO output should be captured per instance
	command = "java -Xmx2500m -jar %s MergeSamFiles %s O=%s SORT_ORDER=coordinate >> stdout_merge 2>> stderr_merge" % (picard_jar, merge_file_list, instance_id_filename)
	#print('combine bam lists ' + command)
	subprocess.check_output(command, shell=True)
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Merge component library bams for a sample ")
	parser.add_argument('-n', "--num_threads", help="size of thread pool", type=int, default =1)
	parser.add_argument('-r', "--reference", help="size of thread pool", default='hg19')
	parser.add_argument("picard_jar", help="jar file for the Broad Institute Picard toolset, used to merge bams and deduplicate")
	parser.add_argument("instance_id", help="")
	parser.add_argument("bam_list", help="component libraries with experiment and paths")
	
	args = parser.parse_args()

	pool = Pool(processes=args.num_threads)
	results = []
	with open(args.bam_list) as f:
		library_ids = []
		experiments = []
		bam_paths = []
		for line in f:
			fields = line.split()
			library_ids.append(fields[0])
			experiments.append(fields[1])
			bam_paths.append(fields[2])
		results.append(pool.apply_async(merge_bam, args=(args.instance_id, library_ids, experiments, bam_paths, args.reference, args.picard_jar)))
	pool.close()
	pool.join()
	for result in results:
		result.get()
