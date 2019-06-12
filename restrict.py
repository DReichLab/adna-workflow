import argparse
import sys

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Restrict an input file to those lines containing at least one line (for example, library ids) from a second file. ")
	parser.add_argument('-r', "--restrict", help="File(s) containing lines to restrict with", nargs='+', required=True)
	parser.add_argument('-t', '--target', help="Restrict this file to lines containing ")
	parser.add_argument('-u', '--unused', help="Print unused restrict values to stderr", action='store_true')
	args = parser.parse_args() 

	# build restrict list from all restrict input files
	restrict_list = []
	for restrict_file in args.restrict:
		with open(restrict_file) as f:
			for line in f:
				restrict_list.append(line.strip())
	
	# filter target file using restrict list
	used = {}
	with open(args.target) as f:
		for line in f:
			for candidate in restrict_list:
				if candidate in line:
					print(line, end='')
					used[candidate] = 1
					break

	unused = {}
	for x in restrict_list:
		if x not in used:
			unused[x] = 1
	if args.unused:
		for x in unused:
			print('{}\tunused'.format(x), file=sys.stderr)
