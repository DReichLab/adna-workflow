import argparse
import sys

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Restrict an input file to those lines containing at least one line (for example, library ids) from a second file. ")
	parser.add_argument('-r', "--restrict", help="File(s) containing lines to restrict with", nargs='+', required=True)
	parser.add_argument('-t', '--target', help="Restrict this file to lines containing ")
	parser.add_argument('-u', '--unused', help="Print unused restrict values to stderr", action='store_true')
	parser.add_argument('-v', '--invert', help="Invert match to print only lines that do not match", action='store_true')
	parser.add_argument('-p', '--preserve_header', help="preserve header in output", action='store_true')
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
		if args.preserve_header:
			line = f.readline()
			print(line, end='')
		for line in f:
			match = False
			for candidate in restrict_list:
				if candidate in line:
					match = True
					used[candidate] = 1
					break
			if match != args.invert: #xor
				print(line, end='')

	unused = {}
	for x in restrict_list:
		if x not in used:
			unused[x] = 1
	if args.unused:
		for x in unused:
			print('{}\tunused'.format(x), file=sys.stderr)
