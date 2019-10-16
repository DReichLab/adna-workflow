import argparse
import sys

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Print the desired number of lines from a file")
	parser.add_argument('-n', "--num_lines", help="number of lines to output", default=10, type=int)
	parser.add_argument('-s', "--skip", help="number of initial lines to output", default=0, type=int)
	parser.add_argument('inputfile', nargs='?', type=argparse.FileType('r'))
	args = parser.parse_args() 
	
	#with open(args.inputfile) as f:
	for n in range(args.skip):
		args.inputfile.readline()
	for n in range(args.num_lines):
		print(args.inputfile.readline(), end='')
