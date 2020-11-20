import pysam
import os
import sys
import re
import optparse


def PolyT(seq, pattern):
	"""
	input sequence infortmation, and return the length of polyT
	"""
	dT = re.findall(pattern, seq)

	return len(dT[0])


def main():
	usage = """munibam.py -i <Bam> -o <outBam> -b <int>
			"""
	fmt = optparse.IndentedHelpFormatter(max_help_position=50, width=100)
	parser = optparse.OptionParser(usage=usage, formatter=fmt)
	group  = optparse.OptionGroup(parser, 'Query arguments', 
								'These options define search query arguments and parameters.')
	group.add_option("-i", "--in", action = "store", dest = "filein", type = str,
					help = "input bam file")
	group.add_option("-o", "--out", action="store", dest="fileout", type=str,
					help = "output bam file")
	group.add_option("-b", "--barcode", action = "store", dest = 'barcode',type = int,
					help = "The last positions of the barcode and UMI. -b 26[10X], 54[Microwell-seq]")
	parser.add_option_group(group)
	
	options, argvs = parser.parse_args()
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	input_hd = pysam.AlignmentFile(options.filein, 'rb', check_sq=False, )
	outbamfile = pysam.AlignmentFile(options.fileout, 'wb', template = input_hd)

	flags = 77

	count_ind = 0
	MAX_report = 1000000
	pattern = re.compile('^T*')

	for line in input_hd.fetch(until_eof=True):

		count_ind += 1
		if count_ind % MAX_report == 0:
			sys.stdout.flush()
			print(f"Parsed {count_ind} reads.")

		if line.flag == flags:
			dTlen = PolyT(line.query_sequence[options.barcode:], pattern)
			p = line.query_qualities

			add_tags = [('DT', line.query_sequence[options.barcode:options.barcode + dTlen])]

			line.tags += add_tags
		else:
			line.tags += add_tags

		outbamfile.write(line)

	samfile.close()
	outbamfile.close()

if __name__ == '__main__':
	main()
