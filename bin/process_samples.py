#!/opt/conda/bin/python3

import subprocess
import argparse


def main(inFile, fastq_dir, fast_test_mode):
	"""

	:param inFile:
	:param fastq_dir:
	:return:
	"""
	input_file = open(inFile)
	new_sample_sheet = inFile.split('.')[0] + '_new.csv'
	output_file = open(new_sample_sheet, 'w')


	for line in input_file:
		line = line.rstrip()
		line = line.split(',')
		if line[4].split('.')[-1] != 'fq' and line[4].split('.')[-1] != 'fastq' and line[4].split('.')[-1] != 'R1' and line[4].split('.')[-1] != 'gz':
			print('SRA detected')

			if fast_test_mode:
				print(line[5])
				subprocess.run(['fastq-dump', '-X', '5', '--split-files', line[5]])
			else:
				print(line[5])
				subprocess.run(['fastq-dump', '--split-files', line[5]])

			new_line = line[0] + ',' + line[1] + ',' + line[2] + ',' + line[3] + ',' + line[4] + ',' + fastq_dir + line[5] + '_1.fastq' + ',' + fastq_dir + line[5] + '_2.fastq' + '\n'
			print(new_line)
			output_file.write(new_line)

		else:
			new_line = ",".join(map(str, line))
			new_line = new_line + '\n'
			output_file.write(new_line)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='''This tool downloads any reads from the SRA in the sample sheet, and makes a new one with 
		the paths.
		''',  epilog="""Version 0.1""")

	parser.add_argument('-i', '--input_file', type=str, help='Input sample file')
	parser.add_argument('-f', '--fastq_dir', type=str, help='Directory where reads from the SRA will be moved')
	parser.add_argument('-q', '--quick_test_mode', type=bool, default=False, help='Only download first 5 spots from SRA')


	args = parser.parse_args()

	inFile = args.input_file
	fastq_dir = args.fastq_dir
	fast_test_mode = args.quick_test_mode

	main(inFile, fastq_dir, fast_test_mode)

