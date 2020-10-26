from os import listdir
from os.path import isfile, join
import argparse

'''
This can be used when there are a large number of files that need to be processed in a directory
'''


def main(sample_dir_path, outfile_dir, outfile_name, sample_origin):

	#sample_dir_path = './'

	onlyfiles = [f for f in listdir(sample_dir_path) if isfile(join(sample_dir_path, f))]

	#outfile_name = 'sample_sheet.csv'
	#outfile_dir = './'
	#sample_origin = 'genomic'
	sample_dir = sample_dir_path

	out_file = open(outfile_dir + outfile_name, 'w')

	allowed_formats = ['gz', 'fastq', 'fq']

	R2_options = ['2', 'R2', 'r2']

	default_pheno = 'NA'

	replicate = '1'

	out_file.write('number,origin,replicate,isolate,phenotype,R1,R2\n')

	sample_count = 0

	for file in onlyfiles:

		if file.split('.')[-1].lower() in allowed_formats:

			seq_info = file.split('.')[0]
			seq_pair = seq_info.split('_')[-1]
			seq_id = '_'.join(seq_info.split('_')[:-1])
			extension = '.'.join(file.split('.')[1:])

			if seq_pair not in R2_options:
				# Use it.
				out_line_list = list()
				out_line_list.append(str(sample_count))
				out_line_list.append(sample_origin)
				out_line_list.append(replicate)
				out_line_list.append(seq_id)
				out_line_list.append(default_pheno)
				out_line_list.append(sample_dir + file)
				out_line_list.append(sample_dir + seq_id + '_' + seq_pair.replace('1', '2') + '.' + extension)
				out_string = ','.join(out_line_list)

				out_file.write(out_string + '\n')

				sample_count += 1

	out_file.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Create sample sheet from a directory of files')

	# TODO: Allow files to be added to an existing sample file
	# TODO: some phenotype sample mapping

	parser.add_argument('-s', '--sample_dir', type=str, required=True, help='Full path, not relative')
	parser.add_argument('-of', '--outfile_dir', type=str, required=True, help='info')
	parser.add_argument('-on', '--outfile_name', type=str, help='info', default='sample_sheet.csv')
	parser.add_argument('-jp', '--sample_origin', type=str, help='info', default='genomic')


	args = parser.parse_args()

	sample_dir_path = args.sample_dir
	outfile_dir = args.outfile_dir
	outfile_name = args.outfile_name
	sample_origin = args.sample_origin

	main(sample_dir_path, outfile_dir, outfile_name, sample_origin)
