import numpy
import csv
import itertools
import pandas as pd
import os
import shutil
import glob

from venn import venn
from matplotlib.pyplot import show, savefig

inpath = '/Volumes/External/CIDRI_Data/Mel_M_smeg/'

anno_vcf = inpath + 'SnpEff/'

ref = 'default'

sample_file = inpath + 'sample_sheet_test.csv'

comparison = ['Small_col', 'Large_col']

pd.set_option('display.max_colwidth', -1)

#comparison = '19119D-02-03' vs '19119D-02-04'

# Variant analysis
#X vs X samples
#What is common, what is unique - Venn
# for each group, do enrichment

# Phylogenetics


# 19119D-02-05
# 19119D-02-03 Small colony morphology isolate
# 19119D-02-04 Big colony morphology isolate


def input_parser(file_path):
	"""

	:param file_path:
	:return:
	"""

	if file_path[-4:] == ".csv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=',')
		data_matrix = list(data_table)
		result = numpy.array(data_matrix)
		return result

	if file_path[-4:] == ".txt":
		list_of_dicts = []
		reader = csv.DictReader(open(file_path, 'r'), delimiter='\t')
		for row in reader:
			list_of_dicts.append(row)
		return list_of_dicts

	if file_path[-4:] == ".vcf":
		list_of_dicts = []
		# Deal with random info at start
		in_file = open(file_path, 'r')
		entry_label = file_path
		for line in in_file:
			if line[0:2] != "##" and line[0] == "#":
				vcf_headder_line = line.split('\t')
				vcf_headder_line[0] = vcf_headder_line[0][1:]
				vcf_headder_line[-1] = vcf_headder_line[-1].strip()

			if not line.startswith('#'):
				entries = line.split('\t')
				entry_dict = {
					vcf_headder_line[0]: entries[0], vcf_headder_line[1]: entries[1],
					vcf_headder_line[2]: entries[2], vcf_headder_line[3]: entries[3],
					vcf_headder_line[4]: entries[4], vcf_headder_line[5]: entries[5],
					vcf_headder_line[6]: entries[6], vcf_headder_line[7]: entries[7],
					vcf_headder_line[8]: entries[8], vcf_headder_line[9]: entries[9].strip(),
					'ORIGIN': entry_label
				}

				list_of_dicts.append(entry_dict)
		return list_of_dicts


def venn_diagram(dict_of_data_lists):
	'''Creates the needed format for visualization of data as a venn diagram. Takes in a dictionary of lists as input, where the keys are the circle sey labels'''
	# print dict_of_data_lists

	list_of_sets = []

	# Get all possible combinations of input data keys

	list_of_labels = []
	list_of_label_conbinations = []
	list_of_values = []

	for label_string in dict_of_data_lists.keys():
		list_of_labels.append(label_string)
		list_of_values = list_of_values + dict_of_data_lists[label_string]

	# print list_of_labels
	list_of_values = unique_list(list_of_values)
	# print list_of_values

	count = 1

	while count <= len(list_of_labels):
		combinations = itertools.combinations(list_of_labels, count)

		for set_combo in combinations:
			set_combo = list(set_combo)
			list_of_label_conbinations.append(set_combo)
			list_of_sets.append({'sets': set_combo, 'size': 0, 'components': []})

		count += 1

	for set_value in list_of_values:
		set_value_set = []
		for a_set in dict_of_data_lists.keys():
			if set_value in dict_of_data_lists[a_set]:
				# print set_value , a_set
				set_value_set.append(a_set)
		# print set_value_set

		for set_dict in list_of_sets:
			if set_dict['sets'] == set_value_set:
				set_dict['size'] = set_dict['size'] + 1
				set_dict['components'].append(set_value)

	# print list_of_label_conbinations

	return list_of_sets


def unique_list(list):
	new_list = []
	for i in list:
		if i not in new_list:
			new_list.append(i)
	return new_list


def union(a, b):
	c = []
	d = []
	for i in a:
		if isinstance(i, str) is False:
			i = i.gene_ID
			c.append(i)
		else:
			c.append(i)
	for j in b:
		if isinstance(j, str) is False:
			j = j.gene_ID
			d.append(j)
		else:
			d.append(j)
	return list(set(c) | set(d))


def difference(a, b):
	c = []
	d = []
	for i in a:
		if isinstance(i, str) is False:
			i = i.gene_ID
			c.append(i)
		else:
			c.append(i)
	for j in b:
		if isinstance(j, str) is False:
			j = j.gene_ID
			d.append(j)
		else:
			d.append(j)
	return list(set(c) - set(d))


def intersect(a, b):
	# Check if string
	if isinstance(a, list):
		c = []
		d = []
		for i in a:
			if isinstance(i, str) is False:
				i = i.gene_ID
				c.append(i)
			else:
				c.append(i)
		for j in b:
			if isinstance(j, str) is False:
				j = j.gene_ID
				d.append(j)
			else:
				d.append(j)
		return list(set(c) & set(d))
	else:
		# Treat as a list of dicts
		return 'th'


class ReportObject:

	def __init__(self):
		self.title = ''
		self.body_elements = []
		self.style_list = [
			'<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css" />',
			'<link rel="stylesheet" href="assets/report_style.css" />',
		]
		self.style_string = ''
		self.body_html = '<body>'

	start_html = '<!DOCTYPE html> <html>'

	navbar_html = '''
	<div class="topnav"> 
	<a href="#">Variants</a>
	<a href="#">Virulence factors</a>
	<a href="#">RunQC</a>
	</div>
	'''

	footer_html = '</html>'

	def add_title(self, title_text):
		self.title = title_text

	def add_table(self, table_df, table_title):
		table_html = table_df.to_html(escape=False)
		table_html = '<div class="var_table"><h3>' + table_title + '</h3>' + table_html + '</div>'

		self.body_elements.append(table_html)

	def add_figure(self, figure_path, figure_title, figure_legend):

		figure_html = '<div class="figure">'

		figure_html += '<h3>' + figure_title + '</h3>'

		figure_html += '<img src="' + figure_path + '" alt="' + figure_title + '">'

		figure_html += '<p class="legend">' + figure_legend + '</p>'

		figure_html += '</div>'

		self.body_elements.append(figure_html)

	def write_html(self, filename):

		for a_style in self.style_list:
			self.style_string += a_style + '\n'

		header_html = '<head> <meta charset="utf-8" /> ' + self.style_string + '</head>'

		title_html = '<h1>' + self.title + '</h1>'

		self.body_html += title_html

		for an_element in self.body_elements:
			self.body_html = self.body_html + an_element
			self.body_html = self.body_html + '<br>'

		self.body_html += '</body>'

		out_file = open(filename, 'w')

		out_file.write(self.start_html)

		out_file.write(header_html)

		out_file.write(self.navbar_html)

		out_file.write(self.body_html)

		out_file.write(self.footer_html)


def make_dir_structure(working_dir, snpeff=True):
	path = working_dir + 'report'
	print('Creating report in: ' + path)
	try:
		os.mkdir(path)
		os.mkdir(path + '/variants')
		os.mkdir(path + '/assets')
	except OSError:
		print("Creation of the directory %s failed" % path)
	else:
		print("Successfully created the directory %s " % path)

	if snpeff:
		snpeff_dir = '/Volumes/External/CIDRI_Data/Mel_M_smeg/SnpEff'
		target_dir = path + '/variants/'

		for file in glob.glob(snpeff_dir + '/*.html'):
			shutil.copy(file, target_dir)

	# Get CSS code

	try:
		os.mkdir(path + '/assets')
	except OSError:
		print("Creation of the assets directory %s failed" % path)

	shutil.copy('../assets/report_style.css', './report/assets')


# Get sample file
reader = csv.DictReader(open(sample_file))
sample_info = {}
for line in reader:
	vcf_file_name = 'sample_' + line['number'] + '_sorted_dedup_filtered.recode_snpEff.ann.vcf'
	sample_info['sample_' + line['number']] = {
		'file': vcf_file_name,
		'origin': line['origin'],
		'replicate': line['replicate'],
		'isolate': line['isolate'],
		'phenotype': line['phenotype'],
		'vcf': input_parser(anno_vcf + vcf_file_name)
	}


# Venn diagram stuff
ven_pos_dict = {}

for a_sample in sample_info.keys():
	ven_pos_dict[a_sample] = []
	for a_variant_pos in sample_info[a_sample]['vcf']:
		ven_pos_dict[a_sample].append(a_variant_pos['POS'])


res = venn_diagram(ven_pos_dict)

# Create comparison

comparison_dict = {}
for a_comparison in comparison:
	comparison_dict[a_comparison] = []

for a_sample in sample_info.keys():
	if sample_info[a_sample]['phenotype'] in comparison:
		comparison_dict[sample_info[a_sample]['phenotype']].append(a_sample)


# Extract variants

# Any mutations that are found in only one phenotype, or only mutations found in all isolates of one phenotype


def extract_variants_by_phenotype(venn_data, phenotype_dict, mode='union'):
	"""
	:param venn_data:
	:param phenotype_dict:
	:param mode: String - union (All variants) or intersection (only shared)
	:return:
	"""
	result = {}

	for a_phenotype in phenotype_dict.keys():
		result[a_phenotype] = []

	for a_set in venn_data:

		if mode is 'union':

			for a_phenotype in phenotype_dict.keys():

				add_to_pheno = True

				for a_ven_sample in a_set['sets']:
					if a_ven_sample not in phenotype_dict[a_phenotype]:
						add_to_pheno = False

				if add_to_pheno:

					for variant in a_set['components']:

						result[a_phenotype].append(variant)

		if mode is 'intersection':

			for a_phenotype in phenotype_dict.keys():

				if set(phenotype_dict[a_phenotype]) == set(a_set['sets']):

					for variant in a_set['components']:

						result[a_phenotype].append(variant)

	return result


def get_var_summary(sample_info_dict, report_path):

	df_list = []

	for a_sample_info in sample_info_dict.keys():
		sample_list = [a_sample_info]
		print(sample_info_dict[a_sample_info]['file'])

		link_to_SnpEFF_html = 'variants/' + a_sample_info + "_sorted_dedup_filtered.recode_snpEff.html"
		html_link = '<a href="' + link_to_SnpEFF_html + '">' + a_sample_info + '</a>'

		sample_list.append(html_link)

		df_list.append(sample_list)

	a_df = pd.DataFrame.from_records(df_list)

	return a_df


phenotype_var_pos = extract_variants_by_phenotype(res, comparison_dict, mode='intersection')

# Get variant info

phenotype_var_detailed = {}

for a_phenotype in comparison_dict.keys():
	phenotype_var_detailed[a_phenotype] = []
	for a_sample in comparison_dict[a_phenotype]:
		for a_variant in sample_info[a_sample]['vcf']:
			if a_variant['POS'] in phenotype_var_pos[a_phenotype]:
				phenotype_var_detailed[a_phenotype].append(a_variant)



# Format for variant output

this = ReportObject()
this.add_title('M smegmatis report 1')


for pheno in phenotype_var_detailed:

	# properly extract the gene names

	for a_variant in phenotype_var_detailed[pheno]:
		# Annotation
		a_variant['Annotation'] = a_variant['INFO'].split('|')[1]

		# Putative_impact
		a_variant['Impact'] = a_variant['INFO'].split('|')[2]

		# Gene name
		a_variant['Gene_name'] = a_variant['INFO'].split('|')[3]

		# Gene ID
		a_variant['Gene_ID'] = a_variant['INFO'].split('|')[4]

	var_df = pd.DataFrame(phenotype_var_detailed[pheno])

	print(pheno)
	print(var_df)

	try:
		var_df = var_df[['POS', 'Gene_ID', 'Gene_name', 'REF', 'ALT', 'Annotation', 'Impact', 'QUAL', 'FORMAT', 'INFO', 'FILTER']]
	except KeyError:
		print('no variants')
	else:
		this.add_table(var_df, pheno + ' ' + str(comparison_dict[pheno]))

make_dir_structure('')

var_df = get_var_summary(sample_info, 'report/')
this.add_table(var_df, 'SnpEFF summaries')
this.add_figure('./variants/snp_venn_fig.png', 'Variants', 'Small_col includes: sample 19119D-02-03_S87. '
					'Large_col includes: sample 19119D-02-04_S88. '
					'Other includes variants found in samples: 19119D-02-01_S92, 19119D-02-02_S93, 19119D-02-05_S89, 19119D-02-06_S90'
				)

this.write_html('./report/report.html')


def make_snp_venn(venn_grouped_info, selected_comparison, save_path):

	experiment_venn_data = {'Other': set()}

	for phenotype in selected_comparison.keys():
		experiment_venn_data[phenotype] = set()

	comparison_sample_list = set()
	all_samples = set()

	for a_pheno_sample in selected_comparison.keys():
		for a_thing in selected_comparison[a_pheno_sample]:
			comparison_sample_list.add(a_thing)

	for venn_comparison in venn_grouped_info:
		for a_sample in venn_comparison['sets']:
			all_samples.add(a_sample)
		for phenotype in selected_comparison.keys():
			for a_pheno_sample in selected_comparison[phenotype]:
				if a_pheno_sample in venn_comparison['sets']:
					for var_pos in venn_comparison['components']:
						experiment_venn_data[phenotype].add(var_pos)

	other_samples = (all_samples - comparison_sample_list)

	for venn_comparison in venn_grouped_info:
		for other_sample in other_samples:
			if other_sample in venn_comparison['sets']:
				for var_pos in venn_comparison['components']:
					experiment_venn_data['Other'].add(var_pos)

	venn(experiment_venn_data)

	savefig(save_path + 'variants/snp_venn_fig.png', format="png")


make_snp_venn(res, comparison_dict, 'report/')


def link_to_pathways(list_of_genes, organism):

	return 'A graph obj'


def GSEA_on_gene_list(list_of_genes, organism):

	return 'table of some sort'



