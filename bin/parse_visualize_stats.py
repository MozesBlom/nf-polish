#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
	import os
	import sys
	import argparse
	import matplotlib
	matplotlib.use('Agg')
	import pandas as pd
	import seaborn as sns
	import matplotlib.pyplot as plt
except ImportError:
	sys.exit("One of the required modules can't be found...")


#########################
## Input - parse arguments
#########################

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fastqc_stats", help="File with reads stats: raw reads", default=False, action = "store", required=True)
parser.add_argument("-d", "--dedup_stats", help="File with reads stats: post deduplication", default=False, action = "store", required=False)
parser.add_argument("-a", "--adapt_stats", help="File with reads stats: post adapter trimming", default=False, action = "store", required=True)
parser.add_argument("-m", "--merge_stats", help="File with reads stats: post merging", default=False, action = "store", required=True)
parser.add_argument("-q", "--qual_stats", help="File with reads stats: post quality trimming", default=False, action = "store", required=True)
parser.add_argument("-c", "--complex_stats", help="File with reads stats: post low complexity removal", default=False, action = "store", required=False)
parser.add_argument("-p", "--prefix_string", help="Prefix to be used for output fn; often library/indiv name", default=False, action = "store", required=True)
parser.add_argument("-o", "--out_dir", help="Directory to store the output figure, if unspecified assumes current work dir", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
fastqc_fn = args.fastqc_stats
dedup_fn = args.dedup_stats
adapt_fn = args.adapt_stats
merge_fn = args.merge_stats
qual_fn = args.qual_stats
complex_fn = args.complex_stats
prefix = args.prefix_string
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()


#########################
## Function
#########################

# A function to parse the seqkit file and search for file stats based on read type
# File type is a list containing all possible strings that may be in the file name
def parse_seqkit(fn_in, file_type):
	df = pd.read_csv(fn_in, delim_whitespace = True)
	for index, row in df.iterrows():
		for x in file_type:
			if x in row["file"]:
				read_number = row["num_seqs"]
				read_length = row["avg_len"]
		else:
			pass
	return read_number, read_length

#########################
## Automated Analysis
#########################

## Create a pandas dataframe with a summary of all the read stats
df = pd.DataFrame(columns=['Process','Type','Number','Length'])

if os.path.isfile(fastqc_fn):
	F_number, F_length = parse_seqkit(fastqc_fn, ['_R1_', '_R1.', '.R1_', '.R1.'])
	R_number, R_length = parse_seqkit(fastqc_fn, ['_R2_', '_R2.', '.R2_', '.R2.'])
	U_number, U_length = 0, 0
	fastqc_df = pd.DataFrame([['1.raw_reads', 'Forward', F_number, F_length], ['1.raw_reads', 'Reverse', R_number, R_length], ['1.raw_reads', 'Unpaired', U_number, U_length]], columns=['Process','Type','Number','Length'])
	df = pd.concat([df, fastqc_df], ignore_index=True)
else:
	sys.exit("Can't find the read stats file for raw reads! Please check the -f flag")

if (dedup_fn != False) & (os.path.isfile(dedup_fn)):
	F_number, F_length = parse_seqkit(dedup_fn, ['_R1_', '_R1.', '.R1_', '.R1.'])
	R_number, R_length = parse_seqkit(dedup_fn, ['_R2_', '_R2.', '.R2_', '.R2.'])
	U_number, U_length = '0', '0'
	dedup_df = pd.DataFrame([['2.deduplication', 'Forward', F_number, F_length], ['2.deduplication', 'Reverse', R_number, R_length], ['2.deduplication', 'Unpaired', U_number, U_length]], columns=['Process','Type','Number','Length'])
	df = pd.concat([df, dedup_df], ignore_index=True)
elif (dedup_fn == False):
	pass
else:
	sys.exit("Can't find the read stats file for dedup reads! Please check the -d flag")

if os.path.isfile(adapt_fn):
	if (dedup_fn != False):
		step = '3.adapter_trimming'
	else:
		step = '2.adapter_trimming'
	F_number, F_length = parse_seqkit(adapt_fn, ['_R1_', '_R1.', '.R1_', '.R1.'])
	R_number, R_length = parse_seqkit(adapt_fn, ['_R2_', '_R2.', '.R2_', '.R2.'])
	U_number, U_length = parse_seqkit(adapt_fn, ['_U_', '_U.', '.U_', '.U.'])
	adapt_df = pd.DataFrame([[step, 'Forward', F_number, F_length], [step, 'Reverse', R_number, R_length], [step, 'Unpaired', U_number, U_length]], columns=['Process','Type','Number','Length'])
	df = pd.concat([df, adapt_df], ignore_index=True)
else:
	sys.exit("Can't find the read stats file for adapter trimmed reads! Please check the -a flag")

if os.path.isfile(merge_fn):
	if (dedup_fn != False):
		step = '4.merging'
	else:
		step = '3.merging'
	F_number, F_length = parse_seqkit(merge_fn, ['forward'])
	R_number, R_length = parse_seqkit(merge_fn, ['reverse'])
	U_number, U_length = parse_seqkit(merge_fn, ['_U_', '_U.', '.U_', '.U.'])
	merge_df = pd.DataFrame([[step, 'Forward', F_number, F_length], [step, 'Reverse', R_number, R_length], [step, 'Unpaired', U_number, U_length]], columns=['Process','Type','Number','Length'])
	df = pd.concat([df, merge_df], ignore_index=True)
else:
	sys.exit("Can't find the read stats file for adapter trimmed reads! Please check the -a flag")

if os.path.isfile(qual_fn):
	if (dedup_fn != False):
		step = '5.quality_trimming'
	else:
		step = '4.quality_trimming'
	F_number, F_length = parse_seqkit(qual_fn, ['_R1_', '_R1.', '.R1_', '.R1.'])
	R_number, R_length = parse_seqkit(qual_fn, ['_R2_', '_R2.', '.R2_', '.R2.'])
	U_number, U_length = parse_seqkit(qual_fn, ['_U_', '_U.', '.U_', '.U.'])
	qual_df = pd.DataFrame([[step, 'Forward', F_number, F_length], [step, 'Reverse', R_number, R_length], [step, 'Unpaired', U_number, U_length]], columns=['Process','Type','Number','Length'])
	df = pd.concat([df, qual_df], ignore_index=True)
else:
	sys.exit("Can't find the read stats file for adapter trimmed reads! Please check the -a flag")

if os.path.isfile(complex_fn):
	if (dedup_fn != False):
		step = '6.low_complexity'
	else:
		step = '5.low_complexity'
	F_number, F_length = parse_seqkit(complex_fn, ['_R1_', '_R1.', '.R1_', '.R1.'])
	R_number, R_length = parse_seqkit(complex_fn, ['_R2_', '_R2.', '.R2_', '.R2.'])
	U_number, U_length = parse_seqkit(complex_fn, ['_U_', '_U.', '.U_', '.U.'])
	complex_df = pd.DataFrame([[step, 'Forward', F_number, F_length], [step, 'Reverse', R_number, R_length], [step, 'Unpaired', U_number, U_length]], columns=['Process','Type','Number','Length'])
	df = pd.concat([df, complex_df], ignore_index=True)
else:
	sys.exit("Can't find the read stats file for adapter trimmed reads! Please check the -a flag")


## Store the read stat results in a single file
tsv_fn = os.path.join(output_path, (prefix + '_polishing_table.tsv'))
df.to_csv(tsv_fn, sep='\t', index=False)

## Make sure that the numbers can be interpreted correctly by seaborn
df = df.replace(',','', regex=True)
df = df.astype({'Number': 'int64', 'Length':'float'})


## Now let's plot how the read number and length changed during polishing
## Specify the figure aesthetics
sns.set()
sns.set(font_scale=1, font='sans-serif')
sns.set_style("whitegrid", {'axes.grid' : False})

## Subplots: One for read number, the other for read length

# Create a stacked figure
fig, ax = plt.subplots(nrows=(2), sharex=True, figsize=(8,8))
fig.subplots_adjust(wspace=0.2)

# Create two cat plots -- Note: it is essential to include plt.close when using seaborn cat plot in subplots
sns_plot1 = sns.catplot(x="Process", y="Number", hue="Type",
                        capsize=.2, palette="YlGnBu_d",
                        kind="point", ax=ax[0], data=df)
plt.close()

sns_plot2 = sns.catplot(x="Process", y="Length", hue="Type",
                        capsize=.2, palette="YlGnBu_d",
                        kind="point", ax=ax[1], data=df)
plt.close()

## Improve figure
plt.xticks(rotation=45, fontweight='bold')
ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)
ax[0].set_ylabel('Read number', fontweight="bold", fontsize=14, labelpad=20)
ax[0].set_xlabel('')
ax[0].set_title("Polishing overview for library %s" % (prefix), fontweight="bold", fontsize=24, pad=40)
ax[0].get_legend().remove()
ax[1].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].set_ylabel('Read length', fontweight="bold", fontsize=14, labelpad=20)
ax[1].set_xlabel('Clean up step', fontweight="bold", fontsize=14, labelpad=20)
plt.legend(bbox_to_anchor=(1.30,1.38), prop={'size': 14}, frameon=False)

## Store output
fig_fn = os.path.join(output_path, (prefix + '_polishing_overview.png'))
plt.savefig(fig_fn, bbox_inches = "tight", dpi=800)


