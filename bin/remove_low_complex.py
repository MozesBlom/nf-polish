#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
	import os
	import sys
	import gzip
	import argparse
except ImportError:
	sys.exit("One of the required modules can't be found...")


#########################
## Input - parse arguments
#########################

parser = argparse.ArgumentParser()

parser.add_argument("-1", "--forward_reads", help="Path to forward reads in read pair", action = "store")
parser.add_argument("-2", "--reverse_reads", help="Path to reverse reads in read pair", action = "store")
parser.add_argument("-u", "--unpaired_reads", help="Path orphan reads", action = "store")
parser.add_argument("-c", "--cut_off", help="Cut-off value to be recognizes as low-complexity", action = "store")
parser.add_argument("-p", "--prefix", help="Prefix ID to be used for output fastq files", action = "store")


args = parser.parse_args()

freads_path = args.forward_reads
rreads_path = args.reverse_reads
ureads_path = args.unpaired_reads
lc_val = args.cut_off
prefix_id = args.prefix


#########################
## Function
#########################

def flagLowComplexity(fq_gz_file, cut_off):
	read_set = set()
	with gzip.open(fq_gz_file, 'rt') as fq_reads:
		while 1:
			read_name = fq_reads.readline().strip().split(' ')[0]
			if not read_name: break
			seq = fq_reads.readline().strip()
			fq_reads.readline()
			fq_reads.readline()
			seq_len = float(len(seq))
			for base in ['A','C','T','G','N']:
				if float(seq.count(base))/seq_len >= float(cut_off):
					read_set.add(read_name)
	return read_set

#########################
## Automated Analysis
#########################

low_complex_1 = flagLowComplexity(freads_path, lc_val)
low_complex_2 = flagLowComplexity(rreads_path, lc_val)
low_complex_u = flagLowComplexity(ureads_path, lc_val)
lowcomp_both = low_complex_1 & low_complex_2
lowcomp_1_only = low_complex_1 - low_complex_2
lowcomp_2_only = low_complex_2 - low_complex_1
# New FastQ files
out1_P = open('%s_R1.fastq' % prefix_id, 'w')
out2_P = open('%s_R2.fastq' % prefix_id, 'w')
out_U_new = open('%s_U.fastq' % prefix_id, 'w')
out_low_complex = open('%s_low_complex.fastq' % prefix_id, 'w')
# Open Old FastQ files
in1_P = gzip.open(freads_path, 'rt')
in2_P = gzip.open(rreads_path, 'rt')
inU = gzip.open(ureads_path, 'rt')
# Now go over the old file and remove low duplicate sequences
while 1:
	# Store read data
	nameline_1 = in1_P.readline()
	seq_1 = in1_P.readline()
	del_1 = in1_P.readline()
	qual_1 = in1_P.readline()
	nameline_2 = in2_P.readline()
	seq_2 = in2_P.readline()
	del_2 = in2_P.readline()
	qual_2 = in2_P.readline()
	if not nameline_1: break
	rname_1 = nameline_1.strip().split(' ')[0]
	rname_2 = nameline_2.strip().split(' ')[0]
	# Compare to list
	if rname_1 in lowcomp_both:
		#continue
		out_low_complex.write('%s\n%s%s%s' % (rname_1, seq_1, del_1, qual_1))
		out_low_complex.write('%s\n%s%s%s' % (rname_2, seq_2, del_2, qual_2))
	elif rname_1 in lowcomp_1_only:
		out_U_new.write('%s\n%s%s%s' % (rname_2, seq_2, del_2, qual_2))
		out_low_complex.write('%s\n%s%s%s' % (rname_1, seq_1, del_1, qual_1))
	elif rname_2 in lowcomp_2_only:
		out_U_new.write('%s\n%s%s%s' % (rname_1, seq_1, del_1, qual_1))
		out_low_complex.write('%s\n%s%s%s' % (rname_2, seq_2, del_2, qual_2))
	else:
		out1_P.write('%s\n%s%s%s' % (rname_1, seq_1, del_1, qual_1))
		out2_P.write('%s\n%s%s%s' % (rname_2, seq_2, del_2, qual_2))
while 1:
	nameline_u = inU.readline()
	seq_u = inU.readline()
	del_u = inU.readline()
	qual_u = inU.readline()
	if not nameline_u: break
	rname_u = nameline_u.strip().split(' ')[0]
	if rname_u in low_complex_u:
		out_low_complex.write('%s\n%s%s%s' % (rname_u, seq_u, del_u, qual_u))
	else:
		out_U_new.write('%s\n%s%s%s' % (rname_u, seq_u, del_u, qual_u))
in1_P.close()
in2_P.close()
inU.close()
out1_P.close()
out2_P.close()
out_U_new.close()
out_low_complex.close()



