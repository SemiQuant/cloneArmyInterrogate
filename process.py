#!/bin/python
import sys
import os
import pysam
from operator import attrgetter
from collections import defaultdict
import pandas as pd
from collections import Counter


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


# base_qual=20
# file_in="/Users/semiquant/Downloads/tmp/Sample-2.bam"
# ref_len = samfile.get_reference_length("rv0678")
# fasta = pysam.FastaFile("/Users/semiquant/Bioinformatics/Projects/BDQmut/refs/ref.fasta")
# fasta = fasta.fetch("rv0678")

base_qual = sys.argv[1]
file_in = sys.argv[2]
ref_name = sys.argv[3]
fasta = sys.argv[4]


samfile = pysam.AlignmentFile(file_in, "rb")
out_seqs = list()
total_reads = 0
mut = list()
ref_len = samfile.get_reference_length(ref_name)
fasta = fasta.fetch(ref_name)

for read1, read2 in read_pair_generator(samfile):
    read1_pos = read1.get_reference_positions()
    read1_qual = read1.query_qualities
    read1_seq = read1.get_reference_sequence()
    read2_pos = read2.get_reference_positions()
    read2_qual = read2.query_qualities
    read2_seq = read2.get_reference_sequence()
    pos = 0
    tmp_read = str()
    for r in fasta:
        #if pos in either read, add, otherwise add in the reference sequence
        if pos in read1_pos and read1_qual[read1_pos.index(pos)] > base_qual:
            tmp_read += read1_seq[read1_pos.index(pos)]
        elif pos in read2_pos and read2_qual[read2_pos.index(pos)] > base_qual:
            tmp_read += read2_seq[read2_pos.index(pos)]
        else:
            tmp_read += r
        pos += 1
    out_seqs.append(tmp_read)
    total_reads += 1

# open file in write mode
with open(file_in + ".seqs.txt", 'w') as f:
    for item in out_seqs:
        # write each item on a new line
        f.write("%s\n" % item)
    f.close()

out_seqs_up = [x.upper() for x in out_seqs]
res = Counter(out_seqs_up)

out_seqs_up = [x.islower() for x in out_seqs]
res = Counter(out_seqs_up)


for i in out_seqs:
    mut.append(sum([x>='a'and x<='z' for x in i]))
        
tmp_df=pd.Series(mut).value_counts()


for i in out_seqs:
    mut.append(sum([x>='a'and x<='z' for x in i]))
        
with open(file_in + ".muts.txt", 'w') as f:
    for s in mut:
        f.write(str(s) +"\n")
    f.close()

