import argparse
import util

###############################################################################################
# Purpose: input adjacency list, paired mapped reads, kmer size, contig table and             #
# outputs variable edge coverage                                                               #
###############################################################################################

parser = argparse.ArgumentParser(description='write the fasta sequences of putative cycles')
parser.add_argument('sample', metavar='<minimum contig size to keep>', type=str, help='min contig size')
parser.add_argument('ifn_assembly', metavar='<minimum contig size to keep>', type=str, help='min contig size')
parser.add_argument('ifn_read_stats', metavar='<contig fasta input>', type=str, help='contig fasta input')
parser.add_argument('ofn_summary', metavar='<contig fasta input>', type=str, help='contig fasta input')

args = parser.parse_args()


num_reads = None
with open(args.ifn_read_stats) as reads:
    next(reads)
    read_stats = util.split(reads.readline())
    num_reads = sum([int(num) for num in read_stats])

contig_lengths = []
with open(args.ifn_assembly) as assembly:
    for line in assembly:
        if line[0] != ">":
            contig_lengths.append(len(line.rstrip()))

    ofn_sum = open(args.ofn_summary,"w+")
    total_size = sum(contig_lengths)
    cum_length = 0
    no_n50 = True
    no_n90 = True
    for index, length in enumerate(contig_lengths):
        cum_length += length
        if cum_length > total_size * 0.5 and no_n50:
            no_n50 = False
            N50 = length
            L50 = index + 1
        if cum_length > total_size * 0.9 and no_n90:
            no_n90 = False
            N90 = length
            L90 = index + 1
    ofn_sum.write(util.write_line("sample","num_reads","assembly_size","num_contigs","N50","L50","N90","L90"))
    ofn_sum.write(util.write_line(args.sample, num_reads, total_size, len(contig_lengths), N50, L50, N90, L90))
    

