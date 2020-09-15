import argparse
import util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Get contig coverages from parsed bwa map file')
parser.add_argument('ifn_paired_table', metavar='<contig table>', type=str, help='contig table w/ lengths')
parser.add_argument('read_stats', metavar='<contig table>', type=str, help='contig table w/ lengths',)
parser.add_argument('graph', metavar='<contig table>', type=str, help='contig table w/ lengths', default=None)
parser.add_argument('read', metavar='<fastq read file>', type=str, help='contig table w/ lengths',default=None)
parser.add_argument('assembly', metavar='<assembly file>', type=str, help='contig table w/ lengths',default=None)
parser.add_argument('paired', metavar='<paired boolean>', type=str, help='contig table w/ lengths')

args = parser.parse_args()

if args.paired == "True":
    read_length = util.get_read_length(args.read_stats)
    insert_dist = util.get_insert(args.ifn_paired_table)
    k = util.get_k(args.read_stats)
    insert_dist = sorted(insert_dist)
    max_distance = np.percentile(insert_dist, 99.5) + 200
    med_insert = np.median(insert_dist)
    hist = plt.hist(insert_dist, bins=1000)
    plt.title('Read Length Distribution')
    plt.xlabel('Sequence Molecule Size')
    plt.ylabel('Count')
    plt.xlim(left=0,right=1250)
    ofn_graph = args.graph
    plt.savefig(ofn_graph)
    plt.close()
    ofile = open(args.read_stats,"w+")
    ofile.write(util.write_line("read_length", "k", "med_insert", "max_distance"))
    ofile.write(util.write_line(read_length, k, med_insert, max_distance))

else:
    ofile = open(args.read_stats, "w+")
    with open(args.read) as read:
        count = 0
        lengths = []
        for line in read:
            if count % 4 == 1:
                lengths.append(len(line.rstrip()))
            count += 1
            if count == 1000:
                break
        ofile.write(str(int(sum(lengths) / len(lengths))) + "\n")
    with open(args.assembly) as assembly:
        assembly_header = assembly.readline().split(" ")[0]
        k = assembly_header.split("_")[0][2:]
        ofile.write(k + "\n")