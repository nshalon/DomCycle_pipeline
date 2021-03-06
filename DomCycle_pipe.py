import argparse
import os
import subprocess
import datetime


def run_steps():
    parser = argparse.ArgumentParser(description='Run DomCycle')
    parser.add_argument('-g', "--graph", type=str, help='Path to FASTG file created from assembly')
    parser.add_argument('-1', "--r1", metavar='R1', type=str, help='read1.fastq')
    parser.add_argument('-2', "--r2", metavar='R2', type=str, help='read2.fastq')
    parser.add_argument('-k', "--k", metavar='KMER SIZE', type=int, help='kmer size used for assembly')
    parser.add_argument("--mega", metavar='MEGAHIT DIRECTORY', type=str, help='the full path to megahit')
    parser.add_argument("--bwa", metavar='BWA DIRECTORY', type=str, help='the full path to bwa')
    parser.add_argument("--mega_tools", metavar='MEGAHIT_TOOLKIT', type=str, help='the full path to megahit toolkits')
    # parser.add_argument('-r', "--read_len", metavar='AVG (PRE-TRIMMED) READ LENGTH', type=float, help='input the average read length in the input library', default=150)
    parser.add_argument('-o', "--outdir", metavar='OUTPUT_DIR', type=str, help='Parent output directory', default="/".join([os.getcwd(),"output"]))
    parser.add_argument('-a', "--alpha", metavar='OUT COVERAGE DAMPENENING FACTOR', type=float, help='Out coverage dampening factor for finding cycles in graph', default=0)
    parser.add_argument('-p', "--maxpval", metavar='P-VALUE', type=float, help='Minimum p-value used to classify dominant cycles', default=0.01)
    parser.add_argument("--minscore", metavar='MIN SCORE', type=float, help='Minimum score for a dominant cycle', default=1)
    parser.add_argument('--sample', metavar='SAMPLE NAME', type=str, help='your custom sample name', default="metagenome")
    parser.add_argument('--min_quality', metavar='MIN QUALITY THRESHOLD', type=int, help='min map quality threshold', default=0)
    parser.add_argument('--min_match_len', metavar='MIN MATCH READ LENGTH', type=int, help='minimum match read length for filtering mapped reads', default=50)
    parser.add_argument('--max_edit_distance', metavar='MAX MISMATCH', type=int, help='maximum mismatches tolerated for filtering mapped reads', default=1)
    parser.add_argument('--steps', metavar='STEPS TO RUN', type=str, help='maximum mismatches tolerated for filtering mapped reads', default="1-16")
    parser.add_argument('--plot', metavar='PLOT (T/F)', type=bool, help='run plots or not', default=True)

    args = parser.parse_args()

    steps = parse_steps(args.steps)
    odir = args.outdir
    py_dir = os.path.join(__file__, "md/Python")
    pl_dir = os.path.join(__file__, "md/pl")
    r_dir = os.path.join(__file__, "md/R")

    if not os.path.exists(odir):
        os.mkdir(odir)
    else:
        print("WARNING!\nOutput directory already exists - consider changing the directory to avoid overwriting any previous results.")

    if args.outdir[0] != "/":
        args.outdir = os.path.join(os.getcwd(), args.outdir)
    if args.graph[0] != "/":
        args.graph = os.path.join(os.getcwd(), args.graph)
    if args.r1[0] != "/":
        args.r1 = os.path.join(os.getcwd(), args.r1)
    if args.r2[0] != "/":
        args.r2 = os.path.join(os.getcwd(), args.r2)

    if 1 in steps:
        print("Running the assembler...")
        subprocess.call(["megahit", "-m", "0.84", "-o", os.path.join(odir,"megahit"), "--min-contig-len", str(3*int(args.k)),
                         "--k-min", "27", "--k-max", str(args.k), "--k-step", "10", "--merge-level", "1000,0.95", "-1", args.r1,
                         "-2", args.r2, "-t", "40"],
                        cwd=args.mega)

    if 2 in steps:
        print("Getting read stats...")
        subprocess.call(["python", "get_read_stats.py", os.path.join(odir, "filter_paired_table"), os.path.join(odir, "read_stats.txt"),
                         os.path.join(odir, "read_distribution.pdf"), args.r1, os.path.join(odir, "megahit", "final.contigs.fa"),
                         "False"], cwd=py_dir)

    if 3 in steps:
        print("Creating FASTG...")
        subprocess.run(["megahit_toolkit", "contig2fastg", str(args.k), os.path.join(odir, "megahit", "final.contigs.fa")],
                         cwd=args.mega, stdout=os.path.join(odir, "k" + str(args.k) + ".fastg"))

    if 4 in steps:
        print("\n\nParsing assembly graph...")
        subprocess.call(["python", "parse_fastg.py", args.graph, os.path.join(odir, "adjacency_list"),
                         os.path.join(odir, "contig_rename_map"), os.path.join(odir, "renamed_final_contigs.fa"),
                         os.path.join(odir, "adjacency_matrix")], cwd=py_dir)

    if 5 in steps:
        print("Indexing assembly...")
        subprocess.call(["bwa", "index", os.path.join(odir, "renamed_final_contigs.fa")], cwd=args.bwa)

    if 6 in steps:
        print("Trimming reads...")
        subprocess.call(["perl", "trim_fastq.pl", args.r1, "0", "50", os.path.join(odir, "trimR1.fastq")], cwd=pl_dir)
        subprocess.call(["perl", "trim_fastq.pl", args.r2, "0", "50", os.path.join(odir, "trimR2.fastq")], cwd=pl_dir)

    if 7 in steps:
        print("Mapping reads...")
        subprocess.run(["bwa", "mem", "-t", "40", os.path.join(odir, "renamed_final_contigs.fa"), os.path.join(odir, "trimR1.fastq")],
                       cwd=args.bwa, stdout=os.path.join(odir, "R1_map.sam"))
        subprocess.run(["bwa", "mem", "-t", "40", os.path.join(odir, "renamed_final_contigs.fa"), os.path.join(odir, "trimR2.fastq")],
            cwd=args.bwa, stdout=os.path.join(odir, "R2_map.sam"))

    if 8 in steps:
        print("\n\nParsing R1...")
        subprocess.call(["perl", "parse_bwa_sam.pl", args.r1, os.path.join(odir, "R1table"),
                         os.path.join(odir, ".R1table_stats")], cwd=pl_dir)
        print("Parsing R2...")
        subprocess.call(["perl", "parse_bwa_sam.pl", args.r2, os.path.join(odir, "R2table"),
                         os.path.join(odir, ".R2table_stats")], cwd=pl_dir)
    if 9 in steps:
        print("\n\nFiltering reads: min quality", args.min_quality, "; min match length", args.min_match_len, "; max edit distance", args.max_edit_distance)
        print("Filtering R1...")
        subprocess.call(["perl", "filter_map.pl", os.path.join(odir, "R1table"), str(args.min_quality),
                         str(args.min_match_len), str(args.max_edit_distance), os.path.join(odir, "filter_R1table"),
                         os.path.join(odir, ".filter_R1table_stats")], cwd=pl_dir)
        print("Filtering R2...")
        subprocess.call(["perl", "filter_map.pl", os.path.join(odir, "R2table"), str(args.min_quality),
                         str(args.min_match_len), str(args.max_edit_distance), os.path.join(odir, "filter_R2table"),
                         os.path.join(odir, ".filter_R2table_stats")], cwd=pl_dir)
        subprocess.call(["python", "assembly_stats.py", args.sample, os.path.join(odir, "renamed_final_contigs.fa"),
                         os.path.join(odir, ".filter_R1table_stats"), os.path.join(odir, "assembly_stats")], cwd=py_dir)

    print("\n\nBUILDING GRAPH...")

    if 10 in steps:
        print("\n\nCalculating contig (internal edge) coverages...")
        subprocess.call(["python", "get_contig_cov.py", os.path.join(odir, "renamed_final_contigs.fa"),
                         os.path.join(odir, "filter_R1table"), os.path.join(odir, "filter_R2table"),
                         os.path.join(odir, "read_stats.txt"), os.path.join(odir, "contig_table")], cwd=py_dir)

    if 11 in steps:
        print("\n\nPairing reads and calculating read statistics...")
        subprocess.call(["python", "pair_reads.py", os.path.join(odir, "filter_R1table"),
                         os.path.join(odir, "filter_R2table"), os.path.join(odir, "filter_singleton_table"),
                         os.path.join(odir, "filter_paired_table")], cwd=py_dir)
        subprocess.call(["python", "get_read_stats.py", os.path.join(odir, "filter_paired_table"),
                         os.path.join(odir, "read_stats.txt"),
                         os.path.join(odir, "read_distribution.pdf"), args.r1,
                         os.path.join(odir, "megahit", "final.contigs.fa"),
                         "True"], cwd=py_dir)

    if 12 in steps:
        print("\n\nFinding all external edges...")
        subprocess.call(["python", "find_new_variable_edges.py", os.path.join(odir, "filter_paired_table"),
                         os.path.join(odir, "adjacency_list"), os.path.join(odir, "contig_table"),
                         os.path.join(odir, "read_stats.txt"), os.path.join(odir, "edge_summary")], cwd=py_dir)

    print("\n\nGRAPH BUILT!")

    if 13 in steps:
        print("\n\nFinding cycles in the graph...")
        subprocess.call(["python", "dominant_cycles.py", os.path.join(odir, "edge_summary"),
                         os.path.join(odir, "read_stats.txt"), str(args.alpha),
                         os.path.join(odir, "cycle_contig_table")], cwd=py_dir)

    if 14 in steps:
        print("\n\nCalculating cycle coverage statistics in the cycle space...")
        subprocess.call(["python", "cycle_coverages_2.py", args.sample, os.path.join(odir, "cycle_contig_table"),
                         os.path.join(odir, "filter_paired_table"), os.path.join(odir, "filter_singleton_table"),
                         os.path.join(odir, "contig_table"), os.path.join(odir, "read_stats.txt"),
                         os.path.join(odir, "cycle_covs_long"), os.path.join(odir, "cycle_cov_summary"),
                         os.path.join(odir, "out_contigs")], cwd=py_dir)

    if 15 in steps:
        print("\n\nCreating cycle fastas...")
        subprocess.call(["python", "cycle_fastas.py", os.path.join(odir, "cycle_contig_table"),
                         os.path.join(odir, "renamed_final_contigs.fa"), os.path.join(odir, "read_stats.txt"),
                         os.path.join(odir, "cycles.fasta")], cwd=py_dir)

    if 16 in steps:
        if not os.path.exists(os.path.join(odir, "dominant_cycles")):
            os.mkdir(os.path.join(odir, "dominant_cycles"))
        print("\n\nIdentifying Dominant Cycles...")
        subprocess.call(["python", "extract_p.py", os.path.join(odir, "cycle_cov_summary"),
                         os.path.join(odir, "cycles.fasta"), os.path.join(odir, "cycle_contig_table"),
                         str(args.maxpval), "0", "0", str(args.minscore), os.path.join(odir, "cycle_stats"), os.path.join(odir, "dominant_cycles", "cycle_stats"),
                         os.path.join(odir, "dominant_cycles", "cycles.fasta"), os.path.join(odir, "dominant_cycles", "cycle_contig_table")],
                         cwd=py_dir)

    if args.plot:
        if not os.path.exists(os.path.join(odir, "cycle_covs")):
            os.mkdir(os.path.join(odir, "cycle_covs"))
        subprocess.call(["Rscript", "plot_cycle_cov.R", "-i", os.path.join(odir, "cycle_covs_long"),
                         "-p", os.path.join(odir, "dominant_cycles", "cycle_contig_table"),
                         "-s", os.path.join(odir, "dominant_cycles", "cycle_stats"),
                         "-o", os.path.join(odir, "cycle_covs")], cwd=r_dir)


def parse_steps(steps):
    if "-" in steps: # sequential list of steps
        start_step = int(steps.split("-")[0])
        stop_step = int(steps.split("-")[1])
        return [(step + start_step) for step in range(stop_step)]
    elif "," in steps: # list of steps
        return [int(step) for step in steps.split(",")]
    else: # singular step
        return int(steps)

if __name__ == "__main__":
    print("Running DomCycle!")
    print("Start time:", datetime.datetime.now())
    run_steps()
    print("\n\nFinished running DomCycle!")
    print("Stop time:", datetime.datetime.now())