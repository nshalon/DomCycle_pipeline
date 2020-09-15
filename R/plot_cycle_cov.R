library(optparse)
library(dplyr)
library(ggplot2)
library(reshape2)

#take in arguments
option_list = list(
  make_option(c("-i", "--inputfile"), action="store", default=NA, type='character',
              help="cycle coverages"),
  make_option(c("-p", "--cyclepaths"), action="store", default=NA, type='character',
              help="cycle contig paths"),
  make_option(c("-s", "--stats"), action="store", default=NA, type='character',
              help="cycle stats"),
  make_option(c("-o", "--outputdirectory"), action="store", default=NA, type='character',
              help="output directory")
)
#parse arguments to list
opt = parse_args(OptionParser(option_list=option_list))

covs = read.table(opt$inputfile, header=T)
contigs = read.table(opt$cyclepaths, header=T)
cycles = unique(contigs$cycle)
stats = read.table(opt$stats, header=T)

#iterate over cycles
for (cyc in cycles) {
  cyc_coverage = covs[which(covs$cycle == cyc),]
  cycle_stats = stats[stats$cycle == cyc,]
  # pval_pass = match(T, as.numeric(c(cycle_stats$pval_1, cycle_stats$pval_2, cycle_stats$pval_3)) < 0.01)
  df = melt(cyc_coverage, id.vars = "base", measure.vars = c("support","in_cov","out_cov","in_paired_weird","out_paired_weird","in_singleton","out_singleton"))
  df$variable = factor(df$variable, levels = c("support","out_cov","out_paired_weird","out_singleton","in_cov","in_paired_weird","in_singleton"))
  colnames(cyc_coverage) = colnames(covs)
  cycle_contigs = contigs[contigs$cycle == cyc, ]
  pass = match(T, c(cycle_stats$pval_1, cycle_stats$pval_2, cycle_stats$pval_3, 10) > 0.01) - 1
  plt = ggplot(df, aes(x = base, y = value, color = variable)) +
    geom_line() +
    scale_color_manual(values = c("grey1","red1","green3","lightpink2", "blue1","lightseagreen", "darkmagenta")) + 
    guides(color = guide_legend(override.aes = list(size = 0.24))) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 4)) +
    geom_vline(xintercept = (cycle_contigs$cum_sum), colour="black", linetype = "dotted", size=0.5) +
    annotate("text", x=(cycle_contigs$cum_sum + 0.5*cycle_contigs$contig_len), y=c(-5), label=(cycle_contigs$renamed_contig), size=1.5) +
    ggtitle(paste(cyc,"pval_thresh", pass,  "LBC", cycle_stats$lower_bound_cov)) + 
    ylab("coverage")
  ggsave(filename = paste(cyc,"coverage.pdf"), plot=plt, path=opt$outputdirectory, width=7,height=3.5)
}





