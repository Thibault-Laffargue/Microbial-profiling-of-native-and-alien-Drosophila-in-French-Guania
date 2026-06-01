setwd("~/Mbiome_16s_and_OL4/basecalling/Mbiome_OL4_decompressed.fastq_1731093916/cleaned_seq2")

stats_reads=read.table("quality_filter_stats.tsv",sep="\t",header=T)

colnames(stats_reads)

# remove NA 

stats_reads.2=stats_reads[-which(is.na(stats_reads$sample_name)),]

# remove pools <=400 reads (except negative controls)

stats_reads.3=stats_reads.2[!(grepl("Pool", stats_reads.2$sample_name) &stats_reads.2$cleaned_reads < 200),]


# Statistics

    # Keep only Pools
stats_reads_pools=stats_reads.2[grep("Pool",stats_reads.2$sample_name),]

#stats_reads_controls=stats_reads.3[grep("Pool",stats_reads.3$sample_name,invert=TRUE),]
#stats_reads_control_negative=stats_reads_controls[-1,] #removal of positive control

pdf("read_quality_plots.pdf", width=8, height=6)

hist(stats_reads_pools$cleaned_reads,col="grey",nclass=40, xlim=c(0, 10000))
abline(v=400,col="red")

boxplot(stats_reads_pools$cleaned_reads,col="grey",ylim=c(0,9000))
abline(h=400,col="red")

dev.off()
