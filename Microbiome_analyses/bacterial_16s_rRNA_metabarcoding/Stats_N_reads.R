dir="~/Mbiome_16s_and_OL4/basecalling/Mbiome_16s__1730125745/cleaned_seq2"
setwd(dir)

stats_reads=read.table("quality_filter_stats.tsv",sep="\t",header=T)

colnames(stats_reads)

# remove NA 

stats_reads.2=stats_reads[-which(is.na(stats_reads$sample_name)),]

# remove pools <=1000 reads (except negative controls)

stats_reads.3=stats_reads.2[!(grepl("Pool", stats_reads.2$sample_name) &stats_reads.2$cleaned_reads < 1000),]


# Statistics

    # Keep only Pools
stats_reads_pools=stats_reads.2[grep("Pool",stats_reads.2$sample_name),]

#stats_reads_controls=stats_reads.3[grep("Pool",stats_reads.3$sample_name,invert=TRUE),]
#stats_reads_control_negative=stats_reads_controls[-1,] #removal of positive control

pdf("read_quality_plots.pdf", width=8, height=6)

hist(stats_reads_pools$cleaned_reads, col="grey", nclass=20)
abline(v=1000, col="red")

boxplot(stats_reads_pools$cleaned_reads, col="grey", ylim=c(0,30000))
abline(h=1000, col="red")

dev.off()

#x11()
#hist(stats_reads_pools$cleaned_reads,col="grey",nclass=20)
#abline(v=1000,col="red")

#x11()
#boxplot(stats_reads_pools$cleaned_reads,col="grey",ylim=c(0,30000))
#abline(h=1000,col="red")



