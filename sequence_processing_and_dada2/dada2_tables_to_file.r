library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

seqtab.nochim = readRDS("intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps = readRDS("intermediate_RDS/taxa_w_bootstraps.rds")
#taxa.print <- taxa.w_bootstraps$tax  # Removing sequence rownames for display only
#rownames(taxa.print) <- NULL

system("mkdir dada2_out")
################################
#Create workable objects (on hd)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "dada2_out/ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "dada2_out/ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa.w_bootstraps$tax
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "dada2_out/ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
#bootstraps
asv_tax_boot = taxa.w_bootstraps$boot
row.names(asv_tax_boot) <- sub(">", "", asv_headers)
write.table(asv_tax_boot, "dada2_out/ASVs_taxonomy_bootstrapVals.tsv", sep="\t", quote=F, col.names=NA)

##################################################
#tidyverse based operations performed after dada2#

require(tidyr)
require(dplyr)
require(ggplot2)
require(RColorBrewer)

source("~/ggplot_theme.txt")

##########################
#read processing tracking#

get.sample.name <- function(fname) strsplit(basename(fname), "_L00")[[1]][1]

#modified from the original dada2 tracking table so can handle different numbers of rows (samples that are dropped)
#out.filtN = readRDS("intermediate_RDS/read_filtering_read_counts.filtN.rds")
#out = readRDS("intermediate_RDS/read_filtering_read_counts.rds")
out2 = readRDS("intermediate_RDS/read_filtering_read_counts_2.rds")
#colnames(out.filtN) = c("input", "N_filtered")
#colnames(out) = c ("cutadapt", "quality_filtered")
colnames(out2) = c ("itsxextracted", "gt_10_len")


### The commented rows are not run because the run and dada2 pipeline comparisons were performed after the quality filtering and read extraction steps which are performed on a per read basis
#rownames(out.filtN) = unname(sapply(rownames(out.filtN), get.sample.name))
#rownames(out) = unname(sapply(rownames(out), get.sample.name))
rownames(out2) = unname(sapply(rownames(out2), get.sample.name))

#Track reads through the pipeline

track = full_join(data.frame(out2, sample = rownames(out2)), readRDS("intermediate_RDS/denoisedF.getN.df.rds"), by = "sample") %>%
full_join(., data.frame(readRDS("intermediate_RDS/denoisedR.getN.df.rds")), by = "sample") %>%
full_join(., data.frame(readRDS("intermediate_RDS/mergers.getN.df.rds")), by = "sample") %>%
full_join(., data.frame(nonchim = rowSums(seqtab.nochim), sample = rownames(data.frame(rowSums(seqtab.nochim)))), by = "sample") %>%
replace(., is.na(.), 0)

#alternate join for processing full tables without using sample names (ordered)
track = cbind(data.frame(out2, sample = rownames(out2)), readRDS("intermediate_RDS/denoisedF.getN.df.rds")[,1]) %>%
    cbind(.,readRDS("intermediate_RDS/denoisedR.getN.df.rds")[,1]) %>%
    cbind(.,readRDS("intermediate_RDS/mergers.getN.df.rds")[,1]) %>%
    cbind(.,data.frame(nonchim = rowSums(seqtab.nochim)))
colnames(track) = c ("itsxextracted", "gt_10_len", "sample", "denoisedF", "denoisedR", "merged", "nonchim")



track.long = track %>% arrange(desc(itsxextracted)) %>% gather(., "step", "count", -sample)
track.long$step = factor(track.long$step, levels = c("itsxextracted", "gt_10_len", "denoisedF", "denoisedR", "merged",
"nonchim"))

write.csv(track.long, "dada2_processing_tables_figs/read_processing_tracking.csv")
write.csv(track, "dada2_processing_tables_figs/read_processing_tracking.wide_format.csv")

#Plot read counts per step

p1 = ggplot(track.long %>% filter(step != "gt_10_len" & step != "denoisedR"), aes(x = reorder(sample, -count), y = count + 1, color = step)) +
geom_point(position = position_jitter(width = 0.75, height = 0)) +
scale_y_log10() +
scale_color_manual(values = c(cbPalette, "dark grey", "black")) +
my_gg_theme +
theme(axis.text.x = element_blank()) +
labs(x = "sample")

pdf("dada2_processing_tables_figs/read_counts_processing_steps_stacked.pdf", width = 20, height = 6)
print(p1)
dev.off()

palette_len = track.long %>% filter(step != "gt_10_len" & step != "denoisedR") %>% select(sample) %>% unique %>% rownames %>% length
getPalette = colorRampPalette(brewer.pal(n = 12, name = "Paired"))

p1 = ggplot(track.long %>% filter(step != "gt_10_len" & step != "denoisedR"), aes(x = step, y = count + 1, group = sample, color = sample)) +
geom_line() +
scale_y_log10() +
my_gg_theme +
scale_color_manual(values = getPalette(palette_len), guide = F) +
labs(x = "") +
theme(axis.text.x =element_text(angle = 35, hjust = 1))

pdf("dada2_processing_tables_figs/read_counts_processing_steps_by_sample.pdf", width = 12, height = 6)
print(p1)
dev.off()

p1 = ggplot(track.long %>% filter(step != "gt_10_len" & step != "denoisedR"), aes(x = reorder(sample, -count), y = count)) +
geom_point() +
#scale_y_log10() +
facet_wrap(~step) +
my_gg_theme +
theme(axis.text.x = element_blank()) +
labs(x = "sample")

pdf("dada2_processing_tables_figs/read_counts_processing_steps_faceted.pdf", width = 10)
print(p1)
dev.off()

palette_len = track.long %>% filter(step != "gt_10_len" & step != "denoisedR") %>% select(sample) %>% unique %>% rownames %>% length
getPalette = colorRampPalette(brewer.pal(n = 12, name = "Paired"))

p1 = ggplot(track.long %>% filter(step != "gt_10_len" & step != "denoisedR"), aes(x = step, y = count + 1, group = sample, color = sample)) +
geom_line() +
scale_y_log10() +
my_gg_theme +
scale_color_manual(values = getPalette(palette_len), guide = F) +
labs(x = "") +
theme(axis.text.x =element_text(angle = 35, hjust = 1))

pdf("dada2_processing_tables_figs/read_counts_processing_steps_by_sample.pdf", width = 12, height = 6)
print(p1)
dev.off()

###################
#plot ASV seq lens#

#lengths of consensus sequences
table(nchar(getSequences(seqtab.nochim)))
write.csv(table(nchar(getSequences(seqtab.nochim))), "dada2_processing_tables_figs/asv_lens.csv")

#plot
p1 = ggplot(data.frame(table(nchar(getSequences(seqtab.nochim)))), aes(as.numeric(as.character(Var1)), as.numeric(Freq))) +
geom_point() +
my_gg_theme +
labs(x = "ASV length (bp)", y = "Frequency")

pdf("dada2_processing_tables_figs/ASV_seq_lens.pdf")
print(p1)
dev.off()
