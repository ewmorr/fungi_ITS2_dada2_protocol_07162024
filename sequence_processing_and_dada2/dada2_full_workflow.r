#This routine is from https://benjjneb.github.io/dada2/ITS_workflow.html with modifications
#itsxpress is run after primer removal and quality filtering but before dada2 core algorithm
#itsxpress is run outside of R
#this script is the core dada2 algorithm run after itsxpress

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

seqDir = "/Users/ericmorrison/GARNAS_neonectria_barcoding_012220/R1_R2_switched/itsxpress"
list.files(seqDir)

#parse and sort file names, adjust regex as needed
itsFs <- sort(list.files(seqDir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
itsRs <- sort(list.files(seqDir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#itsxpress outputs 0 len reads. filter for len
#make files
path.len <- file.path(seqDir, "gt_10")
if(!dir.exists(path.len)) dir.create(path.len)
itsFs.len <- file.path(path.len, basename(itsFs))
itsRs.len <- file.path(path.len, basename(itsRs))

#filter
out2 <- filterAndTrim(itsFs, itsFs.len, itsRs, itsRs.len, maxN = 0, maxEE = c(2, 2),
truncQ = 2, minLen = 10, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out2)

# sort filtered read files
itsFs.len <- sort(list.files(path.len, pattern = "_R1_001.fastq.gz", full.names = TRUE))
itsRs.len <- sort(list.files(path.len, pattern = "_R2_001.fastq.gz", full.names = TRUE))


#Vis read quality of its-extracted reads
#If running on a large sample set should index the filename object to [1:25] otherwise will be unreadable
pdf("read_quality_post_its_extraction.pdf")
print(plotQualityProfile(itsFs.len[1:25]))
print(plotQualityProfile(itsRs.len[1:25]))
dev.off()

#This is the start of the core algorithm pipeline
#At this point the tutorial at https://benjjneb.github.io/dada2/tutorial.html is likely more informative than the ITS specific tutorial

#Learn the error rates
errF <- learnErrors(itsFs.len, multithread = TRUE)
errR <- learnErrors(itsRs.len, multithread = TRUE)

#Viz
pdf("error_rate_graphs.pdf")
print(plotErrors(errF, nominalQ = TRUE))
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

#derep
#it seems like the derep step is not strictly necessary and is wrapped in the dada2 alg. (it's no longer a separate step in the main tutorial)
derepFs <- derepFastq(itsFs.len, verbose = TRUE)
derepRs <- derepFastq(itsRs.len, verbose = TRUE)

get.sample.name <- function(fname) strsplit(basename(fname), "_L00")[[1]][1]
sample.names <- unname(sapply(itsFs.len, get.sample.name))

names(derepFs) <- sample.names
names(derepRs) <- sample.names


#DADA2 alogorithm
#pooling samples is not default, but increases sensitivity to low abundance sequences shared across samples
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = T)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = T)

#merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, maxMismatch = 0)

#make seq table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, file = "dada2_seq_table_no_chim.rds")

#lengths of consensus sequences
table(nchar(getSequences(seqtab.nochim)))
write.csv(table(nchar(getSequences(seqtab.nochim))), "asv_lens.csv")

###
#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))

#####
#assign taxonomy against unite with RDB
print("taxonomy assignment with naive bayesian classifier")
unite.ref <- "/Users/ericmorrison/blast_dbs/sh_general_release_dynamic_02.02.2019/sh_general_release_dynamic_02.02.2019.fasta"
taxa.w_bootstraps <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE, outputBootstraps = T) #With latest R and tidyverse update this step can result in segfault and memalloc errors with tidyverse loaded. Need to load tidyverse after this step if needed (currently using for plotting and table manipulations
taxa.print <- taxa.w_bootstraps$tax  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL


################################
#Create workable objects (on hd)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa.w_bootstraps$tax
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
#bottstraps
asv_tax_boot = taxa.w_bootstraps$boot
row.names(asv_tax_boot) <- sub(">", "", asv_headers)
write.table(asv_tax_boot, "ASVs_taxonomy_bootstrapVals.tsv", sep="\t", quote=F, col.names=NA)

#####################
#tidyverse based operations performed after dada2

require(tidyr)
require(dplyr)
require(ggplot2)
source("../ggplot_theme.txt")

#read processing tracking

#modified from the original dada2 tracking table so can hadle different numbers of rows (samples that are dropped)
out = readRDS("read_filtering_read_counts.rds")
colnames(out) = c ("input", "filtered")
colnames(out2) = c ("itsxextracted", "gt_10_len")

rownames(out) = unname(sapply(rownames(out), get.sample.name))
rownames(out2) = unname(sapply(rownames(out2), get.sample.name))

track = full_join(data.frame(out, sample = rownames(out)), data.frame(out2, sample = rownames(out2)), by = "sample") %>%
full_join(., data.frame(denoisedF = sapply(dadaFs, getN), sample = rownames(data.frame(sapply(dadaFs, getN)))), by = "sample") %>%
full_join(., data.frame(denoisedR = sapply(dadaRs, getN), sample = rownames(data.frame(sapply(dadaRs, getN)))), by = "sample") %>%
full_join(., data.frame(merged = sapply(mergers, getN), sample = rownames(data.frame(sapply(mergers, getN)))), by = "sample") %>%
full_join(., data.frame(nonchim = rowSums(seqtab.nochim), sample = rownames(data.frame(rowSums(seqtab.nochim)))), by = "sample")

track.long = track %>% arrange(desc(input)) %>% gather(., "step", "count", -sample)
track.long$step = factor(track.long$step, levels = c("input", "filtered", "itsxextracted", "gt_10_len", "denoisedF", "denoisedR", "merged",
"nonchim"))

write.csv(track.long, "read_processing_tracking.csv")

#Plot read counts per step

p1 = ggplot(track.long %>% filter(step != "gt_10_len" & step != "denoisedR" & count > 0), aes(x = reorder(sample, -count), y = count, color = step)) +
geom_point() +
scale_y_log10() +
scale_color_manual(values = c(cbPalette, "dark grey", "black")) +
my_gg_theme +
theme(axis.text.x = element_blank()) +
labs(x = "sample")

pdf("read_counts_processing_steps_stacked.pdf", width = 20, height = 6)
print(p1)
dev.off()

p1 = ggplot(track.long %>% filter(step != "gt_10_len" & step != "denoisedR"), aes(x = reorder(sample, -count), y = count)) +
geom_point() +
#scale_y_log10() +
facet_wrap(~step) +
my_gg_theme +
theme(axis.text.x = element_blank()) +
labs(x = "sample")

pdf("read_counts_processing_steps_faceted.pdf", width = 10)
print(p1)
dev.off()

#plot ASV seq lens

p1 = ggplot(data.frame(table(nchar(getSequences(seqtab.nochim)))), aes(as.numeric(as.character(Var1)), as.numeric(Freq))) +
geom_point() +
my_gg_theme +
labs(x = "ASV length (bp)", y = "Frequency")

pdf("ASV_seq_lens.pdf")
print(p1)
dev.off()

#Look at Hypocreales
data.frame(taxa.print) %>% filter(Order == "o__Hypocreales")
