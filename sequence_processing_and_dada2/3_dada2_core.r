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

args = commandArgs(trailingOnly=TRUE)

workdir = args[1]
seqDir = file.path(workdir, "itsxpress_out")
outDir = file.path(workdir, "dada2_core")
if(!dir.exists(outDir)) dir.create(outDir)
checksDir = file.path(workdir, "dada2_processing_tables_figs")

#list.files(seqDir)

#parse and sort file names, adjust regex as needed
itsFs <- sort(list.files(seqDir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
itsRs <- sort(list.files(seqDir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#itsxpress outputs 0 len reads. filter for len as well as normal filters
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
n_seqs = ifelse(length(fnFs.filtN) >= 25, 25, length(fnFs.filtN) )

pdf(file.path(checksDir, "read_quality_post_qual_filter.pdf"))
print(plotQualityProfile(itsFs.len[1:n]))
print(plotQualityProfile(itsRs.len[1:n]))
dev.off()

#This is the start of the core algorithm pipeline
#At this point the tutorial at https://benjjneb.github.io/dada2/tutorial.html is likely more informative than the ITS specific tutorial

#Learn the error rates
errF <- learnErrors(itsFs.len, multithread = TRUE)
errR <- learnErrors(itsRs.len, multithread = TRUE)

#Viz
pdf(file.path(checksDir, "error_rate_graphs.pdf"))
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
saveRDS(seqtab.nochim, file = file.path(outDir, "dada2_seq_table_no_chim.rds"))

#lengths of consensus sequences
table(nchar(getSequences(seqtab.nochim)))
write.csv(table(nchar(getSequences(seqtab.nochim))), file.path(outDir, "asv_lens.csv"))

###
#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(itsFs, getN), sapply(itsRs, getN), sapply(itsFs.len, getN), sapply(itsRs.len, getN), sapply(mergers, getN),
    rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("itsxpressF","itsxpressR", "filteredF", "filterR" "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, file.path(checksDir, "pre_primerTrim_primer_check.csv"))

################################
#Create workable objects (on hd)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(outDir, "ASVs.fa"))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(outDir, "ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)

