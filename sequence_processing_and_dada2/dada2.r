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

seqDir = "R1_R2_switched/itsxpress"
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
saveRDS(out2, "intermediate_RDS/read_filtering_read_counts_2.rds")

# sort filtered read files
itsFs.len <- sort(list.files(path.len, pattern = "_R1_001.fastq.gz", full.names = TRUE))
itsRs.len <- sort(list.files(path.len, pattern = "_R2_001.fastq.gz", full.names = TRUE))


#Vis read quality of its-extracted reads
#If running on a large sample set should index the filename object to [1:25] otherwise will be unreadable
pdf("dada2_processing_tables_figs/read_quality_post_its_extraction.pdf")
print(plotQualityProfile(itsFs.len[1:25]))
print(plotQualityProfile(itsRs.len[1:25]))
dev.off()

#This is the start of the core algorithm pipeline
#At this point the tutorial at https://benjjneb.github.io/dada2/tutorial.html is likely more informative than the ITS specific tutorial

#Learn the error rates
errF <- learnErrors(itsFs.len, multithread = TRUE)
errR <- learnErrors(itsRs.len, multithread = TRUE)

#Viz
pdf("dada2_processing_tables_figs/error_rate_graphs.pdf")
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
saveRDS(seqtab.nochim, file = "intermediate_RDS/dada2_seq_table_no_chim.rds")

#These files saved only for sequence counting purposes in "dada2_tables_to_file
getN <- function(x) sum(getUniques(x))

denoisedF.getN = data.frame(denoisedF = sapply(dadaFs, getN), sample = rownames(data.frame(sapply(dadaFs, getN))))
saveRDS(denoisedF.getN, "intermediate_RDS/denoisedF.getN.df.rds")
denoisedR.getN = data.frame(denoisedR = sapply(dadaRs, getN), sample = rownames(data.frame(sapply(dadaRs, getN))))
saveRDS(denoisedR.getN, "intermediate_RDS/denoisedR.getN.df.rds")
mergers.getN = data.frame(merged = sapply(mergers, getN), sample = rownames(data.frame(sapply(mergers, getN))))
saveRDS(mergers.getN, "intermediate_RDS/mergers.getN.df.rds")

