library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

args = commandArgs(trailingOnly=TRUE)

workDir = args[1]
seqDir = file.path(workDir, "cutadapt")
#list.files(seqDir)

fnFs.cut <- sort(list.files(seqDir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs.cut <- sort(list.files(seqDir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

FWD = "AGCCTCCGCTTATTGATATGCTTAART" #28S primer
REV = "AACTTTYRRCAAYGGATCWCT" #5.8S primer

#reverse, complement, and RC the primers
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    library(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
    RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

checksDir = file.path(workDir, "dada2_processing_tables_figs")

x = rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
    write.csv(x, file.path(checksDir, "post_primerTrim_primer_check.csv"))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_L00")[[1]][1]
sample.names <- unname(sapply(fnFs.cut, get.sample.name))
head(sample.names)

n_seqs = ifelse(length(fnRs.cut) >= 25, 25, length(fnRs.cut) )

#Vis read quality
pdf(file.path(checksDir, "read_quality_pre_dada2_qual_filtering.pdf"))
print(plotQualityProfile(fnFs.cut[1:n_seqs]) )
print(plotQualityProfile(fnRs.cut[1:n_seqs]) )
dev.off()

########################
#Quality filter

path.qual <- file.path(seqDir, "qual_filter")
if(!dir.exists(path.qual)) dir.create(path.qual)
fnFs.qual <- file.path(path.qual, basename(fnFs.cut))
fnRs.qual <- file.path(path.qual, basename(fnRs.cut))

out.qual <- filterAndTrim(fnFs.cut, fnFs.qual, fnRs.cut, fnRs.qual, maxN = 0, maxEE = c(2, 2),
truncQ = 2, minLen = 10, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
saveRDS(out.qual, file = file.path(checksDir, "qual_read_counts.rds"))
# sort filtered read files
fnFs.qual <- sort(list.files(path.qual, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs.qual <- sort(list.files(path.qual, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Vis read quality of quality filtered reads
#If running on a large sample set should index the filename object to [1:25] otherwise will be unreadable
n_seqs = ifelse(length(itsFs.qual) >= 25, 25, length(itsFs.qual) )

pdf(file.path(checksDir, "read_quality_post_qual_filter.pdf"))
print(plotQualityProfile(fnFs.qual[1:n_seqs]))
print(plotQualityProfile(fnRs.qual[1:n_seqs]))
dev.off()

