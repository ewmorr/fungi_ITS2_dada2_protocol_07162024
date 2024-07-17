library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")


workDir = args[1]

seqDir = file.path(workDir, "reads")
#list.files(seqDir)

#parse and sort file names, adjust regex as needed
fnFs <- sort(list.files(seqDir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(seqDir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Taylor primers
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
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

#remove reads with Ns
seqDir = file.path(workDir, "filtN")
if(!dir.exists(seqDir)) dir.create(seqDir)

fnFs.filtN <- file.path(seqDir, basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(seqDir, basename(fnRs))
out.filtN = filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#list.files(seqDir)
#Reset the filename lists in case files are lost at trimming
#parse and sort file names, adjust regex as needed
fnFs.filtN <- sort(list.files(seqDir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs.filtN <- sort(list.files(seqDir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Count occurences of the primers
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

outDir = file.path(workDir, "dada2_processing_tables_figs")
if(!dir.exists(outDir)) dir.create(outDir)

x = rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
write.csv(x, file.path(outDir, "pre_primerTrim_primer_check.csv"))

#Vis read quality
#If running on a large sample set should index the filename object to [1:25] otherwise will be unreadable

pdf(file.path(outDir, "read_quality_pre_cutadapt.pdf"))
print(plotQualityProfile(fnFs.filtN[1:25]))
print(plotQualityProfile(fnRs.filtN[1:25]))
dev.off()
