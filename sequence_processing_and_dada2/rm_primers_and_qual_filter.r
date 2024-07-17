#This routine is from https://benjjneb.github.io/dada2/ITS_workflow.html with modifications
#itsxpress is run after primer removal and quality filtering but before dada2 core algorithm
#itsxpress is run outside of R

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

seqDir = "/mnt/home/garnas/ewj4/EDRR_patho/reads"
list.files(seqDir)

#parse and sort file names, adjust regex as needed
fnFs <- sort(list.files(seqDir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(seqDir, pattern = "_R2_001.fastq.gz", full.names = TRUE))


#primer sequences for trimming
#The reads are renamed to reverse R1/R2 (so that itsxpress will work), so the FWD and Reverse are actualy reverse of Illumina direction
#REV = "NNNNNNAGCCTCCGCTTATTGATATGCTTAART" #LSU
#FWD = "NNNNNNAACTTTYRRCAAYGGATCWCT" #5.8S

#
REV = "AGCCTCCGCTTATTGATATGCTTAART" #28S primer
FWD = "AACTTTYRRCAAYGGATCWCT" #5.8S primer

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
FWD.orients


#remove reads with Ns
fnFs.filtN <- file.path(seqDir, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(seqDir, "filtN", basename(fnRs))
out.filtN = filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
saveRDS(out.filtN, "intermediate_RDS/read_filtering_read_counts.filtN.rds")
#Reset the filename lists in case files are lost at trimming
seqDir = "/mnt/home/garnas/ewj4/EDRR_patho/reads/filtN"
list.files(seqDir)

#parse and sort file names, adjust regex as needed
fnFs.filtN <- sort(list.files(seqDir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs.filtN <- sort(list.files(seqDir, pattern = "_R2_001.fastq.gz", full.names = TRUE))


#Count occurences of the primers
primerHits <- function(primer, fn) {
	# Counts number of reads in which the primer is found
	nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
	return(sum(nhits > 0))
}


x = rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
write.csv(x, "dada2_processing_tables_figs/pre_trimming_primer_check.csv")


#Vis read quality
#If running on a large sample set should index the filename object to [1:25] otherwise will be unreadable

pdf("dada2_processing_tables_figs/read_quality_pre_cutadapt.pdf")
print(plotQualityProfile(fnFs.filtN[1:25]))
print(plotQualityProfile(fnRs.filtN[1:25]))
dev.off()

#point to cutadapt
#This command should probably be run outside of R for running with slurm (so that can be threaded properly)
#Actually cutadapt does not appear to be multi-threaded so it's fine
cutadapt <- "/Users/ericmorrison/miniconda3/envs/qiime2-2019.4/bin/cutadapt"
system2(cutadapt, args = "--version")

#make files
path.cut <- file.path(seqDir, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN)) #THIS needs to be changed to .filtN!!!
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

#reverse complement primers
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)


# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# Run Cutadapt
#added min len filter bc cutadapt was leaving some empty entries
for(i in seq_along(fnFs.filtN)) {
	system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 libraryd to remove FWD and REV from reads
    "--minimum-length", 50,
	"-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
	fnFs.filtN[i], fnRs.filtN[i]), # input files
    stdout = "cutadapt.out") #this will overwrite for each new set of reads. Can remove this to print all to stdout
}


#recheck occurence
xx = rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
write.csv(xx, "dada2_processing_tables_figs/post_trimming_primer_check.csv")


# sort cutadapted read files
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_L00")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


#Vis read quality
pdf("dada2_processing_tables_figs/read_quality_pre_dada2_qual_filtering.pdf")
print(plotQualityProfile(cutFs[1:25]) )
print(plotQualityProfile(cutRs[1:25]) )
dev.off()

#everything looks good

#perform maxEE based read filtering

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#Should not really need to re-enforce the min len filter...
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)
saveRDS(out, file = "intermediate_RDS/read_filtering_read_counts.rds")

# sort filtered read files
filtFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Vis read quality of filtered reads
pdf("dada2_processing_tables_figs/read_quality_post_dada2_qual_filtering.pdf")
print(plotQualityProfile(filtFs[1:25]))
print(plotQualityProfile(filtRs[1:25]))
dev.off()


