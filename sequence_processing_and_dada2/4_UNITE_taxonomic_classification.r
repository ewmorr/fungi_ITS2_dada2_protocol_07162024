library(dada2)

args = commandArgs(trailingOnly=TRUE)

seqTabPath = args[1]
outDir = args[2]
unite.ref = args[3]

seqtab.nochim = readRDS(file.path(seqTabPath))
taxa.w_bootstraps <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = 8, tryRC = TRUE, outputBootstraps = T)

saveRDS(taxa.w_bootstraps, file.path(outDir, "taxa_w_bootstraps.rds"))


#write results

# tax table:
asv_tax <- taxa.w_bootstraps$tax
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, file.path(outDir, "ASVs_taxonomy.tsv"), sep="\t", quote=F, col.names=NA)
#bottstraps
asv_tax_boot = taxa.w_bootstraps$boot
row.names(asv_tax_boot) <- sub(">", "", asv_headers)
write.table(asv_tax_boot, file.path(outDir, "ASVs_taxonomy_bootstrapVals.tsv"), sep="\t", quote=F, col.names=NA)
