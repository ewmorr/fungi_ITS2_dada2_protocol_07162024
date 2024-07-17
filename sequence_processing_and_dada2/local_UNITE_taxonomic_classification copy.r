library(dada2)


seqtab.nochim = readRDS("data/vol_test/dada2_out/dada2_seq_table_no_chim.rds")

#taxa.w_bootstraps <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = 8, tryRC = TRUE, #outputBootstraps = T)
#saveRDS(taxa.w_bootstraps, file.path(outDir, "taxa_w_bootstraps.rds"))


#write results

#reading the files written on premise to rewrite with asv names
asv_tax = read.csv("data/vol_test/dada2_out/ASVs_taxonomy.tsv", sep = "\t")
asv_tax_boot = read.csv("data/vol_test/dada2_out/ASVs_taxonomy_bootstrapVals.tsv", sep = "\t")

# tax table:
#asv_tax <- taxa.w_bootstraps$tax
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}

row.names(asv_tax) <- sub(">", "", asv_headers)

write.table(asv_tax, "data/vol_test/dada2_out/ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
#bottstraps
#asv_tax_boot = taxa.w_bootstraps$boot
row.names(asv_tax_boot) <- sub(">", "", asv_headers)
write.table(asv_tax_boot, "data/vol_test/dada2_out/ASVs_taxonomy_bootstrapVals.tsv", sep="\t", quote=F, col.names=NA)
