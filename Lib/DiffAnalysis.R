library('ballgown')
library('dplyr')
library('genefilter')
argv=commandArgs(TRUE)
indir <- argv[1]
pheno <- argv[2]
outdir <- argv[3]

raw_data <- arrange(read.table(pheno), V1)
pheno_data <- data.frame(sample = paste(raw_data$V1,raw_data$V2,sep="_"), class = raw_data$V2)
bg <- ballgown(dataDir = indir, samplePattern = "*", pData=pheno_data)
bg_filt <- subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)
FPKM <- gexpr(bg_filt)
needName <- sapply(colnames(FPKM), function(x){strsplit(x,split = "\\.")[[1]][2]})
colnames(FPKM) <- needName

results_transcripts <- arrange(data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt),stattest(bg_filt, feature="transcript",covariate="class", getFC=TRUE, meas="FPKM")),pval)
genes <- stattest(bg_filt, feature="gene",covariate="class", getFC=TRUE, meas="FPKM")
indices <- match(genes$id, texpr(bg, 'all')$gene_id)
gene_names_for_result <- texpr(bg, 'all')$gene_name[indices]
results_genes <- arrange(data.frame(geneNames=gene_names_for_result, genes),pval)

write.table(FPKM, paste(outdir, "FPKM.txt", sep="/"), row.names=TRUE, sep='\t', quote=FALSE, col.names=NA)
write.table(results_transcripts, paste(outdir, "transcript_results.csv", sep="/"), row.names=FALSE, sep='\t', quote=FALSE)
write.table(results_genes, paste(outdir, "results_genes.csv", sep="/"), row.names=FALSE, sep='\t', quote=FALSE)
write.table(subset(results_transcripts, results_transcripts$pval<0.05), paste(outdir, "transcript_results_0.05.csv", sep="/"), row.names=FALSE, sep='\t', quote=FALSE)
write.table(subset(results_genes, results_genes$pval<0.05), paste(outdir, "results_genes_0.05.csv", sep="/"), row.names=FALSE, sep='\t', quote=FALSE)
