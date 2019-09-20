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
FPKM_gene <- data.frame(gexpr(bg_filt))
FPKM_transcript <- data.frame(texpr(bg_filt, meas="FPKM"))
needName <- sapply(colnames(FPKM_gene), function(x){strsplit(x,split = "\\.")[[1]][2]})
geneName <- ballgown::geneNames(bg_filt)
colnames(FPKM_gene) <- needName
#rownames(FPKM_gene) <- geneName
colnames(FPKM_transcript) <- needName
#FPKM_gene <- arrange(data.frame(geneNames=ballgown::geneNames(bg_filt),FPKM_gene))
FPKM_transcript <- merge(data.frame(ballgown::transcriptNames(bg_filt)), FPKM_transcript, by=0, all.y=TRUE)

results_transcripts <- arrange(data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt),transcriptNames=ballgown::transcriptNames(bg_filt),stattest(bg_filt, feature="transcript",covariate="class", getFC=TRUE, meas="FPKM")),pval)
genes <- stattest(bg_filt, feature="gene",covariate="class", getFC=TRUE, meas="FPKM")
indices <- match(genes$id, texpr(bg, 'all')$gene_id)
gene_names_for_result <- texpr(bg, 'all')$gene_name[indices]
results_genes <- arrange(data.frame(geneNames=gene_names_for_result, genes),pval)

write.table(FPKM_gene, paste(outdir, "FPKM_gene.txt", sep="/"), row.names=TRUE, sep='\t', quote=FALSE, col.names=NA)
write.table(FPKM_transcript, paste(outdir, "FPKM_transcript.txt", sep="/"), row.names=F, sep='\t', quote=FALSE)
write.table(results_transcripts, paste(outdir, "transcript_results.csv", sep="/"), row.names=FALSE, sep='\t', quote=FALSE)
write.table(results_genes, paste(outdir, "results_genes.csv", sep="/"), row.names=FALSE, sep='\t', quote=FALSE)
write.table(subset(results_transcripts, results_transcripts$pval<0.05), paste(outdir, "transcript_results_0.05.csv", sep="/"), row.names=FALSE, sep='\t', quote=FALSE)
write.table(subset(results_genes, results_genes$pval<0.05), paste(outdir, "results_genes_0.05.csv", sep="/"), row.names=FALSE, sep='\t', quote=FALSE)
