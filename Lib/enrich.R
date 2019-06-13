#!/usr/bin/env Rscript
library("clusterProfiler")      
library("topGO")                 
library("Rgraphviz")         
library("pathview")                
library("org.Hs.eg.db")
argv=commandArgs(TRUE)
outdir=argv[1]
if(file.exists(outdir)) {
    setwd(outdir)
} else {
    dir.create(outdir)
    setwd(outdir)
}
if(file.exists(argv[2])) {
    result=read.table(argv[2], sep="\t", header = T)
} else {
    result=read.table(file('stdin'), sep="\t", header = F)
}
genes=as.character(result$geneNames)
if(length(argv)==3) {
    genetype=argv[3]
} else {
    genetype="SYMBOL"
}
entrez_id=mapIds(x=org.Hs.eg.db,keys=genes,keytype=genetype,column ="ENTREZID")
entrez_id=na.omit(entrez_id)
erich.go.ALL=enrichGO(gene=entrez_id,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="ALL",pvalueCutoff=0.05,qvalueCutoff=0.05)
write.csv(summary(erich.go.ALL),paste(basename(outdir), ".G-enrich.csv", sep=""),row.names =F)
pdf(paste(basename(outdir), ".G-enrich.pdf", sep=""),height = 10,width = 10)
barplot(erich.go.ALL,drop=TRUE,showCategory = 25)
dev.off()
KEGG=enrichKEGG(gene=entrez_id,pvalueCutoff=0.05,qvalueCutoff=0.1)
write.csv(summary(KEGG),paste(basename(outdir), ".K-enrich.csv", sep=""),row.names =F)
pdf(paste(basename(outdir), ".K-enrich.pdf", sep=""), height = 10,width = 10)
dotplot(KEGG,showCategory = 12) 
dev.off()
