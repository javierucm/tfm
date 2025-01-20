
data <- read.table("TFM_JLR_Long10MonthYoungvsOldTreated.cntTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",3),rep("CGroup",4)))
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)
res <- results(dds)
write.table(res, file="TFM_JLR_Long10MonthYoungvsOldTreated_gene_TE_analysis.txt", sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 0.000000)), ]
write.table(resSig, file="TFM_JLR_Long10MonthYoungvsOldTreated_sigdiff_gene_TE.txt",sep="\t", quote=F)
