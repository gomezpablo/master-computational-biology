## LOADING DATA ##

counts = read.table("counts.txt", header=F, row.names=1) # Load the raw counts table
colnames = c("Normal","Normal","High","High") # names for column names
my.design <- data.frame(row.names = colnames( counts ),
                        group = c("Normal","Normal","High","High")
) # our experiment design for DESeq2 analysis

## installing DESeq2 ## # (if necessary)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# BiocManager::install("DESeq2")

library("DESeq2") # Load the DESeq2 package

#import data matrix from NtcA results to generate DeSeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts, colData = my.design, design = ~ group + group:group)


#DEF analysis
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast=c("group","Normal","High"))
res

res = na.omit (res) 

##Histogram of p-values from the call to nbinomTest.
hist(res$pvalue, breaks=100, col="skyblue", border="slateblue", main="")

##Histogram of p-adj from the call to nbinomTest.
hist(res$padj, breaks=100, col="darkorange1", border="darkorange4", main="")


#plotMA
plotMA(res, main="DESeq2", ylim=c(-8,8))

#identifie genes in the plot
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]


#rlogtranformation for PCA analysis
rld <- rlog(dds)
head(assay(rld), 100)
plotPCA(rld, intgroup=c("Location", "group"))



#heatmap
resOrdered <- rld[order(res$log2FoldChange),]
select_genes<-rownames(subset(resOrdered, res$padj< 0.01))
mat<-assay(rld)[select_genes,]
mat<-mat-rowMeans(mat)
heatmap(mat)




library("genefilter")
rld <- rlog(dds, blind=FALSE)
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),50)
heatmap( assay(rld)[ topVarGenes, ], scale="row")




##filter for significant genes, 
##according to some chosen threshold for the false dicovery rate (FDR)



resSig = res[res$padj < 0.01, ]


## QUESTION 2 ##
filtered_padj <- subset(res, res$padj < 0.01)[c("pvalue", "padj", "log2FoldChange")]
