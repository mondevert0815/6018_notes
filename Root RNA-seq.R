###Load EdgeR
library(edgeR)
###Import file
x = read.delim("bam_read_counts", row.names = 1, stringsAsFactors = FALSE)
head(x)
###Read targets file
targets = read.delim("targets.txt", row.names=1, sep="\t", stringsAsFactors = FALSE)
targets
targets$Treatment
###Add objects to a DGE list
y = DGEList(counts = x[,1:16], group = targets$Treatment, gene = data.frame(Length = x[ ,16]))
colnames(y) = targets$Label
dim(y)
y
###Filter
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep,]
dim(y)
###Recompute lib size
y$samples$lib.size <- colSums(y$counts)
###Normalization
y <- calcNormFactors(y)
y$samples
###Plot biological coefficient of variation
###x is bio var & y is tech var
plotMDS(y) 
###Common dispersion over all genes
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
### Plot of variance
plotBCV(y)  
###Genewise tests for differential expression
###Calculation is second term - first term
et <- exactTest(y, pair=c ("NT", "Conv"))
top <- topTags(et)
top 
###CPM values for top genes
cpm(y)[rownames(top), ]
###Total number of DE genes
### Which genes are up or down regulated, or have no difference
summary(de <- decideTestsDGE(et))  
###Plot log-fold-changes & highlight DE genes
###DE is the counts of -1,0 and 1s
detags <- rownames(y)[as.logical(de)]    
head(de)
plotSmear(et, de.tags=detags) 
abline(h=c(-1, 1), col="blue")
y$pseudo.counts ### counts after normalization 
y$pseudo.lib.size
y$AveLogCPM
