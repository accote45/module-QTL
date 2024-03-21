



####################################################################################
### read in CMC expr, complete preliminary QC

library(tidyverse)
library(biomaRt)
library(plyr)
library(limma)
library(Glimma)
library(edgeR)
library('variancePartition')
library('doParallel')
library(data.table)
library(corrplot)
library(RColorBrewer)

library(MatrixEQTL)

COUNT <- fread("expr_files/txicounts_+HBCC.csv.gz",header=T,sep=",",stringsAsFactors=F)
demos <- read.csv("demos/demos.master.EA.CMC1.amanda.ancestry.only.csv",header=T,stringsAsFactors=F) 

# Create RIN^2 and clustered lib batch covariates
demos$RIN2 <- demos$rnaSeq_isolation.RIN^2

cBatchToBatch=as.matrix(list("0"=c("11/11/13","11/26/13"), "A"=c(17,18,10,"8/28/13"), "B"=c(14,7,11,6,16,25), "C"=c(5,1,3), "D"=c(2,28,8,12,4), "E"=c(20,24,26,27,9), "F"=c(13,15,22,21,23), "G"=c("10/15/13",19), "H"=c("10/9/13")))
key <- as.data.frame(cBatchToBatch) %>% rownames_to_column("clusterLIB")
key$Library_Batch.2 <- sapply(key$V1, paste, collapse=",")
key <- key %>% separate_rows(Library_Batch.2, sep = "," )
key$V1 <- NULL
demos$Library_Batch.2 <- substring(demos$Library_Batch, 3)

demos <- merge(demos, key, by.x="Library_Batch.2")


# Define biomart object
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


# Query biomart
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "percentage_gene_gc_content", "gene_biotype", "chromosome_name", "start_position", "end_position"),
                       filters = "ensembl_gene_id", values = COUNT$V1,
                       mart = mart, uniqueRows = TRUE)

Ensemble2HGNC.dedup <- Ensemble2HGNC[!duplicated(Ensemble2HGNC$ensembl_gene_id),]


COUNT$ensembl_gene_id <- COUNT$V1
COUNT = join(COUNT, Ensemble2HGNC.dedup, by="ensembl_gene_id", type="left") %>% filter(chromosome_name != "X" & chromosome_name != "Y" & chromosome_name !="MT")
COUNT <- as.data.frame(COUNT)

# Keep only individuals with info in demos file
index <- c(demos$RNAseq.Sample_RNA_ID, "ensembl_gene_id")
count.trim <- COUNT[,(names(COUNT) %in% index)] %>% column_to_rownames("ensembl_gene_id")

## To remove participants missing one or more covars
covars <- c("Genotypes.Genotyping_Sample_ID","RNAseq.Sample_RNA_ID",'Dx', 'Sex', 'Age_of_Death', 'Institution', 'PMI_.in_hours.', "rnaSeq_isolation.RIN", "RIN2", "EV.1", "EV.2", "EV.3", "EV.4", "EV.5", "clusterLIB")
demos.fin <- demos[,covars] %>% na.omit()
keep <- demos.fin$RNAseq.Sample_RNA_ID
rownames(demos.fin) <- NULL
count.master <- count.trim[ , (names(count.trim) %in% keep)]


### Normalisation

# Sort files, make DGEList object
demos.sort <- demos.fin[order(demos.fin$RNAseq.Sample_RNA_ID),] 
count.sort <- count.master[ , order(names(count.master))]

y = DGEList( count.sort, genes=rownames(count.sort) )
y <- calcNormFactors(y, method="TMM")
dim(y)

tmm.cpm <- cpm(y, normalized.lib.sizes=TRUE, log=TRUE)


### filter genes based on threshold of >0.5 CPM in >/= 30% of samples
filt.genes.cpm <- rownames(tmm.cpm[rowSums(tmm.cpm>=0.5) >= 0.3*ncol(tmm.cpm),])
tmm.cpm.expr <- tmm.cpm[rownames(tmm.cpm) %in% filt.genes.cpm,]

y <- y[filt.genes.cpm,,keep.lib.sizes=F]

### Removal of outlier samples 
# Find principal components of expression to plot
PC <- prcomp(t(tmm.cpm.expr), scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(RNAseq.Sample_RNA_ID=rownames(PC$x), 
                       PC1=PC$x[,1], 
                       PC2=PC$x[,2])
plotdata <- merge(plotdata, demos.sort, by="RNAseq.Sample_RNA_ID", all.y=TRUE)

pdf("~/www/PCAplot_cmc1_euro_dlpfc_unadjusted.pdf")
p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(colour=factor(Institution), fill=factor(Institution), shape=Dx, size=Age_of_Death)) 
p <- p + scale_shape_manual(values=c(3,22,21,24))
p <- p + theme_bw() %+replace% theme(legend.position="right")
p
dev.off()



# For CMC1 cohort only, removed single outlier based on inspection of pc plot
plotdata <- plotdata[order(plotdata$PC2),]
tail(plotdata)

indremove <- c("PENN_RNA_PFC_56")
y = y[,!colnames(y) %in% indremove]
demos.sort = demos.sort[!demos.sort$RNAseq.Sample_RNA_ID %in% indremove,]

# check for outlier samples by inter array correlation (samples with IAC < 3SDs below mean IAC for the dataset will be removed) ### source = https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/HumanBrainTranscriptome/Identification%20and%20Removal%20of%20Outlier%20Samples%20-%20Illumina.pdf
IAC=cor(tmm.cpm.expr,use="p") 
library(WGCNA)

png("~/www/iac_histogram_cmc1.png")
hist(IAC,sub=paste("mean=",format(mean(IAC[upper.tri(IAC)]),digits=3))) 
dev.off()

sampleTree = hclust(as.dist(1-IAC), method = "average");
png("~/www/iacbased_outlier_check_cmc1.png")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
abline(h=0.1,col="red") 
dev.off()

groups.2 <- cutree(sampleTree,2)
table(groups.2)
groups.2[order(groups.2)]

meanIAC=apply(IAC,2,mean)
sdCorr=sd(meanIAC)
numbersd=(meanIAC-mean(meanIAC))/sdCorr
png('~/www/iac_sd_graph_cmc1.png')
plot(numbersd)
abline(h=-3)
dev.off()

sdout=-3
outliers=dimnames(tmm.cpm.expr)[[2]][numbersd<sdout]
outliers

removevec=c("MSSM_RNA_PFC_139", "MSSM_RNA_PFC_133", "MSSM_RNA_PFC_163", "MSSM_RNA_PFC_299", "MSSM_RNA_PFC_60", "PENN_RNA_PFC_56", "PENN_RNA_PFC_71")
overlap1=is.element(dimnames(tmm.cpm.expr)[[2]],removevec)
datrest2=tmm.cpm.expr[,!overlap1]
dim(datrest2)

IAC=cor(datrest2,use="p")
png("~/www/iac_histogram_afteroutlierremoval_cmc1.png")
hist(IAC,sub=paste("mean=",format(mean(IAC[upper.tri(IAC)]),digits=3))) 
dev.off()

cluster1=hclust(as.dist(1-IAC),method="average")
png("~/www/iacbased_outlier_check_round2_cmc1.png")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(cluster1, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
abline(h=0.1,col="red") 
dev.off()

groups.2 <- cutree(cluster1,2)
table(groups.2)
groups.2[order(groups.2)]

meanIAC=apply(IAC,2,mean)
sdCorr=sd(meanIAC)
numbersd=(meanIAC-mean(meanIAC))/sdCorr
png('~/www/iac_sd_graph_round2_cmc1.png')
plot(numbersd)
abline(h=-3)
dev.off()

sdout=-3
outliers=dimnames(datrest2)[[2]][numbersd<sdout]
outliers


removevec=c("MSSM_RNA_PFC_139", "MSSM_RNA_PFC_133", "MSSM_RNA_PFC_163", "MSSM_RNA_PFC_299", "MSSM_RNA_PFC_60", "PENN_RNA_PFC_56", "PENN_RNA_PFC_71","MSSM_RNA_PFC_304","MSSM_RNA_PFC_348","MSSM_RNA_PFC_89")
overlap1=is.element(dimnames(tmm.cpm.expr)[[2]],removevec)
datrest2=tmm.cpm.expr[,!overlap1]
dim(datrest2)

IAC=cor(datrest2,use="p")
png("~/www/iac_histogram_afteroutlierremoval_round3_cmc1.png")
hist(IAC,sub=paste("mean=",format(mean(IAC[upper.tri(IAC)]),digits=3))) 
dev.off()

cluster1=hclust(as.dist(1-IAC),method="average")
png("~/www/iacbased_outlier_check_round3_cmc1.png")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(cluster1, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
abline(h=0.1,col="red") 
dev.off()

groups.2 <- cutree(cluster1,2)
table(groups.2)
groups.2[order(groups.2)]

meanIAC=apply(IAC,2,mean)
sdCorr=sd(meanIAC)
numbersd=(meanIAC-mean(meanIAC))/sdCorr
png('~/www/iac_sd_graph_round3_cmc1.png')
plot(numbersd)
abline(h=-3)
dev.off()

sdout=-3
outliers=dimnames(datrest2)[[2]][numbersd<sdout]
outliers



removevec=c("MSSM_RNA_PFC_329","MSSM_RNA_PFC_139", "MSSM_RNA_PFC_133", "MSSM_RNA_PFC_163", "MSSM_RNA_PFC_299", "MSSM_RNA_PFC_60", "PENN_RNA_PFC_56", "PENN_RNA_PFC_71","MSSM_RNA_PFC_304","MSSM_RNA_PFC_348","MSSM_RNA_PFC_89")
overlap1=is.element(dimnames(tmm.cpm.expr)[[2]],removevec)
datrest2=tmm.cpm.expr[,!overlap1]
dim(datrest2)

IAC=cor(datrest2,use="p")
png("~/www/iac_histogram_afteroutlierremoval_round4_cmc1.png")
hist(IAC,sub=paste("mean=",format(mean(IAC[upper.tri(IAC)]),digits=3))) 
dev.off()

cluster1=hclust(as.dist(1-IAC),method="average")
png("~/www/iacbased_outlier_check_round4_cmc1.png")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(cluster1, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
abline(h=0.1,col="red") 
dev.off()

groups.2 <- cutree(cluster1,2)
table(groups.2)
groups.2[order(groups.2)]

meanIAC=apply(IAC,2,mean)
sdCorr=sd(meanIAC)
numbersd=(meanIAC-mean(meanIAC))/sdCorr
png('~/www/iac_sd_graph_round4_cmc1.png')
plot(numbersd)
abline(h=-3)
dev.off()

sdout=-3
outliers=dimnames(datrest2)[[2]][numbersd<sdout]
outliers


indremove=c("MSSM_RNA_PFC_329","MSSM_RNA_PFC_139", "MSSM_RNA_PFC_133", "MSSM_RNA_PFC_163", "MSSM_RNA_PFC_299", "MSSM_RNA_PFC_60", "PENN_RNA_PFC_56", "PENN_RNA_PFC_71","MSSM_RNA_PFC_304","MSSM_RNA_PFC_348","MSSM_RNA_PFC_89")
tmm.cpm.expr = tmm.cpm.expr[,!colnames(tmm.cpm.expr) %in% indremove]
demos.sort = demos.sort[!demos.sort$RNAseq.Sample_RNA_ID %in% indremove,]

### Winsorize counts (Gene outliers in samples)
# Set gene counts in specific samples that are deviating 3 sd from other samples to 3 SD limit
y = y[,!colnames(y) %in% indremove]

raw.wins.mat = apply(y$counts, 1, function(x){
  mn = mean(x, na.rm = T)
  std.dev = sd(x, na.rm = T)
  x[x < (mn-3*std.dev)] <- (mn-3*std.dev)
  x[x > (mn+3*std.dev)] <- (mn+3*std.dev)
  return(x)
}) %>% t

y$counts = raw.wins.mat




################################################################
## correlation between cibersortx proportions

index <- colnames(y$counts)
master.sub <- master2[master2$RNAseq.Sample_RNA_ID %in% index,]
mat <- master.sub[,colnames(master.sub) %in% c("Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8","In1","In2","In3","In4","In5","In6","In8","Astrocytes","Endothelial","Microglia","Neurons","OPC","Oligodendrocytes")]

form <- ~ Ex1+	Ex2+	Ex3+	Ex4+	Ex5+	Ex8+	In1+	In2+	In3+	In4+	In5+	In6+	In8+ Neurons + Astrocytes + Endothelial + Microglia + OPC + Oligodendrocytes
# compute canonical correlation analysis between all pairs of variables, output = rho/ sum(rho) (range of values = 0 to 1)
# In7,Ex7 all zero
# Ex6 = zero variance
C = canCorPairs(form, mat)

pdf("~/www/allcovar_cca_cmc.pdf")
plotCorrMatrix(C)
dev.off()

# remove collinear variables, repeat
C <- as.data.frame(C)
write.table(C,"/sc/arion/projects/psychgen/alanna/cibersort/cmc_corrs.txt", quote=FALSE)






################################################################
## proportion of GE variance explained by cibersortx proportions

# specify variables to be included in voom() estimates of uncertainty
design <- model.matrix(~ Ex1 + Ex2, mat)
vobjGenes <- voom(y, design)

form <- ~ Ex1+	Ex2+	Ex3+	Ex4+	Ex8+	In1+	In2+	In3+	In4+	In5+	In6+	In8+ Neurons + Astrocytes + Endothelial + Microglia + OPC + Oligodendrocytes

cl <- makeCluster(8)
registerDoParallel(cl)

varPart <- fitExtractVarPartModel(vobjGenes, form, mat)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )


# Bar plot of variance fractions for the first 10 genes
pdf("~/www/barplot_variancefractions_first10genes_unadjustedGE.pdf")
plotPercentBars( vp[1:10,] )
dev.off()

# violin plot of contribution of each variable to total variance
pdf("~/www/violin_varcontribution_unadjustedGE.pdf")
plotVarPart( vp )
dev.off()


write.table(varPart, "cmc_varpart_output.txt", quote=FALSE)





################################
## proportion of GE variance explained by cibersortx proportions (repeat while including other covariates)

form <- ~ Ex1+	Ex2+	Ex3+	Ex4+	Ex8+	In1+	In2+	In3+	In4+	In5+	In6+	In8+ Neurons + Astrocytes + Endothelial + Microglia + OPC + Oligodendrocytes + (1|Dx) + (1|Sex) + Age_of_Death + (1|Institution) + PMI_.in_hours. + rnaSeq_isolation.RIN + RIN2 + EV.1 + EV.2 + EV.3 + EV.4 + EV.5 + (1|clusterLIB)

master2 <- merge(master.sub,demos.sort,by="RNAseq.Sample_RNA_ID")

# specify variables to be included in voom() estimates of uncertainty
design <- model.matrix(~ Ex1 + Ex2, master2)
vobjGenes <- voom(y, design)

cl <- makeCluster(8)
registerDoParallel(cl)

varPart <- fitExtractVarPartModel(vobjGenes, form, master2)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )


# Bar plot of variance fractions for the first 10 genes
pdf("~/www/barplot_variancefractions_first10genes_unadjustedGE.pdf")
plotPercentBars( vp[1:10,] )
dev.off()

# violin plot of contribution of each variable to total variance
pdf("~/www/violin_varcontribution_unadjustedGE.pdf",height=10,width=20)
plotVarPart( vp )
dev.off()


write.table(varPart, "cmc_varpart_output_cells+techcovs.txt", quote=FALSE)








################################
## repeat above, test only 1) cell populations with estimates >0.05% in >5% of samples


form <- ~ Ex1 + Ex2 + Ex3 + Ex4 + Ex8 + In1 + In2 + In6 + In8 + Astrocytes + Endothelial + Microglia + Neurons + OPC + Oligodendrocytes + (1|Dx) + (1|Sex) + Age_of_Death + (1|Institution) + PMI_.in_hours. + rnaSeq_isolation.RIN + RIN2 + EV.1 + EV.2 + EV.3 + EV.4 + EV.5 + (1|clusterLIB)

master2 <- merge(master.sub,demos.sort,by="RNAseq.Sample_RNA_ID")

keep.cell = colSums(cmc>0.0005) >= 0.05*nrow(cmc)
# based on above exclude Ex6, Ex7, In3, In4, In5, In7


# specify variables to be included in voom() estimates of uncertainty
design <- model.matrix(~ Ex1 + Ex2, master2)
vobjGenes <- voom(y, design)

cl <- makeCluster(8)
registerDoParallel(cl)

varPart <- fitExtractVarPartModel(vobjGenes, form, master2)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )


# Bar plot of variance fractions for the first 10 genes
pdf("~/www/barplot_variancefractions_first10genes_unadjustedGE.pdf")
plotPercentBars( vp[1:10,] )
dev.off()

# violin plot of contribution of each variable to total variance
pdf("~/www/violin_varcontribution_unadjustedGE.pdf",height=10,width=20)
plotVarPart( vp )
dev.off()


write.table(varPart, "cmc_varpart_output_cells+techcovs.txt", quote=FALSE)





###########################
### adjust for all covs (tech + cells)

### Normalisation

## keep covariates that explain >/=1% of expr variation in >/=10% of genes
# FALSE = Dx, Sex, Ex4, Ex8, In2, In8, Microglia, PMI, All Ancestry PCs
# will just remove cell props that don't meet this criteria

colSums(varPart>=0.01) >= 0.1*nrow(varPart)

# design matrix
design.adj <- model.matrix(~0+ Ex1 + Ex2 + Ex3 + In1 + In6 + Astrocytes + Endothelial + Neurons + OPC + Oligodendrocytes+ Dx + Sex + Age_of_Death + Institution + PMI_.in_hours. + rnaSeq_isolation.RIN + RIN2 + clusterLIB + EV.1 + EV.2 + EV.3 + EV.4 + EV.5, master2)

# Estimate voom weights
#tmp = y$counts
#tmp[is.na(tmp)] = 0
VOOM.GENE_EXPRESSION = voom(y, design=design.adj, plot=F)

# Fit linear model using new weights and new design
ADJUSTED.FIT = lmFit(VOOM.GENE_EXPRESSION)
  
# Residuals after normalisation
RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(ADJUSTED.FIT, VOOM.GENE_EXPRESSION$E)
write.table(RESIDUAL.GENE_EXPRESSION, "knowncovar_+cellprop_adj_quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt")




