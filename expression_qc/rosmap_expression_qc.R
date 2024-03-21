

####################################################################################
### read in ROSMAP expr, complete preliminary QC


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
library(dplyr)
library(sva)
library(pamr)

# read in major files

expr <- read.table("/sc/arion/projects/psychgen/DATA/ROSMAP_psychgen/RNA_seq/raw/ROSMAP_all_counts_matrix.txt.gz", header=TRUE)
meta.clinical <- read.csv("/sc/arion/projects/psychgen/DATA/ROSMAP_psychgen/metadata/ROSMAP_Clinical_2019-05_v3.csv") #projid, individualID
meta.rna <- read.csv("/sc/arion/projects/psychgen/DATA/ROSMAP_psychgen/metadata/ROSMAP_assay_RNAseq_metadata.csv") #specimenID
meta.bio <- read.csv("/sc/arion/projects/psychgen/DATA/ROSMAP_psychgen/metadata/ROSMAP_biospecimen_metadata.csv") #individualID ,specimenID
key <- read.csv("/sc/arion/projects/psychgen/DATA/ROSMAP_psychgen/metadata/ROSMAP_IDkey.csv") # gwas_id, rnaseq_id, etc

key$specimenID <- key$rnaseq_id

meta.temp <- join(meta.clinical, meta.bio, by="individualID")
meta.master.temp <- join(meta.temp, meta.rna, by="specimenID") # specimenID, individualID, projid, 

expr$ensembl_gene_id <- sub("\\..*", "", expr$feature)


mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                       filters = "ensembl_gene_id", values = expr$ensembl_gene_id,
                       mart = mart, uniqueRows = TRUE)

# filter for only autosomes
expr = join(expr, Ensemble2HGNC, by="ensembl_gene_id", type="left") %>% filter(chromosome_name != "X" & chromosome_name != "Y" & chromosome_name !="MT")
expr <- expr %>% column_to_rownames("ensembl_gene_id")


# subset to only bulk dlpfc samples
temp2 <- meta.master.temp %>% filter( grepl("dorsolateral",tissue))
meta.master <- temp2 %>% filter(grepl("bulk",nucleicAcidSource))

# subset to only sampel of genotype-derived european ancestry
anc <- read.table("/sc/arion/projects/psychgen/DATA/ROSMAP_psychgen2/genotype/ancestry/rosmap_ancestry_determination.txt", header=TRUE)
anc$gwas_id <- anc$sample
anc.merge <- join(anc, key, by="gwas_id")
anc.merge$specimenID <- anc.merge$rnaseq_id
meta.master<- join(meta.master, anc.merge, by='specimenID')
meta.master <- meta.master[!is.na(meta.master$assignment),]

# Add harmonised case-control status
meta.master$Diagnosis = 'OTHER'
meta.master$Diagnosis[meta.master$cogdx == 1 & meta.master$braaksc <= 3 & meta.master$ceradsc >= 3] = 'CONTROL'
meta.master$Diagnosis[meta.master$cogdx == 4 & meta.master$braaksc >= 4 & meta.master$ceradsc <= 2] = 'AD'

# Keep only individuals with info in demos and expr file
covars <- c("msex", "age_death", "RIN", "libraryBatch", "pmi", "Diagnosis", "cogdx", "dcfdx_lv", "specimenID", "individualID")
demos.sub <- meta.master[,covars]

# keep demos rows with RIN, age of death, pmi info
demos.sub <- demos.sub[!is.na(demos.sub$RIN),]
demos.sub <- demos.sub[!is.na(demos.sub$pmi),]
demos.sub <- demos.sub[!is.na(demos.sub$age_death),]
demos.sub <- distinct(demos.sub)

index <- c(paste0("X",demos.sub$specimenID))
count.trim <- expr[,(names(expr) %in% index)]
both <- substring(intersect(index, names(count.trim)), 2)

demos.fin <- demos.sub[(demos.sub$specimenID %in% both),]
expr.fin <- count.trim[,(names(count.trim) %in% paste0("X", both))]


# sort files, make DGEList object
demos.sort <- demos.fin[order(demos.fin$specimenID),] 
count.sort <- expr.fin[ , order(names(expr.fin))]

### Normalisation

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
plotdata <- data.frame(specimenID=substring(rownames(PC$x), 2), 
                       PC1=PC$x[,1], 
                       PC2=PC$x[,2])
plotdata <- merge(plotdata, demos.sort, by="specimenID", all.y=TRUE)

plotdata$age_death_fin <- gsub("\\+", "", as.character(factor(plotdata$age_death)))

pdf("~/www/PCAplot_rosmap_euro_dlpfc_unadjusted.pdf")
p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(colour=factor(libraryBatch), fill=factor(libraryBatch), shape=factor(cogdx), size=as.numeric(age_death_fin))) 
p <- p + scale_shape_manual(values=c(3,22,21,24,5,9))
p <- p + theme_bw() %+replace% theme(legend.position="right")
p
dev.off()
# For ROSMAP cohort , remove one sample based on inspection of pc plot, specimen ID = 180_120424
plotdata <- plotdata[order(plotdata$PC2),]
tail(plotdata)

indremove <- c("X380_120503")
y = y[,!colnames(y) %in% indremove]
demos.sort = demos.sort[!demos.sort$specimenID %in% c("380_120503"),]

# check for outlier samples by inter array correlation (samples with IAC < 3SDs below mean IAC for the dataset will be removed) ### source = https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/HumanBrainTranscriptome/Identification%20and%20Removal%20of%20Outlier%20Samples%20-%20Illumina.pdf
IAC=cor(tmm.cpm.expr,use="p") 
library(WGCNA)

png("~/www/iac_histogram_rosmap.png")
hist(IAC,sub=paste("mean=",format(mean(IAC[upper.tri(IAC)]),digits=3))) 
dev.off()

sampleTree = hclust(as.dist(1-IAC), method = "average");
png("~/www/iacbased_outlier_check_rosmap.png")
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
png('~/www/iac_sd_graph_rosmap.png')
plot(numbersd)
abline(h=-3)
dev.off()

sdout=-3
outliers=dimnames(tmm.cpm.expr)[[2]][numbersd<sdout]
outliers

removevec=c("X246_120426","X367_120502","X380_120503","X956_131107")
overlap1=is.element(dimnames(tmm.cpm.expr)[[2]],removevec)
datrest2=tmm.cpm.expr[,!overlap1]
dim(datrest2)

indremove=c("X246_120426","X367_120502","X380_120503","X956_131107")
tmm.cpm.expr = tmm.cpm.expr[,!colnames(tmm.cpm.expr) %in% indremove]
index <- gsub("X","",colnames(tmm.cpm.expr))
demos.sort = demos.sort[demos.sort$specimenID %in% index,]

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

# at this point found duplicated specimenID in demos.sort file (specimenID = 492_120515), excluded sample

demos.sort <- demos.sort[!(demos.sort$specimenID=="492_120515"),]
demos.sort$age_death_fin <- as.numeric(gsub("\\+", "", as.character(factor(demos.sort$age_death))))
demos.sort$Diagnosis <- as.character(demos.sort$Diagnosis)
demos.sort$msex <- as.character(demos.sort$msex)
demos.sort$libraryBatch <- as.character(demos.sort$libraryBatch)
demos.sort$cogdx <- as.character(demos.sort$cogdx)
demos.sort$dcfdx_lv <- as.character(demos.sort$dcfdx_lv)




################################################################
## correlation between cibersortx proportions

index <-paste0("X",demos.sort$specimenID)

master <- read.csv("/sc/arion/projects/psychgen/alanna/cibersort/rosmap/CIBERSORTx_Job15_Results.csv")
map <- read.table("/sc/arion/projects/psychgen/DATA/ROSMAP_psychgen/metadata/hannah_map_080420", fill=TRUE,header=TRUE)
map.sub <- map[map$specimenID %in% demos.sort$specimenID,]

master.sub <- master[master$Mixture %in% paste0("X",map.sub$RNA_IDs),]
mat <- master.sub[,colnames(master.sub) %in% c("Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8","In1","In2","In3","In4","In5","In6","In8","Astrocytes","Endothelial","Microglia","Neurons","OPC","Oligodendrocytes")]

form <- ~ Ex1+	Ex2+	Ex3+	Ex4+	Ex5+	Ex8+	In1+	In2+	In3+	In4+	In5+	In6+	In8+ Neurons + Astrocytes + Endothelial + Microglia + OPC + Oligodendrocytes
# compute canonical correlation analysis between all pairs of variables, output = rho/ sum(rho) (range of values = 0 to 1)
# In7,Ex7 all zero
# Ex6 = zero variance
C = canCorPairs(form, mat)

pdf("~/www/allcovar_cca_rosmap.pdf")
plotCorrMatrix(C)
dev.off()

# remove collinear variables, repeat
C <- as.data.frame(C)
write.table(C,"/sc/arion/projects/psychgen/alanna/cibersort/rosmap/rosmap_corrs.txt", quote=FALSE)





################################################################
## proportion of GE variance explained by cibersortx proportions

# specify variables to be included in voom() estimates of uncertainty
y = y[,colnames(y) %in% paste0("X",map.sub$specimenID)]

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


write.table(varPart, "/sc/arion/projects/psychgen/alanna/cibersort/rosmap/rosmap_varpart_output.txt", quote=FALSE)


################################
## proportion of GE variance explained by cibersortx proportions (repeat while including other covariates)

form <- ~ Ex1+	Ex2+	Ex3+	Ex4+	Ex8+	In1+	In2+	In3+	In4+	In5+	In6+	In8+ Neurons + Astrocytes + Endothelial + Microglia + OPC + Oligodendrocytes + (1|msex) + age_death_fin + (1|Diagnosis) + (1|cogdx) + (1|dcfdx_lv) + pmi + RIN + (1|libraryBatch)

temp <- merge(demos.sort,map,by="specimenID")
master.sub$RNA_IDs <- gsub("X","",master.sub$Mixture)
master2 <- merge(master.sub,temp,by="RNA_IDs")

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


write.table(varPart, "/sc/arion/projects/psychgen/alanna/cibersort/rosmap/rosmap_varpart_output_cells+techcovs.txt", quote=FALSE)







################################
## repeat above, test only 1) cell populations with estimates >0.05% in >5% of samples


form <- ~ Ex1+	Ex2+	Ex3+	Ex4+	Ex5 + Ex8+	In1+	In2+	In4+	In5+	In6+	In8+ Neurons + Astrocytes + Endothelial + Microglia + OPC + Oligodendrocytes + (1|msex) + age_death_fin + (1|Diagnosis) + (1|cogdx) + (1|dcfdx_lv) + pmi + RIN + (1|libraryBatch)

keep.cell = colSums(master>0.0005) >= 0.05*nrow(master)
# based on above exclude Ex6, Ex7, In3, In4, In7


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


write.table(varPart, "/sc/arion/projects/psychgen/alanna/cibersort/rosmap/rosmap_varpart_output_cells+techcovs.txt", quote=FALSE)







###########################
### adjust for all covs (tech + cells)

### Normalisation (with ancestry adjustment)

anc <- read.table("/sc/arion/projects/psychgen/DATA/ROSMAP_psychgen2/genotype/ancestry/ROSMAP_gemtools_ancestryPCs.txt")
rownames(anc) <- sub('.*\\_', '', rownames(anc))
anc <- anc %>% rownames_to_column("gwas_id")

key <- read.csv("/sc/arion/projects/psychgen/DATA/ROSMAP_psychgen2/metadata/ROSMAP_IDkey.csv") #projid, rnaseq_id
key.sub <- key[,c(2,5)]
key.sub$specimenID <- key.sub$rnaseq_id

demos.fin <- join(master2, key.sub, by="specimenID")
demos.fin.fin <- join(demos.fin, anc, by="gwas_id")
demos.fin.fin <- demos.fin.fin[!is.na(demos.fin.fin$EV.0),]
demos.fin.fin <- distinct(demos.fin.fin)
index <- paste0("X",demos.fin.fin$specimenID)

y <- y[,colnames(y) %in% index]


## keep covariates that explain >/=1% of expr variation in >/=10% of genes
# FALSE = cogdx, dcfdx_lv, Diagnosis, msex, Ex2, Ex8, In4, In5, OPC, age_daeth_fin, pmi
# will just remove cell props that don't meet this criteria

colSums(varPart>=0.01) >= 0.1*nrow(varPart)

# design matrix
design.adj <- model.matrix(~0+ Ex1+	Ex3+	Ex4+	Ex5+	In1+	In2+	In6+	In8+ Neurons + Astrocytes + Endothelial + Microglia + Oligodendrocytes + msex + age_death_fin + Diagnosis + pmi + RIN + libraryBatch + EV.1 + EV.2 + EV.3 + EV.4, demos.fin.fin)

# Estimate voom weights
#tmp = y$counts
#tmp[is.na(tmp)] = 0
VOOM.GENE_EXPRESSION = voom(y, design=design.adj, plot=F)

# Fit linear model using new weights and new design
ADJUSTED.FIT = lmFit(VOOM.GENE_EXPRESSION)
  
# Residuals after normalisation
RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(ADJUSTED.FIT, VOOM.GENE_EXPRESSION$E)
write.table(RESIDUAL.GENE_EXPRESSION, "/sc/arion/projects/psychgen/alanna/cibersort/rosmap/knowncovar_+cellprop_adj_outlierrem_winsorized_expression_rosmap_euro.txt")


