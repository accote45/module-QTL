
args<-commandArgs(trailingOnly=TRUE)
module=args[1]
chr=args[2]

### Find QTLs of PC's (using matrixEQTL)
#EV.1 + EV.2 + EV.3 + EV.4 + EV.5

library(MatrixEQTL)
library(tidyverse)
library(data.table)
library(ggplot2)
library(factoextra)

## Location of the package with the data files.
## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS


# Gene expression file name
### Subset GE matrix for gene set of interest, run PCA
ge.mat <- read.table("knowncovar_+cellprop_adj_quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt", header=TRUE, stringsAsFactors=FALSE)
dat <- read.table(paste0("/megena/",module),header=T,stringsAsFactors=F)
geneset <- dat[,1]
geneset.mat <- subset(ge.mat, rownames(ge.mat) %in% geneset)

index <- names(geneset.mat)


PC <- prcomp(t(geneset.mat), scale. = TRUE)


pdf(paste0("PCsofExpression_screeplot_",module,".pdf"))
fviz_eig(PC)
dev.off()

pdf(paste0("PCsofExpression_pcaplot_",module,".pdf"))
fviz_pca_ind(PC,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
             )
dev.off()

pdf(paste0("PCsofExpression_loadingsplot_",module,".pdf"))
fviz_pca_var(PC,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
             )
dev.off()


## get the name of the top 10 measurements (genes) that contribute most to pc1
loading_scores <- PC$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

write.table(PC$rotation[top_10_genes,1], paste0(module,".top10genesforpc1_cmc1"),quote=FALSE) ## show the scores (and +/- sign)


 ## get the name of the top 10 measurements (genes) that contribute most to pc2
loading_scores <- PC$rotation[,2]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_100_genes <- names(gene_score_ranked[1:10])

write.table(PC$rotation[top_100_genes,2], paste0(module,".top10genesforpc2_cmc1"),quote=FALSE) ## show the scores (and +/- sign)


expr <- t(PC$x[,1])
write.table(expr, paste0("SigPCsofExpression_",module,"_cmc1"),quote=FALSE)

pvOutputThreshold = 1
#errorCovariance = numeric();

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = " ";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(paste0("SigPCsofExpression_",module,"_cmc1"));


geno <- read.table(paste0("chr",chr,".DOSonly.fin.fin.vcf.gz"), header=TRUE)
geno <- geno[,-c(1:2,4:5)]
key <- read.csv("Release3_SampleID_key_metadata.csv.gz")
cols <- key[match(colnames(geno), key[["Genotypes.Genotyping_Sample_ID"]], nomatch=0), 'RNAseq.Sample_RNA_ID']
cols2 <- c("SNP",as.character(cols))

colnames(geno) <- cols2
rownames(geno) <- geno$SNP
geno$SNP <- NULL

index <- colnames(gene)

geno2 <- geno[,colnames(geno) %in% index]


## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = " ";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
#snps$columnNames = colnames()
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$CreateFromMatrix(as.matrix(geno2));

## Load covariates

demos <- read.csv('demos.forMatrixQTL.anc.transposed.cmc1.confetiadj.csv')
demos <- demos[,colnames(demos) %in% index]

cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$CreateFromMatrix(as.matrix(demos));

me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = paste0(chr,"_",module,"_coexpr_qtl_cmc1_EA_allchr.txt"),
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE)



