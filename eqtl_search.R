

args<-commandArgs(trailingOnly=TRUE)
chr=args[1]


source("Matrix_eQTL_engine.R.gz")


### Find QTLs of PC's (using matrixEQTL)
#EV.1 + EV.2 + EV.3 + EV.4 + EV.5

library(tidyverse)
library(data.table)
library(ggplot2)
library(factoextra)
library(biomaRt)

## Location of the package with the data files.
## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Gene expression file name
### Subset GE matrix for gene set of interest, run PCA
ge.mat <- read.table("knowncovar_+cellprop_adj_quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt.gz", header=TRUE, stringsAsFactors=FALSE)
#ge.mat <- as.data.frame(t(ge.mat))

index <- names(ge.mat)

pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 1

cisDist = 1e6

genepos <- read.table("genepos_cmc.txt.gz", header=TRUE)

gene.fin.fin <- ge.mat[rownames(ge.mat) %in% genepos$geneid,]


## Load gene expression data


gene = SlicedData$new();
gene$fileDelimiter = " ";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$CreateFromMatrix(as.matrix(gene.fin.fin));

geno <- read.table(paste0("chr",chr,".DOSonly.fin.fin.vcf"), header=TRUE)

cmc <- fread("master_autosomes.maf",header=F)
maf <- cmc[cmc$V2>0.05,]
geno <- geno[geno$ID %in% maf$V1,]

geno <- geno[,-c(1:2,4:5)]
key <- read.csv("Release3_SampleID_key_metadata.csv.gz")
cols <- key[match(colnames(geno), key[["Genotypes.Genotyping_Sample_ID"]], nomatch=0), 'RNAseq.Sample_RNA_ID']
cols2 <- c("SNP",as.character(cols))

colnames(geno) <- cols2
rownames(geno) <- geno$SNP
geno$SNP <- NULL

index <- colnames(gene)

geno2 <- geno[,colnames(geno) %in% index]

temp <- as.data.frame(rownames(geno2))
names(temp) <- "V1"
temp2 <- separate(temp, V1, into=c("V2","V3","V4","V5",sep=":"))
snpspos <- data.frame(snp=temp$V1, chr=temp2$V2, pos=temp2$V3)
snpspos$pos <- as.numeric(as.character(snpspos$pos))


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

demos <- read.csv('demos.forMatrixQTL.anc.transposed.cmc1.confetiadj.csv.gz')
index <- colnames(snps)
demos <- demos[,colnames(demos) %in% index]

cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$CreateFromMatrix(as.matrix(demos));

me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = paste0(chr,"_single_eqtl_cmc_EA_trans.txt"),
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    verbose = TRUE,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    output_file_name.cis = paste0(chr,"_single_eqtl_cmc_EA_cis.txt"),
    snpspos=snpspos,
    genepos=genepos,
    cisDist=cisDist,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory=TRUE)
