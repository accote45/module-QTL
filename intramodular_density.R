


library(patchwork)
library(data.table)
library(RColorBrewer)
library(plyr)
library(WGCNA)
library(matrixStats)

# calculate intramodular density
# intramodular connectivity  = sum across all genes in module (sum of pairwise abs pearson correlation for gene with every other gene in module) / (modsize*(modsize-1))
# average across genes within module

# read in gene list files for each module
files2 <- list.files(pattern = "*txt", recursive = TRUE, full.names=TRUE)
files2 <- files2[!grepl(c("genes_formods_testable_in_rosmap.txt.gz"),files2)]
files2 <- files2[!grepl(c("genes.txt.gz"),files2)]
files2 <- files2[!grepl(c("master_module_file.txt"),files2)]
files2 <- files2[!grepl(c("module.table.txt.gz"),files2)]

expr <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/knowncovar_+cellprop_adj_quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt.gz", header=TRUE, stringsAsFactors=FALSE)

density <- list()

	# read in module files
	listOfFiles <- lapply(files2, function(x) read.table(x, header = TRUE)) 

	# for each module, cor.mat
	for(j in 1:length(listOfFiles)){
		expr.sub <- expr[rownames(expr) %in% listOfFiles[[j]][,1],]
		cor.mat <- cor(t(expr.sub))
		cor.mat[cor.mat=="1"] <- NA
		density[j] = sum(rowSums(abs(as.matrix(cor.mat)),na.rm = TRUE))/(nrow(cor.mat)*(nrow(cor.mat)-1))
	}
		names(density) <- files2

density.dat <- ldply(density, rbind)

density.dat <- separate(density.dat, .id, into=c("type","module"),sep="/")
density.dat$module <- gsub('.gz','',density.dat$module)
colnames(density.dat) <- c("type","module","density")


