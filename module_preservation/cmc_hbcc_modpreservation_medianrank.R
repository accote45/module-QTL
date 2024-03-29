
#### Module preservation analysis
# reference : http://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html

library(WGCNA)
library(data.table)
library(tidyverse)
library(dplyr)
library(plyr)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

ge.mat <- read.table("knowncovar_+cellprop_adj_quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt.gz")
cmc.dat <- as.data.frame(t(ge.mat))
ge.mat <- read.table("knowncovar_+cellprop_adj_quantnorm_outlierrem_winsorized_expression_hbcc_euro.txt.gz", header=TRUE, stringsAsFactors=FALSE)
hbcc.dat <- as.data.frame(t(ge.mat))

geno <- read.table(paste0("chr22.vcf"), header=TRUE)
geno <- geno[,-c(2:5)]
key2 <- read.csv('demos.master.EA.amanda.ancestry.csv.gz')
key2$ID <- paste0("X",key2$ID)

cols <- key2[match(colnames(geno), key2[["ID"]], nomatch=0), 'RNAseq.Sample_RNA_ID']
cols2 <- c("SNP",as.character(cols))
colnames(geno) <- cols2
rownames(geno) <- geno$SNP
geno$SNP <- NULL
hbcc.fin <- hbcc.dat[rownames(hbcc.dat) %in% colnames(geno),]


# load module labels, and label for unassigned genes
colorCMC1 = read.table("master_module_file.txt",header=T)
colorCMC1 <- colorCMC1[colorCMC1$gene.names %in% colnames(cmc.dat),]
index <- colnames(cmc.dat)[!(colnames(cmc.dat) %in% colorCMC1$gene.names)]
index2 <- as.data.frame(index)
index2$module <- "unknown"
colnames(index2) <- c("gene.names","module")
colorCMC1 <- rbind(index2,colorCMC1)


list1a <- unstack(colorCMC1)

index <- intersect(colnames(cmc.dat),colnames(hbcc.fin))

hbcc.master <- hbcc.fin[,colnames(hbcc.fin) %in% index]
cmc.master <- cmc.dat[,colnames(cmc.dat) %in% index]

# number of networks used in the consensus network analysis:
nSets = 2
# Vector with descriptive names of the two sets.
setLabels = c("CMC", "HBCC")
shortLabels = c("CMC", "HBCC")
# Define a list whose components contain the data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = cmc.master)
names(multiExpr[[1]]$data) = names(cmc.master)
rownames(multiExpr[[1]]$data) = dimnames(cmc.master)[[1]]
multiExpr[[2]] = list(data = hbcc.master)
names(multiExpr[[2]]$data) = names(hbcc.master)
rownames(multiExpr[[2]]$data) = dimnames(hbcc.master)[[1]]
# Check that the data has the correct format:
exprSize = checkSets(multiExpr)
names(multiExpr) <- c("CMC","HBCC")

Genes=colnames(cmc.master)
nGenes=length(Genes)


for(i in 1:length(list1a))
{
PathwayName=names(list1a[i])
printFlush("=====================================");
printFlush(paste("Now working on pathway =>", PathwayName, "\n"));

GeneSet <- list1a[i]
NumberGenes=length(GeneSet[[1]])
Convergent=  GeneSet[[1]][GeneSet[[1]] %in% Genes]
Divergent= NumberGenes-length(Convergent);
    
printFlush("MATCHING GENES:",length(Convergent), "      FALSE GENES:", Divergent);
printFlush("=====================================\n");
#doModulePreservation = TRUE;
#doClusterRepro = TRUE;
common=length(Convergent)

# >50% of genes must be shared btw primary and secondary dataset
minimumoverlap=0.50*NumberGenes

#Set
if (common > minimumoverlap) {
temp <- colorCMC1
temp$module[!(temp$module==PathwayName)] <- "0"
temp <- distinct(temp[temp$gene.names %in% colnames(cmc.master),])
temp <- temp[order(temp$module,decreasing=T),]
temp <- temp %>% arrange(gene.names) %>% filter(!duplicated(gene.names))
multiColor = list(CMC = temp$module)
names(multiColor[[1]]) = temp$gene.names

# preservation stats calc
mp = modulePreservation(multiExpr, multiColor, dataIsExpr=T,networkType ="unsigned", corFnc = "cor",
referenceNetworks = 1, maxModuleSize =common, maxGoldModuleSize=common,
loadPermutedStatistics = FALSE,
nPermutations = 100,
verbose = 3)


medianrank = mp$preservation$observed$ref.CMC$inColumnsAlsoPresentIn.HBCC
logp = mp$preservation$log.p$ref.CMC$inColumnsAlsoPresentIn.HBCC

write.table(medianrank, file=paste0("MedianRank.ControlRef.vHBCC.",PathwayName,".txt"), sep="\t",quote=F)
write.table(logp, file=paste0("LogP.ControlRef.vHBCC.",PathwayName,".txt"), sep="\t",quote=F)
printFlush(paste("\n Finished (CMC=reference network) for pathway", PathwayName, "\n"));


}

}

