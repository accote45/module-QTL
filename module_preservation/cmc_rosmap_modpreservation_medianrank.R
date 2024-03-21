

##########################################
### CMC module preservation in ROSMAP

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
ge.mat <- read.table("knowncovar_+cellprop_adj_outlierrem_winsorized_expression_rosmap_euro.txt.gz")
rosmap.dat <- as.data.frame(t(ge.mat))


anc <- read.table("rosmap_ancestry_determination.txt", header=TRUE)
imp_ids <- read.table("ROSMAP_arrayGenotype.fam", header=FALSE)
imp_ids$sample <- imp_ids$V2
merged <- merge(imp_ids, anc, by="sample")
merged$final <- paste0(merged$V1,"_",merged$sample)
merged <- merged[!is.na(merged$assignment),]
key <- read.csv("ROSMAP_IDkey.csv") # gwas_id, rnaseq_id, etc
key$sample <- key$gwas_id
master <- merge(merged,key,by="sample")
master$rnaseq_id_fin <- paste0("X", master$rnaseq_id)
rosmap.fin <- rosmap.dat[rownames(rosmap.dat) %in% master$rnaseq_id_fin,]
geno <- fread(paste0("chr22.doseonly.vcf"),header=TRUE)
geno$ID <- paste0(geno$"#CHROM",":",geno$POS)
geno <- as.data.frame(geno) %>% column_to_rownames("ID")
geno <- geno[,-c(1:4)]
geno <- geno[,colnames(geno) %in% master$wgs_id]
master <- master %>% filter(!is.na(master$wgs_id) & wgs_id != "")
master <- master %>% filter(!is.na(master$rnaseq_id_fin) & rnaseq_id_fin != "")
master <- distinct(master[,colnames(master) %in% c("sample","rnaseq_id_fin","wgs_id","final")])
# at this point found samples with single rnaseq ID but multiple wgs IDs, excluded these samples
master <- master[!(master$rnaseq_id_fin %in% c("X08_120410","X201_120424","X216_120425","X542_120516","X637_120524")),]
colnames(geno) <- dplyr::recode(
  colnames(geno),
  !!!setNames(master$rnaseq_id_fin, master$wgs_id)
)
geno.fin <- geno[,colnames(geno) %in% rownames(rosmap.fin)]
geno.fin <- geno.fin[,order(names(geno.fin))]
geno.fin[] <- lapply(geno.fin, function(x) as.numeric(as.character(x)))
geno.fin <- geno.fin[complete.cases(geno.fin),]
rosmap.fin <- rosmap.fin[rownames(rosmap.fin) %in% colnames(geno.fin),]


# load module labels, and label for unassigned genes
colorCMC1 = read.table("master_module_file.txt",header=T)
colorCMC1 <- colorCMC1[colorCMC1$gene.names %in% colnames(cmc.dat),]
index <- colnames(cmc.dat)[!(colnames(cmc.dat) %in% colorCMC1$gene.names)]
index2 <- as.data.frame(index)
index2$module <- "unknown"
colnames(index2) <- c("gene.names","module")
colorCMC1 <- rbind(index2,colorCMC1)

list1a <- unstack(colorCMC1)

index <- intersect(colnames(cmc.dat),colnames(rosmap.fin))

rosmap.master <- rosmap.fin[,colnames(rosmap.fin) %in% index]
cmc.master <- cmc.dat[,colnames(cmc.dat) %in% index]

# number of networks used in the consensus network analysis:
nSets = 2
# Vector with descriptive names of the two sets.
setLabels = c("CMC", "ROSMAP")
shortLabels = c("CMC", "ROSMAP")
# Define a list whose components contain the data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = cmc.master)
names(multiExpr[[1]]$data) = names(cmc.master)
rownames(multiExpr[[1]]$data) = dimnames(cmc.master)[[1]]
multiExpr[[2]] = list(data = rosmap.master)
names(multiExpr[[2]]$data) = names(rosmap.master)
rownames(multiExpr[[2]]$data) = dimnames(rosmap.master)[[1]]
# Check that the data has the correct format:
exprSize = checkSets(multiExpr)
names(multiExpr) <- c("CMC","ROSMAP")

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
nPermutations = 3,
verbose = 3)


medianrank = mp$preservation$observed$ref.CMC$inColumnsAlsoPresentIn.ROSMAP
logp = mp$preservation$log.p$ref.CMC$inColumnsAlsoPresentIn.ROSMAP

write.table(medianrank, file=paste0("MedianRank.ControlRef.vROSMAP.",PathwayName,".txt"), sep="\t",quote=F)
write.table(logp, file=paste0("LogP.ControlRef.vROSMAP.",PathwayName,".txt"), sep="\t",quote=F)

printFlush(paste("\n Finished (CMC=reference network) for pathway", PathwayName, "\n"));


}

}
