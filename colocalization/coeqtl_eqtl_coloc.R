library(coloc)
library(tidyverse)
library(data.table)

args<-commandArgs(trailingOnly=TRUE)
locus.name=args[1]
eqtl.name=args[2]

temp <- strsplit(eqtl.name,split="_",fixed=TRUE)
module=paste0(temp[[1]][5],"_",temp[[1]][6])
gene=temp[[1]][1]

eqtl <- read.table(paste0(eqtl.name,".z"), header=TRUE)
coeqtl <- read.table(paste0(locus.name,"_for_coloc.txt"), header=FALSE)
colnames(coeqtl) <- c('chr','pos','A1','A2','module','pc','beta','t.stat','p.value','FDR','MAF')

eqtl$beta.eqtl <- eqtl$beta
eqtl$var.eqtl<- (eqtl$beta/eqtl$t.stat)^2
eqtl$SNP <- eqtl$SNP

coeqtl$beta.coeqtl <- coeqtl$beta
coeqtl$var.coeqtl <- (coeqtl$beta/coeqtl$t.stat)^2
coeqtl$SNP <- paste0(coeqtl$chr,":",coeqtl$pos,":",coeqtl$A1,":",coeqtl$A2)

pc.sd <- read.table(paste0(module,".txtsigpc_sd.gz"), header=FALSE)
names(pc.sd) <- as.matrix(pc.sd[1,])
pc.sd <- as.data.frame(pc.sd)
pc.sd <- pc.sd[-1, , drop=FALSE]

master <- merge(coeqtl, eqtl, by="SNP")
# remove duplicate snps
master <- master[!duplicated(master$SNP),]

# get variance of eqtl gene
ge.mat <- read.table("knowncovar_+cellprop_adj_quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt.gz", header=TRUE, stringsAsFactors=FALSE)
gene.sd <- sd(ge.mat[rownames(ge.mat) %in% gene,])

dataset1=list(beta=master$beta.coeqtl, varbeta=master$var.coeqtl, n=488, type="quant", sdY=as.numeric(as.character(pc.sd[[1]])), snp=as.character(master$SNP))
dataset2=list(beta=master$beta.eqtl, varbeta=master$var.eqtl, n=488,type="quant",sdY=gene.sd,snp=as.character(master$SNP))

my.res <- coloc.abf(dataset1, dataset2, MAF=master$maf)

write.table(t(my.res[[1]]), paste0("coloc_summary_",eqtl.name), col.names=TRUE,row.names=FALSE, quote=FALSE)
write.table(my.res$results, paste0("coloc_output_",eqtl.name), col.names=TRUE,row.names=FALSE, quote=FALSE)
