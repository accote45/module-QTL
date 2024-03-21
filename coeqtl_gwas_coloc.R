
### example script - coeQTL x GWAS colocalization test using coloc


library(coloc)
args<-commandArgs(trailingOnly=TRUE)
module=args[2]
split=args[1]

gwas <- read.table(paste0("gwas_finemap_",split,"_",module,".z_wchrpos"), header=FALSE)
coeqtl <- read.table(paste0("for_coloc_input_",split,"_",module,".z_wchrpos"), header=FALSE)

gwas$beta.gwas <- gwas$V5
gwas$var.gwas <- (gwas$V6)^2
gwas$SNP <- gwas$V13

coeqtl$beta.qtl <- coeqtl$V7
coeqtl$var.qtl <- (coeqtl$V8)^2
coeqtl$SNP <- coeqtl$V9
coeqtl$maf <- coeqtl$V6

pc.sd <- read.table(paste0("/sc/arion/projects/psychgen/alanna/cibersort/cmc/coeqtl/",module,".txtsigpc_sd"), header=FALSE)
names(pc.sd) <- as.matrix(pc.sd[1,])
pc.sd <- as.data.frame(pc.sd)
pc.sd <- pc.sd[-1, , drop=FALSE]

master <- merge(coeqtl, gwas, by="SNP")
# remove duplicate snps
master <- master[!duplicated(master$SNP),]

dataset1=list(beta=master$beta.qtl, varbeta=master$var.qtl, N=488, type="quant", sdY=as.numeric(as.character(pc.sd[[1]])), snp=as.character(master$SNP))
dataset2=list(beta=master$beta.gwas, varbeta=master$var.gwas,N=mean(master$V8.y),type="quant", snp=as.character(master$SNP))

my.res <- coloc.abf(dataset1, dataset2, MAF=master$maf)


