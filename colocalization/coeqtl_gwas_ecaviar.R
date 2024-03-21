

#####################################
#####################################
## eCaviar for coeQTL-GWAS colocalization

############

library(tidyverse)
library(data.table)

setwd("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/ecaviar/coloc_enigma_surfarea")

gwas <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/ecaviar/gwas/ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429.hg38.txt",header=T)
gwas <- as.data.frame(gwas)

##### remove dup snps from enigma gwas (keep only snps not indels)
gwas <- gwas[gwas$A1 %in% c('a','t','c','g'),]
gwas <- gwas[gwas$A2 %in% c('a','t','c','g'),]

gwas$A1 <- toupper(gwas$A1)
gwas$A2 <- toupper(gwas$A2)


##### read in coeqtl files
gwas <- as.data.frame(gwas)
gwas$ID <-paste0("chr",gwas$hg38)

filenames <- list.files(path="/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/coeqtl_maf0.05",pattern="*for_coloc.txt$",full.names=T)
ldf <- lapply(filenames,read.table)
names(ldf) <- gsub("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/coeqtl_maf0.05/","",filenames)
ldf.fin <- list()
for (i in 1:length(ldf)) {
	colnames(ldf[[i]]) <- c('chr','pos','A1','A2','module','pc','beta','t.stat','p.val','FDR','MAF')
	ldf[[i]]$ID <- paste0(ldf[[i]]$chr,":",ldf[[i]]$pos)
	ldf.fin[[i]] <- ldf[[i]][ldf[[i]]$ID %in% gwas$ID,]
}

#### match coeQTL and GWAS SNPs by their effect and ref allele
#gwas$SNP <-paste0("chr",gwas$hg38,":",toupper(gwas$A1),":",toupper(gwas$A2))
gwas$SNP <-paste0("chr",gwas$hg38,":",gwas$A1,":",gwas$A2)

master.qtl <- bind_rows(ldf.fin)
master.qtl$SNP <- paste0(master.qtl$ID,":",master.qtl$A1,":",master.qtl$A2)
master.qtl <- distinct(master.qtl[,c('chr','pos','A1','A2','SNP')])
master.qtl$ID <- paste0(master.qtl$chr,":",master.qtl$pos)
master.qtl$A12 <- paste0(master.qtl$A1,"_",master.qtl$A2)
master.qtl$A21 <- paste0(master.qtl$A2,"_",master.qtl$A1)
colnames(master.qtl) <- c('chr.qtl','pos.qtl','A1.qtl','A2.qtl','SNP.qtl','ID','A12.qtl','A21.qtl')

gwas$Z <- gwas$BETA1/gwas$SE
#gwas$Z <- log(gwas$OR)/gwas$SE
gwas$A12 <- paste0(gwas$A1,"_",gwas$A2)

master.temp <- merge(master.qtl,gwas,by="ID")
master.temp$snpgroup <- ifelse(master.temp$A12==master.temp$A12.qtl | master.temp$A12==master.temp$A21.qtl,"keep","scrap")
master.temp <- master.temp[master.temp$snpgroup=="keep",]
master.temp$Z.fin <- ifelse(master.temp$A1.qtl==master.temp$A1,master.temp$Z,-(master.temp$Z))
gwas.fin <- master.temp[,c(colnames(gwas),"Z.fin")]

snps.fin <- master.temp$ID

###### write per locus GWAS and coeQTL files
names(ldf.fin) <- names(ldf)
ldf.fin.fin <- list()
for (i in 1:length(ldf.fin)) {
	ldf.fin.fin[[i]] <- ldf.fin[[i]][ldf.fin[[i]]$ID %in% snps.fin,]
	ldf.fin.fin[[i]] <- ldf.fin.fin[[i]][,c('ID','t.stat')]
	ldf.fin.fin[[i]] <- ldf.fin.fin[[i]][order(ldf.fin.fin[[i]]$ID),]
	write.table(ldf.fin.fin[[i]],paste0(names(ldf.fin[i]),".z"),quote=F,row.names=F,col.names=F)
}

##### create per locus GWAS files, subset for shared SNPs
gwas.df <- list()
gwas.df.fin <- list()
for (i in 1:length(ldf.fin)) {
	gwas.df[[i]] <- gwas.fin[gwas.fin$ID %in% ldf.fin.fin[[i]]$ID,]
	gwas.df.fin[[i]] <- gwas.df[[i]][,c('ID','Z.fin')]
	gwas.df.fin[[i]] <- gwas.df.fin[[i]][order(gwas.df.fin[[i]]$ID),]
}

names(gwas.df.fin) <- names(ldf)
for (i in 1:length(gwas.df.fin)) {
	write.table(gwas.df.fin[[i]],paste0("gwas_",names(ldf[i]),".z"),quote=F,row.names=F,col.names=F)
}





############
# create LD matrices
library(tidyverse)
library(data.table)

setwd('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/ecaviar/coloc_enigma_surfarea_TEST')

# read in dosage
dat <- list()
for (i in 1:22){
	dat[[i]] <- read.table(paste0("../chr",i,"_coeqtl_esnp_dosage_subset"), header=TRUE)
}

dat <- Filter(function(x) dim(x)[1] > 0, dat)
dosage <- bind_rows(dat)
dosage$SNP <- paste0(dosage$CHR,":",dosage$POS)

# read in snp files
filenames <- list.files(pattern="gwas*",full.names=F)
ldf <- lapply(filenames,read.table)
names(ldf) <- filenames

gwas.df <- list()
corrs <- list()
for (i in 1:length(ldf)) {
	index <- ldf[[i]]$V1
	gwas.df[[i]] <- dosage[dosage$SNP %in% index,]
	gwas.df[[i]] <- gwas.df[[i]][order(gwas.df[[i]]$SNP),]
	gwas.df[[i]] <- as.data.frame(t(gwas.df[[i]]))
	gwas.df[[i]] <- gwas.df[[i]][-c(1:2,4:5,nrow(gwas.df[[i]])),]
	names(gwas.df[[i]]) <- as.matrix(gwas.df[[i]][1,])
	gwas.df[[i]]<- gwas.df[[i]][-1,]
	indx <- sapply(gwas.df[[i]], is.character)
	gwas.df[[i]][indx] <- lapply(gwas.df[[i]][indx], function(x) as.numeric(as.character(x)))
	corrs[[i]] <- cor(gwas.df[[i]])
}

names <- gsub("gwas_","",filenames)
names <- gsub(".txt.z","",names)

for (i in 1:length(corrs)){
	write.table(corrs[[i]],paste0(names[[i]],".ld"),quote=FALSE,row.names=FALSE,col.names=FALSE)
}





##############################################################################
## run eCaviar (in each subdirectory)
# CLPP threshold of 0.001 is recommended by creator

ml caviar

for locus in $(ls gwas*.z | sed 's/gwas_//g' | sed 's/.txt.z//g'); do
bsub -J surfarea.${locus} -P acc_psychgen -q premium -R rusage[mem=40000] -W 96:00 -e e.${locus} -o o.${locus} eCAVIAR -l ${locus}.ld -z ${locus}.txt.z -l ${locus}.ld -z gwas_${locus}.txt.z -f 1 -c 3 -o ${locus}.output
done








