

######################################################
###### colocalization of coeQTLs and cis/trans eQTLs for module genes

ls *ld | awk -F"_" '{print $0" "$4"_"$5" "$1}' | sed 's/chr//g' | sed 's/_for_coloc.txt.ld//g' | awk '{print "chr"$0}' > job.file.txt

while read coeqtl module chr; do
Rscript --vanilla input_files.R $coeqtl $module $chr
#bsub -J ${coeqtl} -P acc_psychgen -q express -R rusage[mem=40000] -W 12:00 -e e.${coeqtl} -o o.${coeqtl} Rscript --vanilla input_files.R
done < job.file.txt

## create coeQTL and eQTL input files (no need to match effect sizes)
# loop through script for each coeQTL
library(tidyverse)
library(data.table)

args<-commandArgs(trailingOnly=TRUE)
coeqtl=args[1]
module=args[2]
chr=args[3]

setwd('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/ecaviar/coloc_single_eqtl')

coeqtl.dat <- read.table(paste0('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/coeqtl_maf0.05/',coeqtl,'_for_coloc.txt'))
colnames(coeqtl.dat) <- c("CHR","POS","A1","A2","ProbeID","V6","BETA","t.stat","PVAL","FDR","MAF")
coeqtl.dat$SNP <- paste0(coeqtl.dat$CHR,":",coeqtl.dat$POS,":",coeqtl.dat$A1,":",coeqtl.dat$A2)

# find significant cis/trans eQTLs for coeQTL module
eqtl.sig <- read.table('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/single_eqtls/single_gene_eqtl_allmodgenes/master_sig_maf0.05.txt',header=T)
module.genes <- read.table(paste0("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/megena/",module,".txt"),header=T)
eqtl.sig.modgenes <- eqtl.sig[eqtl.sig$gene %in% module.genes[,1],]

# if sig eSNPs for coeQTL and module eQTLs do not overlap, end job
if(length(intersect(coeqtl.dat$SNP,eqtl.sig.modgenes$SNP))==0) {
	stop("Job stopped. No overlapping eSNPs between significant coeQTL and eQTL")
}

# for each eQTL eGene, create eCaviar input file
ciseqtl <- fread(paste0('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/single_eqtls/single_gene_eqtl_allmodgenes/',chr,'_single_eqtl_cmc_EA_cis.txt'))
transeqtl <- fread(paste0('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/single_eqtls/single_gene_eqtl_allmodgenes/',chr,'_single_eqtl_cmc_EA_trans.txt'))
genes <- unique(eqtl.sig.modgenes$gene)

loadings <- read.table(paste0('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/pc_modgene_corrs/',module,".txt_pc_modgene_corrs.txt.gz"))
loadings$sign <- sign(loadings$PC)

eqtl.input <- list()
for (i in 1:length(genes)){
	cis.sub <- ciseqtl[ciseqtl$gene==genes[[i]] & ciseqtl$SNP %in% coeqtl.dat$SNP,]
	trans.sub <- transeqtl[transeqtl$gene==genes[[i]] & transeqtl$SNP %in% coeqtl.dat$SNP,]
	eqtl.input[[i]] <- rbind(cis.sub,trans.sub)
	eqtl.input[[i]] <- subset(eqtl.input[[i]],select=c('SNP','t-stat'))
	eqtl.input[[i]]$'t-stat.fin' <- (loadings[rownames(loadings) %in% genes[[i]],]$sign)*eqtl.input[[i]]$'t-stat'
	write.table(eqtl.input[[i]],paste0(genes[[i]],"_",coeqtl,".z"),quote=F,row.names=F,col.names=F)
}
 
### reformat input files out of R
for i in $(ls ENS*.z); do
	awk '{print $1" "$3}' ${i} > ${i}.fin
done
rm ENS*.z

for i in $(ls ENS*.fin | sed 's/.fin//g'); do
mv ${i}.fin ${i}
done









############################ 
#### create LD matrices in R
library(tidyverse)
library(data.table)

setwd('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/ecaviar/coloc_single_eqtl')

##### read in coeqtl files
filenames <- list.files(path="/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/coeqtl_maf0.05",pattern="*for_coloc.txt",full.names=T)
ldf <- lapply(filenames,read.table)
names(ldf) <- gsub("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/coeqtl_maf0.05/","",filenames)
ldf.fin <- list()
for (i in 1:length(ldf)) {
	colnames(ldf[[i]]) <- c('chr','pos','A1','A2','module','pc','beta','t.stat','p.val','FDR','MAF')
	ldf[[i]]$ID <- paste0(ldf[[i]]$chr,":",ldf[[i]]$pos)
}

# read in dosage
dat <- list()
for (i in 1:22){
	dat[[i]] <- read.table(paste0("../chr",i,"_coeqtl_esnp_dosage_subset"), header=TRUE)
}

dat <- Filter(function(x) dim(x)[1] > 0, dat)
dosage <- bind_rows(dat)
dosage$SNP <- paste0(dosage$CHR,":",dosage$POS)

gwas.df <- list()
corrs <- list()
for (i in 1:length(ldf)) {
	index <- ldf[[i]]$ID
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

for (i in 1:length(corrs)){
	write.table(corrs[[i]],paste0(names(ldf[i]),".ld"),quote=FALSE,row.names=FALSE,col.names=FALSE)
}




############################ 
## run eCaviar
# CLPP threshold of 0.001 is recommended by creator

ls ENS*.z | awk -F"_" '{print $0" "$2"_"$3"_"$4"_"$5"_"$6}' | sed 's/.z//g' > job.file.ecaviar.txt

for i in $(ls /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/coeqtl_maf0.05/*for_coloc.txt | xargs -n1 basename); do
awk '{print $1":"$2":"$3":"$4" "$8}' ${i} > /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/alanna/cibersort/cmc/ecaviar/coloc_single_eqtl/${i}.z
done


ml caviar
while read eqtl locus;do
bsub -J cav.${locus} -P acc_psychgen -q premium -R rusage[mem=40000] -W 96:00 -e e.${eqtl} -o o.${eqtl} eCAVIAR -l ${locus}_for_coloc.txt.ld -z ${eqtl}.z -l ${locus}_for_coloc.txt.ld -z ${locus}_for_coloc.txt.z -f 1 -c 3 -o ${eqtl}.output
done < job.file.ecaviar.txt

