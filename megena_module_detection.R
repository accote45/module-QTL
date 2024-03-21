

library(MEGENA)
library(tidyverse)
library(dplyr)

# input parameters
n.cores <- 15; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = "pearson" # method for correlation. either pearson or spearman.
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples.
module.pval = 0.05 # module significance p-value. Recommended is 0.05.
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs.
hub.perm = 100; # number of permutations for calculating connectivity significance p-value.

# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2
###########

datExpr <- read.table("/sc/arion/projects/psychgen/alanna/cibersort/knowncovar_+cellprop_adj_quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt")
#datExpr <- as.data.frame(t(datExpr))

ijw <- calculate.correlation(datExpr,doPerm = cor.perm,output.corTable = FALSE,output.permFDR = FALSE)


#### register multiple cores if needed: note that set.parallel.backend() is deprecated.
run.par = doPar & (getDoParWorkers() == 1)
if (run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}

##### calculate PFN
el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores,keep.track = FALSE)

g <- graph.data.frame(el,directed = FALSE)

##### perform MCA clustering.
MEGENA.output <- do.MEGENA(g,
 mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
 min.size = 10,max.size = vcount(g)/2,
 doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
 save.output = FALSE)

###### unregister cores as these are not needed anymore.
if (getDoParWorkers() > 1)
{
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


summary.output <- MEGENA.ModuleSummary(MEGENA.output,
    mod.pvalue = module.pval,hub.pvalue = hub.pval,
    min.size = 10,max.size = vcount(g)/2,
    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
    output.sig = TRUE)

if (!is.null(annot.table))
{
  # update annotation to map to gene symbols
  V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
  summary.output <- output[c("mapped.modules","module.table")]
  names(summary.output)[1] <- "modules"
}


for (i in seq_along(summary.output[[1]])) {
    filename=paste0(names(summary.output[[1]][i]),".txt")
    write.table(summary.output[[1]][i], filename, quote=FALSE, row.names=FALSE)
}

write.table(summary.output[[2]], "module.table.txt", quote=FALSE, row.names=TRUE)



