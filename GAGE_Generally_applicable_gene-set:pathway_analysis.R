#to clean global environment (optional)
rm(list=ls())

#installing required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

install.packages("openxlsx")

library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)

#set work directory
directory <- setwd('/path to the dir where the count data is stored/')


#list all files from HTSEq counts in directory
samples <- list.files(recursive = TRUE, pattern = "*_TSMRS.htSeqCounts") #Select the pattern in which the file is stored
samples

htseq.results <- data.frame() #to rearrange data in the required dataframe

# Here’s the magic loop. For each sample, grab the read the input file, extract the TPM column, and append it to the tpm.results data.frame.
for ( i in 1: length(samples)) {
  filename <- samples[i]
  sample <- sub("_FB_Whole_F.*G.*S.*M.*P.*_TSMRS.htSeqCounts", "", filename) #input file pattern
  message("Now running sample: ", sample)
  
  tpm.file <- read.table(filename, header = TRUE)
  names(tpm.file) <- c("EnsemblID", sample)
  tpm.file <- tpm.file[, c("EnsemblID",sample)]
  
  if (nrow(htseq.results) == 0) {
    htseq.results <- tpm.file
  } else {
    htseq.results <- merge(htseq.results, tpm.file, by = "EnsemblID")
  }
}
#converting fies to numeric type

htseq.results$filename without extenstion <-as.numeric(htseq.results$filename without extenstion)
#repeat for all the sample files

library(GenomicAlignments)
require(gageData)

#to get read counts
cnts=htseq.results[first data col num : last data col num] #specifying column numbers. col 1 is serial num so needs to be excluded
cnts

sel.rn=rowSums(cnts) != 0 #selects no-zero rows
sel.rn
cnts=cnts[sel.rn,] #gives counts of non-zero rows
dim(cnts)

#normalizing read counts 
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
range(cnts.norm)
cnts.norm=log2(cnts.norm+8)#cnts.norm is the normalized read count. Addition of 8 throughout to do away with very small values
range(cnts.norm) 

#MA plot to check for processed data variances (optional)

pdf("filename.pdf", width=8, height=10)
op=par(lwd=2, cex.axis=1.5, cex.lab=1.5, mfrow=c(2,1))
plot((cnts.norm[,6]+cnts.norm[,5])/2, (cnts.norm[,6]-cnts.norm[,5]),main="(a) Control vs Control", xlab="mean", ylab="change",ylim=c(-5,5), xlim=c(0,20), lwd=1)
abline(h=0, lwd=2, col="red", lty="dashed")
plot((cnts.norm[,1]+cnts.norm[,5])/2, (cnts.norm[,1]-cnts.norm[,5]), main="(b) Knockdown vs Control", xlab="mean", ylab="change",ylim=c(-5,5), xlim=c(0,20), lwd=1)
abline(h=0, lwd=2, col="red", lty="dashed")
dev.off()

#Gene set test
require(gageData)
library(gage)
ref.idx=first controls col: last controls col #col num of controls
samp.idx=first cases col : last cases col #col num of cases
data(kegg.gs) #pulled pathway data from KEGG

##generates the required gene enrichment data
#Say, cases and control samples are unpaired depending on research
cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = ref.idx,samp = samp.idx, compare ="unpaired")
cnts.kegg.p #lists dysregulated pathways

#export the file and/or object in desired format(optional)
library(openxlsx)
write.xlsx(cnts.kegg.p, "path where it will be stored.xlsx",row.names = TRUE, col.names = TRUE) #To download file in excel format
saveRDS(cnts.kegg.p, file = "path.rds") #to save as RDS object

##To generate heatmap and scatterplot with the output (optional)
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
require(gageData)
library(gage)
for (gs in rownames(cnts.kegg.p$greater)[1:n]) #to obtain plots for the pathways upregulated upto top n
{
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = cnts.norm, ref = ref.idx,
           samp = samp.idx, outname = outname, txt = TRUE, heatmap = TRUE,
           Colv = TRUE, Rowv = TRUE, dendrogram = NULL, limit = 3, scatterplot = TRUE)
}


#differential expression: log2 ratio or fold change, unpaired samples (optional. not required if already ran diff exp analysis separately)
cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])
cnts.d

#up-regulated pathways (top n) visualized by pathview- annotates in KEGG
sel <- cnts.kegg.p$greater[, "q.val"] < selecting the rows with the desired q-values under 'greater' section & !is.na(cnts.kegg.p$greater[,"q.val"]) #excludes not applicable ones out
sel
path.ids <- rownames(cnts.kegg.p$greater)[sel]#assign rownames
path.ids
path.ids2 <- substr(path.ids, 1, n) #n is the number of pathways you want to visualize
path.ids2

library(pathview)
pv.out.list <- sapply(path.ids2[1:n], function(pid) pathview(gene.data =cnts.d, pathway.id = pid,species = "hsa"))#maps the genes in the pathway. creates file in wd
#kegg.sig<-sigGeneSet(cnts.kegg.p, outname="/Volumes/narayanan/abhowmik/F0530_htSeqCounts/heatmap.kegg",pdf.size = c(7,12))#maybe used based on the results. gives the dyregulated pathways on the basis of default q<0.1

#repeat same with down-regulated pathway by looking for q-val> threshold in the $lesser section 


##____________________________________________________________end of KEGG__________________________________________________________________________________________________________##

###Gene Ontology Analysis
library(gageData)
#calls data from GO database
data(go.sets.hs)
data(go.subs.hs)
lapply(go.subs.hs, head)


##For molcular functions
cnts.mf.p <- gage(cnts.norm, gsets = go.sets.hs[go.subs.hs$MF],ref = ref.idx, samp = samp.idx, compare ="unpaired")
#gives dysregulated molecular functions. export object and/ or write it in output dir (Optional)

#Biological Process analysis takes a few minutes if you try it
cnts.bp.p <- gage(cnts.norm, gsets = go.sets.hs[go.subs.hs$BP],
                  ref = ref.idx, samp = samp.idx, compare ="unpaired") #takes some time
#gives dysregulated molecular functions. export object and/ or write it in output dir (Optional)
#job takes time depending on the sample size. Running parallely in multiple cores may be helpful 
#export object and/ or write it in output dir (Optional)



# generate heatmap and scatter plot (optional)
for (gs in rownames(cnts.bp.p$greater/lesser)[1:n]) 
  outname = gsub(" |:|/", "_", substr(gs, 12, 100))
geneData(genes = go.sets.hs[[gs]], exprs = cnts.norm, ref = ref.idx,samp = samp.idx, outname = outname, txt = T, heatmap = T,limit = 3, scatterplot = T)


##______________________________________________________________________End of GO analysis__________________________________________________________________________##

##By default GO performs pair-wise comparisons between samples. however other per gene score from group wise comparison can also be used
#Per gene score choices- 

#Eg: t-test and mean fold change

cnts.t= apply(cnts.norm, 1, function(x) t.test(x[samp.idx], x[ref.idx],alternative = "two.sided", paired = F)$statistic)#gives the t-test stats
range(cnts.t)
cnts.t.kegg.p <- gage(cnts.t, gsets = kegg.gs, ref = NULL, samp = NULL) 
cnts.t.kegg.p

cnts.meanfc= rowMeans(cnts.norm[, samp.idx])-rowMeans(cnts.norm[, ref.idx])#gives the mean fold change stats
range(cnts.meanfc)
cnts.meanfc.kegg.p <- gage(cnts.meanfc, gsets = kegg.gs, ref = NULL, samp = NULL)
cnts.meanfc.kegg.p

