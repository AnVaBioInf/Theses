source("samples_ann_preprocessing.R")
source("create_sajr_object.R")
source("create_DJexpress_object.R")


# install.packages("dynamicTreeCut")
# BiocManager::install("edgeR")
# install.packages("fastcluster")
# install.packages("highcharter")
# BiocManager::install("impute")
# BiocManager::install("preprocessCore")
# install.packages("WGCNA")
# install.packages('MatrixEQTL')
# install.packages('/home/an/bin/DJExpress_0.1.0.tar.gz', repos = NULL, type="source")
library(DJExpress)
library("recount3")
library("dbplyr")
library("tidyr")
library("readr")
library('tibble')

library("dplyr")
library("gtools")
library('gridExtra')
library("tidyr")
library("ggplot2")
#/??

#  metadata for brain and testies samples
samples_bt = samples_meta[samples_meta$tissue %in% c('Brain','Testis') & samples_meta$age_group=="adult",]

b = rownames(samples_bt[samples_bt$tissue=='Brain',])
t = rownames(samples_bt[samples_bt$tissue=='Testis',])


##########################################################################
# SAJR
##########################################################################

#  -----------------
brain_adult_samples = samples_meta[samples_meta$tissue=="Brain" & samples_meta$age_group=="adult",]
testies_adult_samples = samples_meta[samples_meta$tissue=="Testis" & samples_meta$age_group=="adult",]

#----------------??
b = rownames(brain_adult_samples)
t = rownames(testies_adult_samples)
b
t

# the Binomial Generalized Linear Model was applied to test the segnificance of the inclusion change with age
r = t(apply(sajr$ir,1,function(x){ # каждая строчка
  pv =NA
  if(sum(!is.na(x[b]))>1 & sum(!is.na(x[t]))>1)
    pv = wilcox.test(x[b],x[t])$p.value
  c(pv=pv, dpsi=mean(x[t],na.rm=T)-mean(x[b],na.rm=T), l2fc=log2(mean(x[t], na.rm=T)/mean(x[b],na.rm=T)))
}))

r = as.data.frame(r)
r$fdr = p.adjust(r$pv,m='BH')
table(is.na(r[,1]))

plot(r$l2fc,1/r$pv,log='y',pch=16,col=(r$fdr<0.05 & abs(r$dpsi)>0.2) + 1)
plot(r$dpsi,1/r$pv,log='y',pch=16,col=(r$fdr<0.05 & abs(r$dpsi)>0.2) + 1)

sgn = (r$fdr<0.05 & abs(r$dpsi)>0.4)
sgn[is.na(sgn)]  = F
sum(sgn,na.rm=T)
# [1] 22
unique(sajr$seg$gene_names[sgn])
# [1] "WASF3"  "ACTR2"  "ABI2"   "WASF1"  "ACTR3B"




#====================================================================================
# -----------------------------------DJExpress---------------------------------------
#====================================================================================
# Run DJEanalyze
anlz.out <- DJEanalyze(prepare.out = out.file, Group1 = b,
                       FDR = 0, logFC = 0)

dje_out = anlz.out$dje.out

dim(anlz.out$dje.out)
dim(r)/2

strand_dict = c('+'=1, '-'=2, "*"=0)
r$junctionID = rownames(r)
r$junctionID = sapply(r$junctionID, function(x) sub('(\\d+)(-)(\\d+)', '\\1:\\3', x))
r$junctionID = sapply(r$junctionID, function(x) sub('([-*+])(:)(.)+', '\\1', x))
r$junctionID = sapply(r$junctionID, function(x) sub('.$', strand_dict[sub('.*(?=.$)', '', x, perl=T)], x))

m = merge(r, dje_out, by='junctionID')

plot(m$dpsi, m$logFC, log='xy', xlab = "dpsi", ylab = "logFC")
plot(m$pv, m$P.Value, log='xy', xlab = "pv", ylab = "P.Value")

p=na.omit(m[c('pv', 'P.Value')])
a=na.omit(m[c('dpsi', 'logFC')])

cor(d$pv, d$P.Value, method = "spearman")
cor(a$dpsi, a$logFC, method = "spearman")

