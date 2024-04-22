source("samples_ann_preprocessing.R")
source("create_sajr_object.R")
source("create_DJexpress_object.R")

# ?????? что за бардак?

##########################################################################
# SAJR
##########################################################################

#  -----------------
rse = rse.ERP109002.jxn.cytosk.genes[,rse.ERP109002.jxn.cytosk.genes@colData$tissue %in% c('Brain','Testis') & 
                                       rse.ERP109002.jxn.cytosk.genes@colData$age_group=="adult"]

#----------------??
b = rownames(rse@colData[rse@colData$tissue=='Brain',])
t = rownames(rse@colData[rse@colData$tissue=='Testis',])
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

# что то тут сломалось!?
# там посчитал l2fc

#====================================================================================
# -----------------------------------DJExpress---------------------------------------
#====================================================================================
# Run DJEanalyze
dje_out = 
  makeDJexpress(condition.1 = 'Brain', condition.2='Testis', reference.condition='Brain', FDR.threshold=0, logFC.threashold=0)

#dje_out = anlz.out$dje.out

#str(anlz.out)
# Number of junctions in both methods
dim(dje_out)
dim(r)/2

################################################################################
# Comparison
################################################################################
# Merging SAJR and DJexpress results
# formatting junction id column, to use it for merging
strand_dict = c('+'=1, '-'=2, "*"=0)
r$junctionID = rownames(r)
r$junctionID = sapply(r$junctionID, function(x) sub('(\\d+)(-)(\\d+)', '\\1:\\3', x))
r$junctionID = sapply(r$junctionID, function(x) sub('([-*+])(:)(.)+', '\\1', x))
r$junctionID = sapply(r$junctionID, function(x) sub('.$', strand_dict[sub('.*(?=.$)', '', x, perl=T)], x))
m = merge(r, dje_out, by='junctionID')

# Plots
plot(m$dpsi, m$logFC, xlab = "dpsi", ylab = "logFC")
plot(m$pv, m$P.Value, log='xy', xlab = "pv", ylab = "P.Value")
plot(m$l2fc, m$logFC)

# Correlation
p=na.omit(m[c('pv', 'P.Value')])
a=na.omit(m[c('dpsi', 'logFC')])
d = na.omit(m[c('l2fc', 'logFC')])
cor(p$pv, p$P.Value, method = "spearman")
cor(a$dpsi, a$logFC, method = "spearman")
cor(d$l2fc, d$logFC, method = "spearman")


# Number of events not represented on the plot
table(is.na(m$dpsi))
table(is.na(m$logFC))


# ================= DIEGO =====================
dje_out$junctionIDDiego = 
  sapply(dje_out$junctionID, function(x) sub('(\\d{3,})(:)(\\d{3,})', '\\1-\\3', x))
dje_out$junctionIDDiego = 
  sapply(dje_out$junctionIDDiego, function(x) sub('(\\d{3,})(:)(\\d)', '\\1', x))

r$junctionIDDiego = 
  sapply(rownames(r), function(x) sub('(\\d{3,})(:.+)', '\\1', x))

diego_output = read.table('./DIEGO_output', sep = "\t", header = TRUE)
dim(diego_output)

colnames(diego_output)[colnames(diego_output) == "junction"] <- "junctionIDDiego"

m = merge(diego_output, dje_out, by='junctionIDDiego')
cor(m$abundance_change, m$logFC, method = "spearman")

plot(m$abundance_change, m$logFC)

m = merge(diego_output, r, by='junctionIDDiego')
a=na.omit(m[c('dpsi', 'abundance_change')])

cor(a$abundance_change, a$dpsi, method = "spearman")
plot(m$abundance_change, m$dpsi)


head(dje_out)

#junction
