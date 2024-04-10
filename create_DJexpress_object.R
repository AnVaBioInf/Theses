library(DJExpress)

source("data.R")
source("samples_ann_preprocessing.R")

#rownames(samples_meta) = paste0(rownames(samples_meta), "_", samples_meta$tissue)

#  metadata for brain and testis samples
samples_bt = samples_meta[samples_meta$tissue %in% c('Brain','Testis') & samples_meta$age_group=="adult",]


# these are actual read counts and not raw base-pair counts !!! what does this mean?
counts = rse.ERP109002.jxn.cytosk.genes@assays@data$counts[,rownames(samples_bt)] %>%
  as.matrix() 

#=========
# formatting row ids
# replace +/- strand with 1/2 (0: undefined, 1: +, 2: -)
# STAR manual, p.12 on output splice junctions file

strand_dict = c('+'=1, '-'=2, "*"=0)
chain_to_number = function (coordinates) {
  strand_sign = sub('.*(?=.$)', '', coordinates, perl=T) # replace everything aside from the last character with ''. ?= look ahead. It performs the match, but does not capture the match
  coordinates = sub('([0-9]+)(-)([0-9]+)', '\\1:\\3' , coordinates)
  coordinates = sub('.$', strand_dict[strand_sign], coordinates)
}

rownames(counts) = sapply(row.names(counts),
                         chain_to_number)

out.file = list(JunctExprfilt = counts, 
                featureID = rownames(counts),
                groupID = rse.ERP109002.jxn.cytosk.genes@rowRanges$gene_names,
                design = model.matrix(~samples_bt$tissue))

summary(out.file)
out.file$featureID # why only x chr?

#
b = rownames(samples_bt[samples_bt$tissue=='Brain',])

# Run DJEanalyze
anlz.out <- DJEanalyze(prepare.out = out.file, Group1 = b,
                       FDR = 0, logFC = 0)

anlz.out$dje.out

