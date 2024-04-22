#source("data.R")
source("samples_ann_preprocessing.R")
library(DJExpress)

# 
# strand_dict = c('+'=1, '-'=2, "*"=0)
# # replace +/- strand with 1/2 (0: undefined, 1: +, 2: -)
# # STAR manual, p.12 on output splice junctions file
# chain_to_number = function (coordinates) {
#   strand_sign = sub('.*(?=.$)', '', coordinates, perl=T) # replace everything aside from the last character with ''. ?= look ahead. It performs the match, but does not capture the match
#   coordinates = sub('([0-9]+)(-)([0-9]+)', '\\1:\\3' , coordinates)
#   coordinates = sub('.$', strand_dict[strand_sign], coordinates)
# }

# # ??? пока только для тканей, потом для возраста исправить код
# makeDJexpress = function(rse, reference.condition, 
#                               minMean, minVar,
#                               FDR.threshold, logFC.threashold){
#   
#   counts = as.matrix(rse@assays@data$counts)
#   rownames(counts) = sapply(row.names(counts), chain_to_number) # formatting row ids
#   
#   reference.samples.ids =
#     rownames(rse@colData[rse@colData$age_group==reference.condition,])
#   
#   ann.out = list(quant.annotated = counts,
#                  featureID = rownames(counts),
#                  groupID = rse@rowRanges$gene_names)
#   
#   prep.out <- DJEprepare(annotate.out = ann.out, 
#                          Group1 = reference.samples.ids,
#                          minMean = minMean,
#                          maxMean = Inf,
#                          minVar = minVar,
#                          maxVar = Inf) # these are the default values for minMean, maxMean, minVar and MaxVar
#   
#   print(summary(prep.out))
#   
#   # Run DJEanalyze
#   anlz.out <- DJEanalyze(prepare.out = prep.out,
#                          Group1 = reference.samples.ids,
#                          FDR = FDR.threshold,
#                          logFC = logFC.threashold)
#   anlz.out
# }
# 

# same result if preparing data manually

strand_dict = c('+'=1, '-'=2, "*"=0)
# replace +/- strand with 1/2 (0: undefined, 1: +, 2: -)
# STAR manual, p.12 on output splice junctions file
chain_to_number = function (coordinates) {
  strand_sign = sub('.*(?=.$)', '', coordinates, perl=T) # replace everything aside from the last character with ''. ?= look ahead. It performs the match, but does not capture the match
  coordinates = sub('([0-9]+)(-)([0-9]+)', '\\1:\\3' , coordinates)
  coordinates = sub('.$', strand_dict[strand_sign], coordinates)
}

# ??? пока только для тканей, потом для возраста исправить код
makeDJexpress = function(rse, reference.condition,
                         minMean, minVar,
                         FDR.threshold, logFC.threashold){

  # rse = rse[apply(rse@assays@data$counts,1,mean) >= 10 &
  #           apply(rse@assays@data$counts, 1, var) >= 0 ]

  counts = as.matrix(rse@assays@data$counts)
  rownames(counts) = sapply(row.names(counts), chain_to_number) # formatting row ids

  out.file = list(JunctExprfilt = counts,
                  featureID = rownames(counts),
                  groupID = rse@rowRanges$gene_names,
                  design = model.matrix(~rse@colData$age_group))

  print(summary(out.file))

  reference.samples.ids =
    rownames(rse@colData[rse@colData$age_group==reference.condition,])


  # Run DJEanalyze
  anlz.out <- DJEanalyze(prepare.out = out.file,
                         Group1 = reference.samples.ids,
                         FDR = FDR.threshold,
                         logFC = logFC.threashold)
  anlz.out
}

# 
# bt.self.prep = makeDJexpress(rse, reference.condition='Brain', 
#                              minMean=NA, minVar=NA,
#                              FDR.threshold=0.05, logFC.threashold=2)
# bt.self.prep$dje.out$junctionID %in% bt_prep$dje.out$junctionID
# 
# ids.not.in.prep=bt.self.prep$dje.out$junctionID[!bt.self.prep$dje.out$junctionID %in% bt_prep$dje.out$junctionID]
# ids.in.both=bt.self.prep$dje.out$junctionID[bt.self.prep$dje.out$junctionID %in% bt_prep$dje.out$junctionID]

# 
# counts[ids.not.in.prep,]
# counts[ids.in.both,]

# сравнить с каунтами в sajr
