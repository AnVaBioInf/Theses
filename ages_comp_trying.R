#devtools::install_github('iaaka/sajr')
library('SAJR')
source("create_sajr_object.R")
source("create_DJexpress_object.R")
#source("samples_ann_preprocessing.R")
source('plots.R')

# взять одинаковые списки джанкшенов
# ir всегда положительная величина

unique.tissues = unique(rse.ERP109002.jxn.cytosk.genes@colData$tissue) # we excluded ovaries from analysis -> not enought of samples

findPairsInfo = function(tissue){
  samples.ann = rse.ERP109002.jxn.cytosk.genes@colData
  adult.samples.ids = rownames(samples.ann[samples.ann$tissue %in% tissue & samples.ann$age_group == 'adult',])
  fetus.samples.ids = rownames(samples.ann[samples.ann$tissue %in% tissue & samples.ann$age_group == 'fetus',])
  list(adult.samples.ids=adult.samples.ids, fetus.samples.ids=fetus.samples.ids)
}

# NULL???

# find all junction ids that fulfill the criteria
filterJxns = function(tissue){
  jxn.ids.passed = c()
  
  for (tissue in unique.tissues){
    age.groups.ids.tissues = findPairsInfo(tissue)
    fetus.ids = age.groups.ids.tissues$fetus.samples.ids
    adult.ids = age.groups.ids.tissues$adult.samples.ids
    
    sajr.age.pairs = sajr[ ,c(fetus.ids,adult.ids)] # selecting only samples for the tissue
    # finding jxn ids for every age pair, fulfilling the criteria
    sajr.age.pairs =
      sajr.age.pairs[
        apply(sajr.age.pairs$i+sajr.age.pairs$e>=10,1,sum, na.rm = TRUE)>=2 &
        apply(sajr.age.pairs$ir>0,1,sum, na.rm = TRUE)>=2 &
        apply(sajr.age.pairs$ir,1,sd,na.rm=TRUE) > 0,]

    # making a list of all ids for all tissues that passed the criteria
    jxn.ids.passed = unique(c(jxn.ids.passed, rownames(sajr.age.pairs$ir)))
  }
  jxn.ids.passed
}

# sajr jxns, for rse must be converted!
covertSajrtoRse = function(jxns.sajr.ids){
  sapply(jxns.sajr.ids, function(x){
    sub("(:[-+*])(.+)", "\\1", x)
  })
}

prepareData = function(tissue){
  age.groups.ids.tissues = findPairsInfo(tissue)
  fetus.ids = age.groups.ids.tissues$fetus.samples.ids
  adult.ids = age.groups.ids.tissues$adult.samples.ids

  filtered.jxns = filterJxns()
  sajr.filtered = sajr[filtered.jxns, c(fetus.ids,adult.ids)]
  
  rse.filtered.jxn.ids = covertSajrtoRse(filtered.jxns) # unque, because sajr has r/l for each junction
  rse.filtered.jxn.ids = unique(rse.filtered.jxn.ids)
  rse.filtered = rse.ERP109002.jxn.cytosk.genes[rse.filtered.jxn.ids, ]
  
  mod = c(rep("fetus", length(fetus.ids)), rep("adult", length(fetus.ids)))
  mod = list(age=factor(mod)) # ~ model data
  
  list(sajr.filtered=sajr.filtered, mod=mod, rse.filtered=rse.filtered)
}

calculateMetrics = function(sajr.prepared, fetus.ids, adult.ids){
  psi = apply(sajr.age.pairs$ir, 1, function(x) mean(x[adult.ids],na.rm=T)-mean(x[fetus.ids],na.rm=T))
  logFC = apply(sajr.age.pairs$ir, 1, function(x) log2(mean(x[adult.ids],na.rm=T)/mean(x[fetus.ids],na.rm=T)))
  list(psi, logFC)
}

runSAJR = function(){
  sajr.output = list()
  for (tissue in c('Brain')){  # for every tissue
    pair.ids = findPairsInfo(tissue)
    prepared.Data = prepareData(tissue)
    mod = prepared.Data$mod
    sajr.filtered = prepared.Data$sajr.filtered
    
    
    alt.glm = as.data.frame( fitSAGLM(sajr.filtered, terms(x ~ age, keep.order=T),mod,return.pv=T) )
    
    # metrics = calculateMetrics(sajr.filtered, pair.ids$fetus.ids, pair.ids$adult.ids)
    # alt.glm$dPSI = metrics$psi
    # alt.glm$logFC.sajr = metrics$logFC
    alt.glm$FDR.sajr = p.adjust(alt.glm$age,m='BH')
    
    sajr.output[[tissue]] = alt.glm
    print(head(sajr.filtered))
  }
  print(c("output",sajr.output))
  #sajr.output
}

runSAJR()




runDJE = function(){
  dje.output = list()
  for (tissue in unique.tissues){
    
    samples.pairs.annotation = rse.ERP109002.jxn.cytosk.genes[,rse.ERP109002.jxn.cytosk.genes@colData$tissue %in% tissue]
    anlz.out = makeDJexpress(samples.pairs.annotation, reference.condition='fetus', FDR.threshold=0.05, logFC.threashold=2)
    dje.output[[tissue]] = anlz.out$dje.out[,c('logFC', 'P.Value', 'FDR', 'junctionID')]
  }
}

findSignificantEvents = function(tool.output, par, fdr, par.threshold, fdr.threshold=0.05){
  lapply(tool.output, function(x) {
    x[abs(x[,par]) >= par.threshold & x[,fdr] <= fdr.threshold,]
  })
}

sajr_to_djexpress_jxn_coord = function(sjar.pair.output.i) {
  strand_dict = c('+' = 1, '-' = 2, '*' = 0)
  sjar.pair.output.i$junctionID = rownames(sjar.pair.output.i)
  sjar.pair.output.i$junctionID = sub("(\\d{3,})(-)(\\d{3,})(:[-+*])(:).+", "\\1:\\3\\4", sjar.pair.output.i$junctionID)
  sjar.pair.output.i$junctionID = 
    sapply(sjar.pair.output.i$junctionID, function(x) sub('.$', strand_dict[sub('.*(?=.$)', '', x, perl=T)], x))
  sjar.pair.output.i
}

processOutputs = function(){
  # убрать этот колл в run sajr?
  # change junction ids in DJexpress to match sajr
  sajr.output = lapply(sajr.output, function(x) sajr_to_djexpress_jxn_coord(x))
  
  # filtering output
  sajr.significant.events = findSignificantEvents(sajr.output, 'dPSI', 'FDR.sajr', 0.1) 
  dje.significant.events = findSignificantEvents(dje.output, 'logFC', 'FDR', 2) 
}

# убрать ненужное!
groupOutputs = function(){
  outputs.compared <- lapply(unique.tissues, function(tissue) {
    sajr.sign.tissue = sajr.significant.events[[tissue]]
    dje.sign.tissue = dje.significant.events[[tissue]]
    sajr.output.tissue = sajr.output[[tissue]]
    dje.output.tissue = dje.output[[tissue]]
    
    all.common.jxns = merge(sajr.output.tissue, dje.output.tissue, by = "junctionID", all = FALSE)
    sign.both.tools = merge(sajr.sign.tissue, dje.sign.tissue, by = "junctionID", all = FALSE)
    dje.not.common.jxns = dje.output.tissue[!dje.output.tissue$junctionID %in% all.common.jxns$junctionID,]
    sajr.not.common.jxns = sajr.output.tissue[!sajr.output.tissue$junctionID %in% all.common.jxns$junctionID,]
    
    # тут косфк
    list(
      all.common.jxns = all.common.jxns,
      dje.not.common.jxns = dje.not.common.jxns,
      sajr.not.common.jxns = sajr.not.common.jxns,
      
      sign.both.tools = sign.both.tools,
      sign.sajr.only = sajr.sign.tissue[!sajr.sign.tissue$junctionID %in% dje.sign.tissue$junctionID, ],
      sign.dje.only = dje.sign.tissue[!dje.sign.tissue$junctionID %in% sajr.sign.tissue$junctionID, ],
      
      dje.sign.common = dje.sign.tissue[dje.sign.tissue$junctionID %in% all.common.jxns$junctionID,],
      sajr.sign.common = sajr.sign.tissue[sajr.sign.tissue$junctionID %in% all.common.jxns$junctionID,],
      dje.sign.not.common = dje.sign.tissue[dje.sign.tissue$junctionID %in% dje.not.common.jxns$junctionID,],
      sajr.sign.not.common = sajr.sign.tissue[sajr.sign.tissue$junctionID %in% sajr.not.common.jxns$junctionID,]
    )
  })
  # Set names of the list elements to tissue names
  names(outputs.compared) <- names(sajr.significant.events)
}


# graphs
# significant only in djexpress
# significant only in sajr
# not significant in djexpress
# not dignificant in sajr
# significant in both-dj express and sajr

# all.common.jxns = outputs.compared.tissue$all.common.jxns
# dje.only.sign = outputs.compared.tissue$dje.only.sign
# sajr.only.sign = outputs.compared.tissue$sajr.only.sign
# 
# all.common.jxns[all.common.jxns$junctionID %in% sajr.only.sign$junctionID,][,par.1]


