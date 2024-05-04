#devtools::install_github('iaaka/sajr')
library('SAJR')
source("create_sajr_object.R")
source("create_DJexpress_object.R")
source("create_Diego_object.R")

#source("samples_ann_preprocessing.R")
#source('plots.R')

# взять одинаковые списки джанкшенов
# ir всегда положительная величина

unique.tissues = unique(rse.ERP109002.jxn.cytosk.genes@colData$tissue)
tissue.pairs = combn(unique.tissues, 2, simplify = FALSE)


# running sajr for every pair of tissues

# sajr jxns, for rse must be converted!
covertSajrtoRse = function(jxns.sajr.ids){
  sapply(jxns.sajr.ids, function(x){
    sub("(:[-+*])(.+)", "\\1", x)
  })
}

unifyCoordinates = function(programm.output, coordinates) {
  programm.output$unifiedJxnID = coordinates
  programm.output$unifiedJxnID = sub("(\\d{3,})([-:])(\\d{3,}).+", "\\1-\\3", programm.output$unifiedJxnID)
  programm.output
}

findPairsInfo = function(){
  tissue.groups.ids.tissues = list()
  samples.ann = rse.ERP109002.jxn.cytosk.genes@colData
  for (tissue.pair in tissue.pairs){
    tissue1.samples.ids = rownames(samples.ann[samples.ann$tissue == tissue.pair[1] & samples.ann$tissue == 'adult',])
    tissue2.samples.ids = rownames(samples.ann[samples.ann$tissue == tissue.pair[2] & samples.ann$tissue == 'adult',])
    tissue.groups.ids.tissues[[tissue.pair]]$tissue1.samples.ids = tissue1.samples.ids
    tissue.groups.ids.tissues[[tissue.pair]]$tissue2.samples.ids = tissue2.samples.ids
  }
  tissue.groups.ids.tissues
}

# find all junction ids that fulfill the criteria
filterJxns = function(){
  jxns.sajr.ids.passed = c()
  
  for (tissue.pair in tissue.pairs){
    tissue1.ids = tissue.groups.ids.tissues[[tissue.pair]]$tissue1.samples.ids
    tissue2.ids = tissue.groups.ids.tissues[[tissue.pair]]$tissue2.samples.ids
    
    sajr.tissue.pairs = sajr[ ,c(tissue1.ids,tissue2.ids)] # selecting only samples for the tissue
    # finding jxn ids for every age pair, fulfilling the criteria
    sajr.tissue.pairs =
      sajr.tissue.pairs[
        apply(sajr.tissue.pairs$i+sajr.tissue.pairs$e>=10,1,sum, na.rm = TRUE)>=2 &
          apply(sajr.tissue.pairs$ir>0,1,sum, na.rm = TRUE)>=2 &
          apply(sajr.tissue.pairs$ir,1,sd,na.rm=TRUE) > 0,]
    
    # making a list of all ids for all tissues that passed the criteria
    jxns.sajr.ids.passed = unique(c(jxns.sajr.ids.passed, rownames(sajr.tissue.pairs$ir)))
  }
  jxns.ids.passed = covertSajrtoRse(jxns.sajr.ids.passed) # unque, because sajr has r/l for each junction
  names(jxns.ids.passed) = jxns.sajr.ids.passed
  
  jxns.ids.passed
}

prepareData = function(tissue.pair){
  tissue1.ids = tissue.groups.ids.tissues[[tissue.pair]]$tissue1.samples.ids
  tissue2.ids = tissue.groups.ids.tissues[[tissue.pair]]$tissue2.samples.ids
  sajr.jxns.ids.passed = names(jxns.ids.passed)
  rse.jxns.ids.passed = unique(jxns.ids.passed)
  
  sajr.filtered = sajr[sajr.jxns.ids.passed, c(tissue1.ids,tissue2.ids)]
  rse.filtered = rse.ERP109002.jxn.cytosk.genes[rse.jxns.ids.passed, c(tissue1.ids,tissue2.ids)]
  
  mod = c(rep(tissue.pair[1], length(tissue1.ids)),rep(tissue.pair[2], length(tissue2.ids)))
  mod = list(tissue=factor(mod)) # ~ model data
  
  list(sajr.filtered=sajr.filtered, mod=mod, rse.filtered=rse.filtered)
}

calculateMetrics = function(sajr.filtered, tissue1.samples.ids, tissue2.samples.ids){
  dPSI = apply(sajr.filtered$ir, 1, function(x) mean(x[tissue1.samples.ids],na.rm=T)-mean(x[tissue2.samples.ids],na.rm=T))
  #logFC = apply(sajr.filtered$ir, 1, function(x) log2(mean(x[adult.samples.ids],na.rm=T)/mean(x[fetus.samples.ids],na.rm=T)))
  list(dPSI=dPSI) #, logFC=logFC)
}

runSAJR = function(){
  sajr.output = list()
  for (tissue.pair in tissue.pairs){  # for every tissue
    prepared.data = prepareData(tissue.pair)
    mod = prepared.data$mod
    sajr.filtered = prepared.data$sajr.filtered
    
    alt.glm = as.data.frame( fitSAGLM(sajr.filtered, terms(x ~ age, keep.order=T),mod,return.pv=T) )
    metrics = calculateMetrics(sajr.filtered, 
                               tissue.groups.ids.tissues[[tissue.pair]]$tissue1.samples.ids
                               tissue.groups.ids.tissues[[tissue.pair]]$tissue2.samples.ids)
    alt.glm$dPSI = metrics$dPSI
    alt.glm$FDR.sajr = p.adjust(alt.glm$age,m='BH')
    alt.glm = unifyCoordinates(alt.glm, rownames(alt.glm))
    names(alt.glm)[names(alt.glm) == "age"] = "p.value.sajr"
    #alt.glm = alt.glm[!is.na(alt.glm$dPSI) & !is.na(alt.glm$p.value.sajr),]
    sajr.output[[tissue.pair]] = alt.glm
    
    print(c(tissue.pair,dim(alt.glm)))
    # # Identify rows with NA values
    # rows_with_na <- rowSums(is.na(sajr.output[[tissue]])) > 0
    # # Print rows with NA
    # print(sajr.output[[tissue]][rows_with_na, ])
    # 
  }
  sajr.output
}
# logFC 0/0 produces NaN
# NA in all samples of one pair (fetus/adult) produces dPSI NaN and logFC NaN
# Inf in logFC??

runDJE = function(){
  dje.output = list()
  for (tissue.pair in tissue.pairs){
    prepared.data = prepareData(tissue.pair)
    rse.filtered = prepared.data$rse.filtered
    #print(dim(rse.filtered))
    anlz.out = makeDJexpress(rse.filtered, reference.condition=tissue.pair[[2]], 
                             minMean=NA, minVar=NA,
                             FDR.threshold=0.1, logFC.threashold=1.5)
    anlz.out = anlz.out$dje.out[,c('junctionID', 'GeneID', 'logFC', 'P.Value', 'FDR')]
    anlz.out = unifyCoordinates(anlz.out, anlz.out$junctionID)
    
    #anlz.out = anlz.out[!is.na(anlz.out$logFC) & !is.na(anlz.out$P.Value),]
    dje.output[[tissue.pair]] = anlz.out
    print(c(tissue.pair,dim(anlz.out)))
    # почему то 1 джанкшен отваливается, но в sajr он не значимый
    
  }
  dje.output
}

makeDiegoFiles = function(){
  for (tissue.pair in tissue.pairs){
    prepared.data = prepareData(tissue.pair)
    rse.filtered = prepared.data$rse.filtered
    makeDiegoInputFiles(rse.filtered,tissue.pair)
    
    # diego_output = read.table('./DIEGO_output', sep = "\t", header = TRUE)
    # dim(diego_output)
    # 
    # head(diego_output)
    # 
  }
}

readDiegoOutput = function(){
  
  # Create an empty list to store data frames
  diego.output = list()
  # Loop through tissues and read files
  for (tissue.pair in tissue.pairs) {
    # Construct file path
    file_path = paste0("/home/an/DIEGO_output_files/DIEGO_output", tissue.pair)  # Replace with actual path
    
    # Read the file (assuming tab-delimited)
    data <- read.delim(file_path, sep = "\t")
    names(data)[names(data) == "junction"] = "unifiedJxnID"
    data = data[,c('unifiedJxnID','abundance_change', 'p_val', 'significant')]
    
    # Add data frame to the list with tissue name as key
    diego.output[[tissue.pair]] <- data
    
  }
  diego.output
}

# работает ли эта функция?
# сравнить с аутпутом dje
# что происходит с na в sajr??

processData = function(){
  output =list()
  for (tissue.pair in tissue.pairs){
    # merging tools outputs
    sajr.dje.common.jxns = merge(sajr.output[[tissue.pair]], dje.output[[tissue.pair]],
                                 by = "unifiedJxnID", all = TRUE) # inner join, duplicated rows will remain
    
    all.common.jxns = merge(sajr.dje.common.jxns, diego.output[[tissue.pair]],
                            by = "unifiedJxnID", all = TRUE) # inner join, duplicated rows will remain
    
    # -- finding significant
    # significant events
    sajr.all.significant.events = all.common.jxns[abs(all.common.jxns$dPSI)>=0.05 &
                                                    all.common.jxns$FDR.sajr <=0.2,]
    dje.all.significant.events = all.common.jxns[abs(all.common.jxns$logFC)>=1 &
                                                   all.common.jxns$FDR <=0.2,]
    diego.all.significant.events = all.common.jxns[all.common.jxns$significant=='yes',]
    
    # removing NA
    all.common.jxns =
      all.common.jxns[!is.na(all.common.jxns$logFC) &
                        !is.na(all.common.jxns$FDR) &
                        !is.na(all.common.jxns$dPSI) &
                        !is.na(all.common.jxns$FDR.sajr),]
    
    sajr.all.significant.events =
      sajr.all.significant.events[!is.na(sajr.all.significant.events$dPSI) &
                                    !is.na(sajr.all.significant.events$FDR.sajr),]
    dje.all.significant.events =
      dje.all.significant.events[!is.na(dje.all.significant.events$logFC) &
                                   !is.na(dje.all.significant.events$FDR),]
    
    diego.all.significant.events = na.omit(diego.all.significant.events)
    
    
    # searching for sinificant in both tools
    sign.all.tools = all.common.jxns[abs(all.common.jxns$logFC)>=1 &
                                       all.common.jxns$FDR <=0.2 &
                                       abs(all.common.jxns$dPSI)>=0.05 &
                                       all.common.jxns$FDR.sajr <=0.2 &
                                       all.common.jxns$significant=='yes',]
    
    sign.all.tools = sign.all.tools[!is.na(sign.all.tools$significant),]
    
    
    sajr.dje.significant.events = sajr.all.significant.events[(sajr.all.significant.events$unifiedJxnID %in%
                                                                 dje.all.significant.events$unifiedJxnID) &
                                                                !(sajr.all.significant.events$unifiedJxnID %in%
                                                                    sign.all.tools$unifiedJxnID),]
    
    sajr.diego.significant.events = sajr.all.significant.events[(sajr.all.significant.events$unifiedJxnID %in%
                                                                   diego.all.significant.events$unifiedJxnID) &
                                                                  !(sajr.all.significant.events$unifiedJxnID %in%
                                                                      sign.all.tools$unifiedJxnID),]
    
    dje.diego.significant.events =  dje.all.significant.events[(dje.all.significant.events$unifiedJxnID %in%
                                                                  diego.all.significant.events$unifiedJxnID) &
                                                                 !(dje.all.significant.events$unifiedJxnID %in%
                                                                     sign.all.tools$unifiedJxnID),]
    
    sajr.only.significant.events =
      sajr.all.significant.events[!(sajr.all.significant.events$unifiedJxnID %in%
                                      sajr.dje.significant.events$unifiedJxnID) &
                                    !(sajr.all.significant.events$unifiedJxnID %in%
                                        sajr.diego.significant.events$unifiedJxnID) &
                                    !(sajr.all.significant.events$unifiedJxnID %in% 
                                        sign.all.tools$unifiedJxnID),]
    
    dje.only.significant.events =
      dje.all.significant.events[!(dje.all.significant.events$unifiedJxnID %in%
                                     sajr.dje.significant.events$unifiedJxnID) &
                                   !(dje.all.significant.events$unifiedJxnID %in%
                                       sajr.diego.significant.events$unifiedJxnID) &
                                   !(dje.all.significant.events$unifiedJxnID %in% 
                                       sign.all.tools$unifiedJxnID),]
    
    diego.only.significant.events =
      diego.all.significant.events[!(diego.all.significant.events$unifiedJxnID %in%
                                       dje.diego.significant.events$unifiedJxnID) &
                                     !(diego.all.significant.events$unifiedJxnID %in%
                                         sajr.diego.significant.events$unifiedJxnID) &
                                     !(diego.all.significant.events$unifiedJxnID %in% 
                                         sign.all.tools$unifiedJxnID),]
    
    # ordering
    sign.all.tools =
      sign.all.tools[order(-abs(sign.all.tools$logFC),
                           -abs(sign.all.tools$dPSI),
                           sign.all.tools$FDR,
                           sign.all.tools$FDR.sajr), ]
    
    # sajr.only.significant.events =
    #   sajr.only.significant.events[order(-abs(sajr.only.significant.events$dPSI),
    #                                      sajr.only.significant.events$FDR.sajr), ]
    # dje.only.significant.events =
    #   dje.only.significant.events[order(-abs(dje.only.significant.events$logFC),
    #                                     dje.only.significant.events$FDR), ]
    
    # saving to the list
    output[[tissue.pair]]$all.common.jxns = all.common.jxns
    output[[tissue.pair]]$sajr.dje.significant.events = sajr.dje.significant.events
    output[[tissue.pair]]$sajr.diego.significant.events = sajr.diego.significant.events
    output[[tissue.pair]]$dje.diego.significant.events = dje.diego.significant.events
    
    output[[tissue.pair]]$sajr.only.significant.events = sajr.only.significant.events
    output[[tissue.pair]]$dje.only.significant.events = dje.only.significant.events
    output[[tissue.pair]]$diego.only.significant.events = diego.only.significant.events
    output[[tissue.pair]]$sign.all.tools = sign.all.tools
    
  }
  output
}

age.groups.ids.tissues = findPairsInfo()
jxns.ids.passed = filterJxns()
sajr.output = runSAJR()
dje.output = runDJE()
diego.output = readDiegoOutput()

output = processData()

# проверить что в output все работает как надо!

#output$Brain

#dje.output$Brain$dje.sig
# 83 sig
# сравнить с тем что фильтрую я
# что с sajr с сердцем, посмотреть на фильтрации, что с сердцем?

#original.dje.output = dje.output

#dje.output$Brain

#significant.events
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


