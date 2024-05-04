#devtools::install_github('iaaka/sajr')
library('SAJR')
source("create_sajr_object.R")
source("create_DJexpress_object.R")
#source("samples_ann_preprocessing.R")

sajr_to_djexpress_jxn_coord = function(sjar.pair.output.i) {
  strand_dict = c('+' = 1, '-' = 2, '*' = 0)
  sjar.pair.output.i$junctionID = rownames(sjar.pair.output.i)
  sjar.pair.output.i$junctionID = sub("(\\d{3,})(-)(\\d{3,})(:[-+*])(:).+", "\\1:\\3\\4", sjar.pair.output.i$junctionID)
  sjar.pair.output.i$junctionID = 
    sapply(sjar.pair.output.i$junctionID, function(x) sub('.$', strand_dict[sub('.*(?=.$)', '', x, perl=T)], x))
  sjar.pair.output.i
}

unique.tissues = unique(rse.ERP109002.jxn.cytosk.genes@colData$tissue)
tissue.pairs = combn(unique.tissues, 2, simplify = FALSE)
# running sajr for every pair of tissues

sajr.pairs.output = list()
for (pair in tissue.pairs){
  samples.pairs.annotation = 
    rse.ERP109002.jxn.cytosk.genes@colData[rse.ERP109002.jxn.cytosk.genes@colData$tissue %in% pair
                                           & rse.ERP109002.jxn.cytosk.genes@colData$age_group == 'adult',]
  
  tissue.1.sample.ids = rownames(samples.pairs.annotation[samples.pairs.annotation$tissue == pair[1],])
  tissue.2.sample.ids = rownames(samples.pairs.annotation[samples.pairs.annotation$tissue == pair[2],])
  
  mod = list(tissue=factor(samples.pairs.annotation$tissue))
  sajr.tissue.pairs = sajr[,samples.pairs.annotation$external_id]
  
  sajr.tissue.pairs = 
    sajr.tissue.pairs[apply(sajr.tissue.pairs$i+sajr.tissue.pairs$e>=10,1,sum)>=2 & apply(sajr.tissue.pairs$ir,1,sd,na.rm=TRUE) > 0,]
  
  #print(sajr.tissue.pairs$ir)
  #print(dim(sajr.tissue.pairs$seg))
  
  alt.glm = as.data.frame( fitSAGLM(sajr.tissue.pairs,terms(x ~ tissue, keep.order=T),mod,return.pv=T) )
  # calculating psi
  alt.glm$dpsi = apply(sajr.tissue.pairs$ir, 1, function(x) mean(x[tissue.1.sample.ids],na.rm=T)-mean(x[tissue.2.sample.ids],na.rm=T))
  alt.glm$logFC.sajr = apply(sajr.tissue.pairs$ir, 1, function(x) log2(mean(x[tissue.1.sample.ids],na.rm=T)/mean(x[tissue.2.sample.ids],na.rm=T)))
  
  alt.glm$fdr = p.adjust(alt.glm$tissue,m='BH')
  #mean(x[,tissue.1.sample.ids]-mean(x[,tissue.2.sample.ids]))
  sajr.pairs.output[[paste0(pair, collapse = "")]] = alt.glm
  

  #sajr.pairs.output[[paste0(pair, collapse = "")]] = alt.glm[!is.na(alt.glm$tissue) &
  #                                                            alt.glm$tissue <=0.05,]
  
}

dje.pairs.output = list()
# running DJexpress for every pair of tissues
for (pair in tissue.pairs){
  samples.pairs.annotation = 
    rse.ERP109002.jxn.cytosk.genes[,rse.ERP109002.jxn.cytosk.genes@colData$tissue %in% pair
                                           & rse.ERP109002.jxn.cytosk.genes@colData$age_group == 'adult']
  anlz.out = makeDJexpress(samples.pairs.annotation, reference.condition=pair[2], FDR.threshold=0.05, logFC.threashold=2)
  dje.pairs.output[[paste0(pair, collapse = "")]] = anlz.out$dje.out
}

# change junction ids in DJexpress to match sajr
sajr.pairs.output = lapply(sajr.pairs.output, function(x) sajr_to_djexpress_jxn_coord(x))

# Merge element-wise using Map
# only junction names

#length(unique(sajr.pairs.output$BrainHeart$junctionID[sajr.pairs.output$BrainHeart$junctionID %in% dje.pairs.output$BrainHeart$junctionID]))

#length(unique(dje.pairs.output$BrainHeart$junctionID))

#dje_j_id = dje.pairs.output$BrainHeart[dje.pairs.output$BrainHeart$junctionID %in% sajr.pairs.output$BrainHeart$junctionID,]$junctionID

 
                                                        

# filtering output
sajr.significant = lapply(sajr.pairs.output, function(x) {
  x=na.omit(x[,c('tissue','dpsi','logFC.sajr','fdr','junctionID')])
  x[abs(x$dpsi) >=0.1 &
    x$fdr<=0.05,]
  #x[x$tissue<=0.05,]
  })

djexpress.significant = lapply(dje.pairs.output, function(x) {
  x[x$FDR<=0.05 &
    abs(x$logFC) >=2,]  #? fdr?
})


# only for sajr
common.jxns.sign= Map(function(x, y) {
  y <- y[, c('logFC', 'P.Value', 'junctionID')]  # Select desired columns from y
  merge(x, y, by = "junctionID", all = FALSE)
}, sajr.significant, djexpress.significant)

# only for DJexpress1
unique.jxns.sajr.sign = Map(function(x, y) {
  y <- y[, c('logFC', 'P.Value', 'junctionID')]  # Select desired columns from y
  x[!x$junctionID %in% y$junctionID, ]  # Rows in x but not in y
}, sajr.significant, djexpress.significant)

unique.jxns.djexpress.sign = Map(function(x, y) {
  y <- y[, c('logFC', 'P.Value', 'junctionID')]  # Select desired columns from y
  y[!(y$junctionID %in% x$junctionID), ]  # Rows in y but not in x
}, sajr.significant, djexpress.significant)

common.jxns.all= Map(function(x, y) {
  y <- y[, c('logFC', 'P.Value', 'junctionID')]  # Select desired columns from y
  merge(x, y, by = "junctionID", all = FALSE)
}, sajr.pairs.output, dje.pairs.output)


# graphs
# significant only in djexpress
# significant only in sajr
# not significant in djexpress
# not dignificant in sajr
# significant in both-dj express and sajr

# stacked barplot
par(mar = c(2, 2, 2, 2))  # Set the margins to smaller values

par(mfrow = c(7, 4))

Map(function(all.dje, all.sajr, common.sign, dj.only.sign, sajr.only.sign, pair.name, common.all) {
  sales <- data.frame(
    software = c("DJexpress", "SAJR"),
    common.sign.jxns = c(nrow(common.sign)+1, nrow(common.sign)+1),
    #common.all = c(nrow(common.all)-nrow(common.sign)+1, nrow(common.all)-nrow(common.sign)+1),
    
    only.one.tool.sgn = c(nrow(dj.only.sign)+1, nrow(sajr.only.sign)+1),
    
    not.sgn = c((nrow(all.dje)-nrow(common.sign)-nrow(dj.only.sign)+1), 
                (nrow(all.sajr)-nrow(common.sign)-nrow(dj.only.sign)+1))
  )
  
  print(sales)
  # Create the vector of colors for each quarter
  colors <- c("green", "red", "grey")
  
  words <- unlist(strsplit(pair.name, split = "(?<=[a-z])(?=[A-Z])", perl = TRUE))
  # Join the words with a dash
  result <- paste(words, collapse = " - ")
  
  # Create the stack bar chart
  barplot(
    t(as.matrix(sales[, -1])), 
    main = result,
    xlab = "",
    ylab = "# jxns",
    col = colors,
    beside = F,
    log='y',
    ylim = c(0+1, max(nrow(all.sajr),nrow(all.dje))),
    names.arg = sales$software,  # Specify the labels for each subplot
    legend = FALSE  # Disable the legend creation within each subplot
  )
  
  # intersecting 
  # 
  
  # dpsi vs logFC
  plot(common.all$dpsi,
       common.all$logFC,
       xlim=c(-1,1),
       ylim=c(-4,4))
  
  # Highlight points belonging to "common"
  points(common.sign$dpsi,
         common.sign$logFC,
         col = "green",
         pch = 16)
  
  # Highlight points belonging to
  points(common.all[common.all$junctionID %in% sajr.only.sign$junctionID,]$dpsi,
         common.all[common.all$junctionID %in% sajr.only.sign$junctionID,]$logFC,
         col = "blue",
         pch = 16)
  
  # Highlight points belonging to
  points(common.all[common.all$junctionID %in% dj.only.sign$junctionID,]$dpsi,
         common.all[common.all$junctionID %in% dj.only.sign$junctionID,]$logFC,
         col = "violet",
         pch = 16)

  # logFC vs logFC
  plot(common.all$logFC.sajr,
       common.all$logFC,
       xlim=c(-4,4),
       ylim=c(-4,4))
  
  # Highlight points belonging to "common"
  points(common.sign$logFC.sajr,
         common.sign$logFC,
         col = "green",
         pch = 16)
  
  # Highlight points belonging to
  points(common.all[common.all$junctionID %in% sajr.only.sign$junctionID,]$logFC.sajr,
         common.all[common.all$junctionID %in% sajr.only.sign$junctionID,]$logFC,
         col = "blue",
         pch = 16)
  
  # Highlight points belonging to
  points(common.all[common.all$junctionID %in% dj.only.sign$junctionID,]$logFC.sajr,
         common.all[common.all$junctionID %in% dj.only.sign$junctionID,]$logFC,
         col = "violet",
         pch = 16)
  
  
  plot(common.all$tissue,
       common.all$P.Value, 
       log='xy')
  
  
  # Highlight points belonging to "common"
  points(common.sign$tissue,
         common.sign$P.Value,
         col = "green",
         pch = 16)
  
  # Highlight points belonging to
  points(common.all[common.all$junctionID %in% sajr.only.sign$junctionID,]$tissue,
         common.all[common.all$junctionID %in% sajr.only.sign$junctionID,]$P.Value,
         col = "blue",
         pch = 16)
  
  # Highlight points belonging to
  points(common.all[common.all$junctionID %in% dj.only.sign$junctionID,]$tissue,
         common.all[common.all$junctionID %in% dj.only.sign$junctionID,]$P.Value,
         col = "violet",
         pch = 16)
  
  
  },
  dje.pairs.output,
  sajr.pairs.output,
  common.jxns.sign, 
  unique.jxns.djexpress.sign,
  unique.jxns.sajr.sign, 
  names(sajr.pairs.output),
  common.jxns.all)





# Начать с возрастов



################################
# tutorial
################################

# data = loadSAData(ann.gff='example/a.gff',c(1,2))
# data = setSplSiteTypes(data,'example/a.gff')
# data.f = data[data$seg$type %in% c('ALT','INT') & data$seg$position %in% c('LAST','INTERNAL','FIRST') & apply(data$i+data$e>=10,1,sum)==2 & apply(data$ir,1,sd,na.rm=TRUE) > 0,]
# mod = list(f=factor(c('treatment','control')))
# data.f.glm = fitSAGLM(data.f,terms(x ~ f),mod,return.pv=T)
# data.f.glm



# сколько джанкшенов проходит фильтрацию?

# Список фич
# 

# Значимые список фич

# константные интроны? тот у которого левый и правый сайт вырезаются только одним способом 
# sajr смотрит только на альтернативные интроны


# Сделать возрастные изменения
# Что брать за референс в тканях?

# График
# общие сер
# красные значимые


# корреляция
# значимые только в одном синим
# Значимые только в другом
# подписатьь количество

# Взять одну за референс
# поставить порог - те которые нам интересны для экзаменац
# большая амплитуда в одном методе
# 

# графики покрытия для каждого из методов
# для каждого из случаев

# для sajr посмотреть logFC

# для графиков покрытия из рекаунта загрузить покрытия ридами (в формате виг)
# чтобы выявить артифакты

# подписать самые выдающиеся гены
# топ 5 из каждого метода  - подписать эти гены


# графики для 

# посмотреть документ!
# сравнить с раковыми (есть табличка)

# пятница в 4 по москве
# разместить на одной страницы - 1 строка для кажд ткани
