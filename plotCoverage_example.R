source('plotCoverage.R')
#source('samples_ann_preprocessing.R')
source('tumor_ann.R')
#library(RColorBrewer)


# +
get.gene.region = function(gene){
  # selecting significant junctions for the GENE in both, development and cancer
  sign.gene.jxns.df = sign.jxns.info.dev.and.cancer[sign.jxns.info.dev.and.cancer$GeneID==gene,]
  
  # region of gene where significant junctions of gene are located
  sign.gene.jxns.coords = strsplit(sign.gene.jxns.df$junctionID,'[:-]')
  gene.region.coords = c(min(unlist(lapply(sign.gene.jxns.coords, function(jxn.coord) jxn.coord[2]))),
                         max(unlist(lapply(sign.gene.jxns.coords, function(jxn.coord) jxn.coord[3]))))
  gene.region.coords = as.integer(gene.region.coords)
  gene.region.coords
}


get.sample.ids = function(rse,tissue){
  # sample ids
  # gene annotation for a tissue
  ann.tissue = rse@colData[rse@colData$tissue==tissue,]
  # tissue samples ids
  all.samples.ids.tissue = rownames(ann.tissue)
  # adult and fetus sample ids
  adult.samples.ids = rownames(ann.tissue[ann.tissue$age_group=='adult',])
  list(adult.samples.ids=adult.samples.ids, all.samples.ids.tissue=all.samples.ids.tissue)
}


get.covs = function(gene){
  # importing rse file
  rse.ERP109002.jxn.cytosk.genes = readRDS('rse.ERP109002.jxn.cytosk.genes.rds', refhook = NULL)
  
  # gene.grange is needed by rtracklayer to filter bw files
  gene.grange = 
    rse.ERP109002.jxn.cytosk.genes@rowRanges[rse.ERP109002.jxn.cytosk.genes@rowRanges$gene_names==gene,]
  rse.gene = 
    rse.ERP109002.jxn.cytosk.genes[rse.ERP109002.jxn.cytosk.genes@rowRanges$gene_names==gene,]
  
  sample.ids = get.sample.ids(rse.gene,tissue)
  all.samples.ids = sample.ids[['all.samples.ids.tissue']]
  adult.samples.ids = sample.ids[['adult.samples.ids']]

  # -- coverages
  # covearge for each sample, output is a list of lists with read coverages, start:end positions, juncs df
  # gene.cov.all.samples.list - named list for each tissue sample
  gene.cov.all.samples.list = 
    lapply(all.samples.ids, function(samples.id){getRecountCov(samples.id, rse.gene, gene.grange)}) 
  # assigning elements of the list sample ids names
  names(gene.cov.all.samples.list) = all.samples.ids
  
  # --- merging
  # sum coverage in each condition
  fetus.covs.summed.gene = 
    sumCovs(gene.cov.all.samples.list[!(names(gene.cov.all.samples.list) %in% adult.samples.ids)])
  
  adult.covs.summed.gene =
    sumCovs(gene.cov.all.samples.list[adult.samples.ids])
  
  list(fetus.covs.summed.gene=fetus.covs.summed.gene, adult.covs.summed.gene=adult.covs.summed.gene)
}

# --
assign.colors = function(sign.jxns.info.dev.and.cancer){
  sign.colors = c("#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  not.sign.colors = "#377EB8"
  
  # creating a vector of colors assigned to significant junctions
  # unique junction = unique color
  all.sign.jxns = unique(sign.jxns.info.dev.and.cancer$unifiedJxnID)
  # Generate unique colors for each junction

  jxn.colors = sign.colors[1:length(all.sign.jxns)]
  names(jxn.colors) = all.sign.jxns
  jxn.colors
}

set.colors = function(sign.jxns.info.dev.and.cancer, jxn, covs){
  # to set color to sign jxn in the tissue
  jxn.colors = assign.colors(sign.jxns.info.dev.and.cancer)
  
  sign.gene.jxns.tissue.df = 
    sign.jxns.info.dev.and.cancer[sign.jxns.info.dev.and.cancer$tissue==tissue &
                                  sign.jxns.info.dev.and.cancer$unifiedJxnID==jxn , , drop=F]
  # colors
  sign.jxn.col = jxn.colors[jxn]
  print(sign.gene.jxns.tissue.df)
  
  cols = sign.jxn.col[sub(':.$', '', rownames(covs$juncs))]
  cols = ifelse(is.na(cols),"#377EB8", cols)
  covs$juncs$cols = cols
  print(table(unname(cols)))
  covs
}


# for every gene
for (jxn in unique(sign.jxns.info.dev.and.cancer$unifiedJxnID)){
  # setting number of plot rows to number of tissues where selected gene junctions were significant, but no more than 3
  # selecting significant junctions for the GENE in both, development and cancer
  sign.jxn.df = sign.jxns.info.dev.and.cancer[sign.jxns.info.dev.and.cancer$unifiedJxnID==jxn,]
  gene = unique(sign.jxn.df$GeneID)
  par(mfrow = c(min(length(unique(sign.jxn.df$tissue)),3),2), bty='n')
  print(c("here, here", jxn))
  
  # for every tissue where any of selected junctions are significant
  for (tissue in (unique(sign.jxn.df$tissue))){ 
    print(tissue)
    print(gene)
    print(jxn)
    
    covs.summed.gene = get.covs(gene)
    fetus.covs.summed.gene =  covs.summed.gene[['fetus.covs.summed.gene']]
    adult.covs.summed.gene = covs.summed.gene[['adult.covs.summed.gene']]
    gene.region.coords = strsplit(jxn,'[:-]')
    gene.region.coords = as.integer(c(gene.region.coords[[1]][2], gene.region.coords[[1]][3]))
    gene.region.coords = c(gene.region.coords[1], gene.region.coords[2])
    
    print(gene.region.coords)
    # same?
    fetus.covs.summed.gene = set.colors(sign.jxns.info.dev.and.cancer, jxn, fetus.covs.summed.gene)
    adult.covs.summed.gene = set.colors(sign.jxns.info.dev.and.cancer, jxn, adult.covs.summed.gene)
  
    
    sign.jxn.tissue = sign.jxn.df[sign.jxn.df$tissue==tissue,]
    
    
    plotReadCov(fetus.covs.summed.gene,
                junc.col = fetus.covs.summed.gene$cols,
                xlim=gene.region.coords,
                plot.junc.only.within = F,
                min.junc.cov.f = 0.05,
                sub='Before birth',
                main=paste(tissue,gene, jxn)
                )
    plotReadCov(adult.covs.summed.gene,
                junc.col = fetus.covs.summed.gene$cols,
                xlim=gene.region.coords,
                plot.junc.only.within = F,
                min.junc.cov.f = 0.05,
                sub='After birth',
                main=paste('dPSI.sajr=',round(sign.jxn.tissue$dPSI, digits=2), 
                           'FDR.sajr=', round(sign.jxn.tissue$dPSI, digits=2),
                           'logFC.dje=', round(sign.jxn.tissue$logFC, digits=2),
                           'FDR.dje=', round(sign.jxn.tissue$FDR, digits=2),
                          'sign.diego=', sign.jxn.tissue$significant
                ))
  }
}
