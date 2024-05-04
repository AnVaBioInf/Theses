source('plotCoverage.R')
#source('samples_ann_preprocessing.R')
source('tumor_ann.R')

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


set.colors = function(sign.jxns.info.dev.and.cancer, gene, covs){
  # creating a vector of colors assigned to significant junctions
  # unique junction = unique color
  all.sign.jxns = unique(sign.jxns.info.dev.and.cancer$unifiedJxnID)
  # Generate unique colors for each junction
  jxn.colors = rainbow(length(all.sign.jxns))
  names(jxn.colors) = all.sign.jxns
  
  # to set color to sign jxn in the tissue
  sign.gene.jxns.tissue.df = 
    sign.jxns.info.dev.and.cancer[sign.jxns.info.dev.and.cancer$tissue==tissue &
                                  sign.jxns.info.dev.and.cancer$GeneID==gene , , drop=F]
  
  # colors
  sign.jxn.col = jxn.colors[sign.gene.jxns.tissue.df$unifiedJxnID]
  print(sign.jxn.col)
  print(table(ifelse(sub(':.$', '', rownames(covs$juncs)) %in% names(sign.jxn.col),
         sign.jxn.col[sub(':.$', '', rownames(covs$juncs))],
         'blue')))
  ifelse(sub(':.$', '', rownames(covs$juncs)) %in% names(sign.jxn.col),
           sign.jxn.col[sub(':.$', '', rownames(covs$juncs))],
           'blue'
   )
}


# for every gene
for (gene in unique(sign.jxns.info.dev.and.cancer$GeneID)){
  # setting number of plot rows to number of tissues where selected gene junctions were significant, but no more than 3
  # selecting significant junctions for the GENE in both, development and cancer
  sign.gene.jxns.df = sign.jxns.info.dev.and.cancer[sign.jxns.info.dev.and.cancer$GeneID==gene,]
  par(mfrow = c(min(length(unique(sign.gene.jxns.df$tissue)),3),2), bty='n')
  
  # for every tissue where any of selected junctions are significant
  for (tissue in (unique(sign.gene.jxns.df$tissue))){ 
    print(tissue)
    print(gene)
    
    covs.summed.gene = get.covs(gene)
    fetus.covs.summed.gene =  covs.summed.gene[['fetus.covs.summed.gene']]
    adult.covs.summed.gene = covs.summed.gene[['adult.covs.summed.gene']]
    gene.region.coords = get.gene.region(gene)
    # same?
    jxn.colors.adult = set.colors(sign.jxns.info.dev.and.cancer, gene, adult.covs.summed.gene)
    jxn.colors.fetus = set.colors(sign.jxns.info.dev.and.cancer, gene, fetus.covs.summed.gene)
    
    plotReadCov(fetus.covs.summed.gene,
                junc.col = jxn.colors.fetus,
                xlim=gene.region.coords,
                plot.junc.only.within = F,
                min.junc.cov.f = 0.05,
                sub='Before birth',
                main=paste(tissue,gene)
                )
    plotReadCov(adult.covs.summed.gene,
                junc.col = jxn.colors.adult,
                xlim=gene.region.coords,
                plot.junc.only.within = F,
                min.junc.cov.f = 0.05,
                sub='After birth'
                )
  }
}
