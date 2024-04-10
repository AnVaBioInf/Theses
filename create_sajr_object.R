source("data.R")

# functions ---------------------
makeSites = function(junctions.info){ # junxtions of a gene
  junctions.info = unique(as.data.frame(junctions.info)) # choosing only unique rows nothing changes
  if(nrow(junctions.info)==0) 
    return(NULL) # if there are no junctions, return null
  
  #--
  # adding extra columns
  junctions.info$rightmost = junctions.info$leftmost =  NA  # adding 2 columns for rightmost and leftmost coordinates of a junction
  junctions.info$id = rownames(junctions.info) # adding a column with coordinates of each junction for gene i ("chrX:71910743-71913070:+" "chrX:71910743-71958862:+")

  #--
  # dublicating df
  junctions.info.start.fixated = junctions.info.end.fixated = junctions.info  # duplicating jxns dataframe
  junctions.info.start.fixated$side = 'l' # adding a column 'side'
  junctions.info.end.fixated$side = 'r' 
  
  # in duplicated dfs formatting rownames
  rownames(junctions.info.start.fixated) = 
    paste0(rownames(junctions.info.start.fixated),':',junctions.info.start.fixated$gene_name,':l')
  rownames(junctions.info.end.fixated) = 
    paste0(rownames(junctions.info.end.fixated),':',junctions.info.end.fixated$gene_name,':r')

  #----------
  # forming a list, where key is a junction, and corresponding element is all junctions with same start/end coordinate
  junctions.same.coordinate.list = list() 
  
  # filling in rigtmost and leftmost columns and the list
  for(junction in 1:nrow(junctions.info)){ # number of rows in junctions.info
    # same start (including the junction!)
    # finding junctions with the start coordinate same to the junction
    junctions.same.start.indicators = junctions.info$start == junctions.info$start[junction]  # true/false vector for the junction
    # adding to the list [the junction coordinate] as a key, and all junctions that have the same start as the junction
    junctions.same.coordinate.list[[rownames(junctions.info.start.fixated)[junction]]] = 
      rownames(junctions.info)[junctions.same.start.indicators] # assigning junction.info rownames! They fill further be used to select rows from counts df. 
    # --
    # filling in leftmost column = start coordinate
    junctions.info.start.fixated$leftmost[junction]  = 
      unique(junctions.info$start[junctions.same.start.indicators]) 
    # searching for the rightmost coordinate amongst junctions
    junctions.info.start.fixated$rightmost[junction] = 
      max(junctions.info$end[junctions.same.start.indicators])

    # same end (including the junction!)
    junctions.same.end.indicators = junctions.info$end==junctions.info$end[junction]
    junctions.same.coordinate.list[[rownames(junctions.info.end.fixated)[junction]]] = 
      rownames(junctions.info)[junctions.same.end.indicators]
    
    # --
    junctions.info.end.fixated$leftmost[junction]  = 
      min(junctions.info$start[junctions.same.end.indicators])
    junctions.info.end.fixated$rightmost[junction] = 
      unique(junctions.info$end[junctions.same.end.indicators]) 
  }
  
  # concatenating dataframes
  junctions.info.start.end.fixated = rbind(junctions.info.start.fixated, junctions.info.end.fixated) 
  
  list(junctions.info.start.end.fixated = junctions.info.start.end.fixated, 
       junctions.same.coordinate.list = junctions.same.coordinate.list[rownames(junctions.info.start.end.fixated)]) # elements to junctions.same.coordinate.list were asigned l after r. However, in df junctions.info.start.end.fixated first there is information for all l, than for all r. To make it corresponding, we will reorder the list
}


makeSAJR = function(rse.ERP109002.jxn.cytosk.gene,min.cov=10){ # min.cov was set based on binomial dispersion
  junctions.info.gene = rse.ERP109002.jxn.cytosk.gene@rowRanges
  counts.gene = as.matrix(rse.ERP109002.jxn.cytosk.gene@assays@data$counts) 
  
  junctions.sites.list = makeSites(junctions.info.gene) # list, contaning jxns_s and js2j  
  if(is.null(junctions.sites.list))
    return(NULL)
  junctions.info.start.end.fixated = junctions.sites.list$junctions.info.start.end.fixated
  junctions.same.coordinate.list = junctions.sites.list$junctions.same.coordinate.list
  
  # inclusion ratio calculation
  # -- creating count matrixes
  
  # inclusion segments - reads theat map to the segment itself (cassette exon)
  # exclusion segments - reads that map to junction between upstream and downstream segments
  # all segments - reads that map to junction between upstream and downstream segments + reads that map to the junction
  
  junction.counts.inclusion = junction.counts.all.same.coordinate = matrix(0, nrow=nrow(junctions.info.start.end.fixated), ncol=ncol(counts.gene)) 
  colnames(junction.counts.inclusion) = colnames(junction.counts.all.same.coordinate) = colnames(counts.gene) # samples
  rownames(junction.counts.inclusion) = rownames(junction.counts.all.same.coordinate) = rownames(junctions.info.start.end.fixated)
  
  # -- filling in count matrixes
  for(junction in 1:nrow(junctions.info.start.end.fixated)){    
    # counts for each sample are stored in counts.gene. New matrix's columns are asigned names same to colnames(counts.gene). 
    # when filling in matrix, counts are taken from counts.gene df, so counts for samples in matrix and df counts.gene are filled in correctly. 
    
    # selecting all raw of counts for junctions (number of inclusion)
    junction.counts.inclusion[junction,] = counts.gene[junctions.info.start.end.fixated$id[junction], ] # id column contains same junction ids as in count matrix
    
    junction.counts.all.same.coordinate[junction,] = apply(counts.gene[junctions.same.coordinate.list[[junction]], , drop=F], 2, sum)
    # drop=F prevents R from converting single column to vector. But here is not a single column?
    # rows from count.gene df are selected, and we sum counts for excluded (alternative) junctions for each column (sample) - 2 (apply to column), sum (function to apply)
  }
  
  # inclusion ratio - the proportion of transcripts that contains a segment
  junction.counts.inclusion.ratio = junction.counts.inclusion/junction.counts.all.same.coordinate # matrix
  junction.counts.inclusion.ratio[junction.counts.all.same.coordinate < min.cov] = NA # frequency of inclusion is not defined (min.cov = 10, set based on binomial distribution dispersion)
  
  # exclusion counts. Substract included junctions from all
  junction.counts.exclusion = junction.counts.all.same.coordinate - junction.counts.inclusion 
  
  sajr.gene = list(seg=junctions.info.start.end.fixated, i=junction.counts.inclusion, e=junction.counts.exclusion, ir=junction.counts.inclusion.ratio) 
  class(sajr.gene)='sajr'
  sajr.gene
}

# почитать про psi!

# -------------------------------
sajr = NULL

# for each cytoskeleton gene
for(gene in 1:nrow(genes.annotaion)){ 
  # gene id
  gene.id = rownames(genes.annotaion)[gene] # 
  # selecting information only for the gene
  rse.ERP109002.jxn.cytosk.gene = rse.ERP109002.jxn.cytosk.genes[rse.ERP109002.jxn.cytosk.genes@rowRanges$gene_id==gene.id,]  
  # making sajr object for the geme
  sajr.gene =  makeSAJR(rse.ERP109002.jxn.cytosk.gene) 

  if(is.null(sajr.gene))
    next
  
  # those columns already exist in seg
  # adding columns to seg (=junctions.info.start.end.fixated)
  #sajr.gene$seg$gene_id = gene.id # what's gene id????????????????????????????????????????
  #sajr.gene$seg$gene_names = genes.annotaion$gene_name[gene] # name of the gene
  
  # making combined sajr object for all genes
  if(is.null(sajr))
    sajr = sajr.gene
  else{
    for(element in names(sajr)){ # combining seg with seg, e with e, etc for each element, for every gene
      sajr[[element]] = rbind(sajr[[element]],sajr.gene[[element]])
      
    }
  }
}

# нет gene names!
sajr$seg


