library(recount3)

#source('src/rcount.as.utils.R')


# load counts ######
human_projects <- recount3::available_projects()
proj_info <- subset(
  human_projects,
  project == "ERP109002" & project_type == "data_sources"
)

rse.ERP109002 = create_rse(proj_info)  # genes
rse.ERP109002.jxn <- create_rse(proj_info,
                                type="jxn") # junctions

# write R object to a file
#saveRDS(rse.ERP109002,'rse.ERP109002.rds')
#saveRDS(rse.ERP109002.jxn,'rse.ERP109002.jxn.rds') 

genes = rbind(data.frame(gene_name=c('CYFIP1','CYFIP2','NCKAP1','NCKAP1L','ABI1','ABI2','ABI3','WASF1','WASF2','WASF3','BRK1'),group='WAVE'),
              data.frame(gene_name=c('NHS','NHSL1','NHSL2','KIAA1522'),group='NHS'),
              data.frame(gene_name=c('ARPC1A','ARPC1B','ARPC2','ARPC3','ARPC4','ARPC5','ACTR2','ACTR3','ACTR3B'),group='Arp2/3'))

# leaving only citoskeleton genes in rse.ERP109002
genes_annotaion = as.data.frame(rse.ERP109002@rowRanges[rse.ERP109002@rowRanges$gene_name %in% genes$gene_name,])

# adding group names
genes_annotaion$group = genes$group[match(genes_annotaion$gene_name,genes$gene_name)]

# leaving only genes on main chr
genes_annotaion = genes_annotaion[startsWith(as.character(genes_annotaion$seqnames),'chr'),]
#saveRDS(genes,'rds/sel.genes.rds')

#-----------------------------------------

# junxtion "annotation"
# intersecting ranges
overlaps = findOverlaps(query = rse.ERP109002.jxn@rowRanges, 
                        subject = rse.ERP109002@rowRanges[genes_annotaion$gene_id], 
                        type='any') # or within?
# посмотреть на те которы не входят в within?

#?
#rse.ERP109002.jxn@rowRanges[64106]
#rse.ERP109002@rowRanges[genes_annotaion$gene_id][1]
# junction is larger than gene?

rse.ERP109002.jxn.cytosk = rse.ERP109002.jxn[overlaps@from,]  
rse.ERP109002.jxn.cytosk@rowRanges$gene_id = genes_annotaion$gene_id[overlaps@to]
rse.ERP109002.jxn.cytosk@rowRanges$gene_name = genes_annotaion$gene_name[overlaps@to]


# functions ---------------------
makeSites = function(jxns){ # junxtions of a gene
  jxns = unique(as.data.frame(jxns)) # chosing only unique rows nothing changes
  if(nrow(jxns)==0) # nrow returns number of rows (junxtions)
    return(NULL) # if there are no junctions, return null
  
  #--
  # merging junctions by leftmost and rightmost coordinate
  jxns$rightmost = jxns$leftmost =  NA  # adding 2 columns for rightmost and leftmost coordinates of a junction
  jxns$jxn.id = rownames(jxns) # adding a column with coordinates of each junction for gene i ("chrX:71910743-71913070:+" "chrX:71910743-71958862:+")

  #--
  jxns_l = jxns_r = jxns  # duplicating jxns
  jxns_l$side = 'l' # adding a column 'side'
  jxns_r$side = 'r' # start and end coordinates of a junction
  rownames(jxns_l) = paste0(rownames(jxns_l),':',jxns_l$gene_name,':l')
  rownames(jxns_r) = paste0(rownames(jxns_r),':',jxns_r$gene_name,':r')

  #----------
  js2j = list() #?
  for(i in 1:nrow(jxns_l)){ # number of rows in jl (which is equal to nrow(j))
    # finding same start coordinates
    jxns_same_start_tf = jxns$start==jxns$start[i]  # true/false vector  for junction i
    
    # l - общее начало у джанкшенов
    # list js2s: junction = all junctions with the same start/end junction
    # берем имена строк jxns_l - они были изменены! (добавлено имя гена итд)
    js2j[[rownames(jxns_l)[i]]] = rownames(jxns)[jxns_same_start_tf] # listing all jxn with the same start coord
    # добавляем самый левый старт (они же должны ьыть одинаковые?)
    jxns_l$leftmost[i]  = min(jxns$start[jxns_same_start_tf]) # the lefmost start coordinate, but all coordinates are the same, so it doesn't do anything
    # добавляем самый последний конец соответсвующий этому началу
    jxns_l$rightmost[i] = max(jxns$end[jxns_same_start_tf])

    # r - общий конец
    # right (doing the same thing for the end coordinate)
    jxns_same_end_tf = jxns$end==jxns$end[i]
    js2j[[rownames(jxns_r)[i]]] = rownames(jxns)[jxns_same_end_tf]
    # добавляем самое первое начало соответствующий этому концу?? или
    jxns_r$leftmost[i]  = min(jxns$start[jxns_same_end_tf])
    jxns_r$rightmost[i] = max(jxns$end[jxns_same_end_tf])  # end rightmost coordinate
  }
  jxns_s = rbind(jxns_l,jxns_r) # binding dfs
  list(jxns_s=jxns_s, js2j=js2j[rownames(jxns_s)]) # что такое rownames(..)?
  # jxns_s - junction ... leftmost and rightmost isoforms ?
  # js2j - jxns_s to jsnx - list of junctions with the same start for every junction
}

# for every gene separately finding starts and ends (they can be on different chromosomes and start/end can intersect)

makeSAJR = function(jxns,min.cov=10){ # min.cov=10? filtering?
  js = makeSites(jxns@rowRanges) # list, contaning jxns_s and js2j  

  if(is.null(js))
    return(NULL)

  counts = as.matrix(jxns@assays@data$counts) 
  
  # i - intron
  # e - exon
  # ir - intron retention
  i = e = matrix(0, nrow=nrow(js$jxns_s), ncol=ncol(counts)) # making 2 matrixes ??
  colnames(i) = colnames(e) = colnames(counts) # samples
  rownames(i) = rownames(e) = rownames(js$jxns_s) # ??????????????????????????????????
  
  for(j in 1:nrow(js$jxns_s)){    # заполняем матрицу
    i[j,] = counts[js$jxns_s$jxn.id[j], ] # из матрицы каунтов выбираем
    e[j,] = apply(counts[js$js2j[[j]], , drop=F], 2, sum) # сумма ..?
  }
  
  ir = i/e   # intron retention
  ir[e<min.cov] = NA
  e = e - i 
  r = list(seg=js$jxns_s, i=i, e=e, ir=ir)   # seg here
  class(r)='sajr'
  r
}


# -------------------------------
sajr = NULL


# на вход даем гены! это правильно?
for(i in 1:nrow(genes_annotaion)){ # for each gene
  gene_id=rownames(genes_annotaion)[i] # gene id
  jxns = rse.ERP109002.jxn.cytosk[rse.ERP109002.jxn.cytosk@rowRanges$gene_id==gene_id,] # RangedSummarizedExperiment for a gene 

  # t?
  t =  makeSAJR(jxns) # passing rse object for gene i to function makeSAJR

  if(is.null(t))
    next
  
  t$seg$gene_id = gene_id
  t$seg$gene_names = genes_annotaion$gene_name[i]
  if(is.null(sajr))
    sajr = t

  else{
    for(j in names(sajr)){
      sajr[[j]] = rbind(sajr[[j]],t[[j]])
    }
  }
}

sajr

