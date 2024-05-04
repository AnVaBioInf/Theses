library(recount3)
# 
# # load counts ######
# human_projects <- recount3::available_projects()
# proj_info <- subset(
#   human_projects,
#   project == "ERP109002" & project_type == "data_sources"
# )
# 
# rse.ERP109002 = create_rse(proj_info)  # genes
# rse.ERP109002.jxn <- create_rse(proj_info,
#                                 type="jxn") # junctions
# 
# # write R object to a file
# saveRDS(rse.ERP109002,'rse.ERP109002.rds')
# saveRDS(rse.ERP109002.jxn,'rse.ERP109002.jxn.rds')
# 
# genes = rbind(data.frame(gene_name=c('CYFIP1','CYFIP2','NCKAP1','NCKAP1L','ABI1','ABI2','ABI3','WASF1','WASF2','WASF3','BRK1'),group='WAVE'),
#               data.frame(gene_name=c('NHS','NHSL1','NHSL2','KIAA1522'),group='NHS'),
#               data.frame(gene_name=c('ARPC1A','ARPC1B','ARPC2','ARPC3','ARPC4','ARPC5','ACTR2','ACTR3','ACTR3B'),group='Arp2/3'))
# 
# # leaving only citoskeleton genes in rse.ERP109002
# genes.annotaion = as.data.frame(rse.ERP109002@rowRanges[rse.ERP109002@rowRanges$gene_name %in% genes$gene_name,])
# 
# # adding group names
# genes.annotaion$group = genes$group[match(genes.annotaion$gene_name,genes$gene_name)]
# 
# # leaving only genes on main chr
# genes.annotaion = genes.annotaion[startsWith(as.character(genes.annotaion$seqnames),'chr'),]
# 
# # sorting
# genes.annotaion = genes.annotaion[order(genes.annotaion$seqnames, genes.annotaion$start), ]
# 
# saveRDS(genes.annotaion,'genes.annotaion.rds')
# 
# # -----
# 
rse.ERP109002 = readRDS('rse.ERP109002.rds', refhook = NULL)
rse.ERP109002.jxn = readRDS('rse.ERP109002.jxn.rds', refhook = NULL)
genes.annotaion = readRDS('genes.annotaion.rds', refhook = NULL)
# 
# #-----------------------------------------
# 
# # --junxtion "annotation"
# # intersecting ranges
# overlaps = findOverlaps(query = rse.ERP109002.jxn@rowRanges, 
#                         subject = rse.ERP109002@rowRanges[genes.annotaion$gene_id], 
#                         type='any') # or within?
# 
# #rse.ERP109002.jxn@rowRanges[64106]
# #rse.ERP109002@rowRanges[genes_annotaion$gene_id][1]
# # junction is larger than gene?
# 
# rse.ERP109002.jxn.cytosk.genes = rse.ERP109002.jxn[overlaps@from,]  
# rse.ERP109002.jxn.cytosk.genes@rowRanges$gene_id = genes.annotaion$gene_id[overlaps@to]
# rse.ERP109002.jxn.cytosk.genes@rowRanges$gene_names = genes.annotaion$gene_name[overlaps@to]
# 
# #?????? перенести этот код до оверлепинга!? Что за ворнинг?
# # -- removing duplicates
# counts = as.matrix(rse.ERP109002.jxn.cytosk.genes@assays@data$counts) 
# # duplicated junctions
# rownames(counts)[duplicated(rownames(counts))]
# table(duplicated(rownames(counts)))
# # adding column with junction ids (to make sure that it finds duplicates only for the same junctions)
# counts = cbind(ids=rownames(counts), counts)
# # dropping duplicated rows
# rse.ERP109002.jxn.cytosk.genes = rse.ERP109002.jxn.cytosk.genes[which(!duplicated(counts)),]
# 
# 
# 
# rse.ERP109002.jxn.cytosk.genes = saveRDS(rse.ERP109002.jxn.cytosk.genes,'rse.ERP109002.jxn.cytosk.genes.rds')
rse.ERP109002.jxn.cytosk.genes = readRDS('rse.ERP109002.jxn.cytosk.genes.rds', refhook = NULL)

# если что то сломается- искать здесь!!!!
#rm(overlaps, counts, rse.ERP109002, rse.ERP109002.jxn, genes.annotaion)

# checking
#counts = as.matrix(rse.ERP109002.jxn.cytosk.genes@assays@data$counts) 
#table(duplicated(rownames(counts)))


# preprocessing annotation
#source('samples_ann_preprocessing.R')

# почему не сохранить все в rse.ERP109002.jxn.cytosk.genes.rse?