### In this file project rse data are downloaded, filtered by genes (only cytoskeleton genes are left)
### and gene names are added to rse@rowRanges 

# library(recount3)
# 
# # -- loading project rse
# human_projects <- recount3::available_projects()
# proj_info <- subset(
#   human_projects,
#   project == "ERP109002" & project_type == "data_sources"
# )
# rse.gene = create_rse(proj_info) 
# rse.jxn <- create_rse(proj_info,
#                                 type="jxn") # junctions
# 
# # -- gene rse only with cytoskeleton genes --
# # cytoskeleton genes
# cytoskeleton.genes = rbind(data.frame(gene_name=c('CYFIP1','CYFIP2','NCKAP1','NCKAP1L','ABI1','ABI2','ABI3','WASF1','WASF2','WASF3','BRK1'),group='WAVE'),
#                            data.frame(gene_name=c('NHS','NHSL1','NHSL2','KIAA1522'),group='NHS'),
#                            data.frame(gene_name=c('ARPC1A','ARPC1B','ARPC2','ARPC3','ARPC4','ARPC5','ACTR2','ACTR3','ACTR3B'),group='Arp2/3'))
# # 24 genes
# 
# rse.gene.cytosk = rse.gene[rse.gene@rowRanges$gene_name %in% cytoskeleton.genes$gene_name,]
# 
# # adding group names
# rse.gene.cytosk@rowRanges$group = cytoskeleton.genes$group[
#   match(rse.gene.cytosk@rowRanges$gene_name,cytoskeleton.genes$gene_name)]
# 
# # leaving only genes on main chr
# rse.gene.cytosk@rowRanges = 
#   rse.gene.cytosk@rowRanges[startsWith(as.character(rse.gene.cytosk@rowRanges@seqnames),'chr'),]
# 
# # sorting
# rse.gene.cytosk@rowRanges = 
#   rse.gene.cytosk@rowRanges[order(rse.gene.cytosk@rowRanges@seqnames, rse.gene.cytosk@rowRanges@ranges), ]
# 
# # subsetting seqinfo 
# new_seqinfo <- Seqinfo(seqnames = as.character(unique(rse.gene.cytosk@rowRanges@seqnames)),
#                        seqlengths = NA,
#                        isCircular = NA)
# seqlevels(rse.gene.cytosk@rowRanges) <- seqlevels(new_seqinfo)
# 
# # Assign the new Seqinfo to the GRanges object
# seqinfo(rse.gene.cytosk@rowRanges) <- new_seqinfo
# 
# # --making junxtion "annotation"
# # intersecting ranges
# # selecting only cytoskeleton genes in rse.jxn and adding gene column to it (making jxn "annotation")
# # The function compares each interval in the query object with the intervals in the subject object and identifies any overlaps.
# overlaps = findOverlaps(query = rse.jxn@rowRanges,
#                         subject = rse.gene.cytosk@rowRanges,
#                         type='within') # the start and end positions of the query range must fall within the start and end positions of the subject range.
# 
# # check
# rse.jxn@rowRanges[64374]
# rse.gene.cytosk@rowRanges[1]
# 
# # selecting only jxns of cytoskeleton genes and assing gene info
# rse.jxn.cytosk = rse.jxn[overlaps@from,]
# rse.jxn.cytosk@rowRanges$gene_id = rse.gene.cytosk@rowRanges$gene_id[overlaps@to]
# rse.jxn.cytosk@rowRanges$gene_names = rse.gene.cytosk$gene_name[overlaps@to]
# 
# # -- removing duplicates
# # checking
# table(duplicated(rse.jxn.cytosk@rowRanges@ranges@NAMES))
# nms = unique(rse.jxn.cytosk@rowRanges@ranges@NAMES)
# rse.jxn.cytosk = rse.jxn.cytosk[nms,]
# 
# # removing unneeded data
# rm(overlaps, rse.gene, rse.jxn)
# 
# # -- saving files
# saveRDS(rse.gene.cytosk,'rse.gene.cytosk.rds')
# saveRDS(rse.jxn.cytosk,'rse.jxn.cytosk.rds')

# -- reading files
rse.gene.cytosk = readRDS('rse.gene.cytosk.rds', refhook = NULL)
rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)


