# the project goal is to search for AS events in development
# DIEGO, DJexpress and SAJR tools are used for that (work with junction counts)
# in this file functions to make input for those tools are gathered
library(SAJR)
library(DJExpress)
source('rseAnnotationPreprocessing.R')

rse2countDf = function(rse){
  counts = as.matrix(assay(rse, "counts"))
  counts = as.data.frame(counts)
  counts
}

# проверить!
filterRse = function(rse, tissue, mim_row_sum=10, min_variance=0){
  rse = rse[,rse@colData$tissue==tissue]
  rse = rse[apply(rse@assays@data$counts, 1, sum) >= mim_row_sum &
            apply(rse@assays@data$counts, 1, var) >= min_variance, ]
  rse
}

findConditionIds = function(rse_filtered, age_group){
  rownames(rse_filtered@colData[rse_filtered@colData$age_group == age_group,])
}

# 
# #=================================DIEGO=======================================
# # http://legacy.bioinf.uni-leipzig.de/Software/DIEGO/
# # -- a file (table of splice junction supports per sample)
# makeAFile = function(rse.filtered, tissue, path_input){
#   junction_table = as.matrix(assay(rse.filtered, "counts"))
#   junction_table = as.data.frame(junction_table)
#   sample_ids = colnames(junction_table)
#   # location of the splice junction or exon, which must have the string:number-number format
#   junction_table$junction = sub('(:[-+*])', '' , rownames(junction_table))
#   #  the type of splice junction
#   junction_table$type = 'N_w'
#   junction_table$geneID = rse.filtered@rowRanges$gene_id
#   junction_table$geneName = rse.filtered@rowRanges$gene_name
#   junction_table = junction_table[,c('junction', 'type', sample_ids, 'geneID', 'geneName')]
#   junction_table = junction_table[order(junction_table$geneID),]
#   write.table(junction_table,
#               paste0(path_input, "/junction_table_",tissue, ".txt"),
#               sep = "\t",
#               row.names = FALSE, col.names = TRUE, quote = FALSE)
# }
# 
# # --b file (condition to sample relation in the format: condition tab-delimiter sampleName)
# makeBFile = function(rse.filtered, tissue, path_input){
#   group_table = rse.filtered@colData[,'age_group',drop=F]
#   group_table$sample_id = rownames(group_table)
#   group_table = group_table[order(group_table$age_group), ]
#   write.table(group_table,
#               paste0(path_input, "/group_table_", tissue, ".txt"),
#               sep = "\t",
#               row.names = FALSE, col.names = FALSE, quote = FALSE)
# }
# 
# # ++++++++++++++++++++++++++++++
# makeDiegoBashFile = function(tissue,
#                              path_input, path_output,
#                              reference_condition='adult', 
#                              min_support = 1, # minimum jxn count for a splice site to be considered.
#                              min_samples = 1, #  minimum number of samples that must show the minimum support
#                              FDR_threashold = 0.05, # adjusted p-value threshold
#                              fold_change = 1 # ? ratio of read counts of a splice junction in one condition compared to another
#                              ){
#   # Define the Bash script content
#   diego_bash_file_content = paste0(
#   '#!/bin/bash
#   python DIEGO/diego.py \\\
#   -a ', path_input, '/junction_table_', tissue, '.txt \\\
#   -b ', path_input, '/group_table_', tissue, '.txt \\\
#   -x ', reference_condition, ' \\\
#   --minsupp ', min_support,' \\\
#   --minsamples ', min_samples,' \\\
#   --significanceThreshold ', FDR_threashold,' \\\
#   --foldchangeThreshold ', fold_change,' \\\
#   > ', path_output, '/DIEGO_output_', tissue, '.txt')
#   filename = paste0(path_input,'/run_diego_', tissue,'.sh')
#   writeLines(diego_bash_file_content, filename)
#   system(paste("chmod +x", filename))
# }
# 
# makeDiegoInputFiles = function(rse, tissue,
#                               path_input = '/home/an/DIEGO_input_files',
#                               path_output = '/home/an/DIEGO_output_files'){
#   rse.filtered = filterRse(rse, tissue)
#   makeAFile(rse.filtered, tissue, path_input)
#   makeBFile(rse.filtered, tissue, path_input)
#   makeDiegoBashFile(tissue, path_input, path_output)
# }
# 
# readDiegoOutput = function(path_output = '/home/an/DIEGO_output_files', tissue){
#   file_path = paste0(path_output, '/DIEGO_output_', tissue, '.txt')  
#   diego_output = read.delim(file_path, sep = "\t")
#   diego_output
# }
# 
# 
# # to run DIEGO in terminal
# # conda create -n DIEGO_1 numpy=1.9 scipy matplotlib
# # Packages were reinstalled because -e (for drawing dendrograms) wasn't working (numpy 1.9 installed instead, and than other packages reinstalled, in the above code I specified numpy version, idk if it will help)
# # conda activate DIEGO_1
# # cp /home/an/DIEGO_input_files/a.input.file /home/an/DIEGO_input_files/b_file /home/an/anaconda3/envs/DIEGO_1
# # wget http://legacy.bioinf.uni-leipzig.de/Software/DIEGO/DIEGO.tar.gz
# # tar -xzf DIEGO.tar.gz
# # python DIEGO/diego.py -a a.input.file -b b_file -x fetus --minsupp 1 -d 1 -q 1.0 -z 1.0  > DIEGO_output
# 
# 
# 
# #===================================DJexpress==================================
# # replace +/- strand with 1/2 (0: undefined, 1: +, 2: -)
# # STAR manual, p.12 on output splice junctions file
# makeDjeCoordinates <- function(coordinates_vector) {
#   strand_to_numb_dict <- c("+" = 1, "-" = 2, "*" = 0)
#   strand_signs <- sub(".*(?=.$)", "", coordinates_vector, perl = TRUE)
#   coordinates_vector =  sub("([0-9]+)(-)([0-9]+)", "\\1:\\3", coordinates_vector)
#   coordinates_vector = sapply(seq_along(coordinates_vector), function(i) {
#                                 sub(".$", strand_to_numb_dict[strand_signs[i]], coordinates_vector[i])
#                               })
#   coordinates_vector
# }
# 
# # instead of DJEimport() etc, because input data differ (recount3 jxns instead of STAR raw out file)
# makePrepOutObj = function(filtered_rse, tissue, reference_condition){
#   JunctExprfilt = rse2countDf(rse_filtered)
#   rownames(JunctExprfilt) = makeDjeCoordinates(rownames(JunctExprfilt))
#   featureID = rownames(JunctExprfilt)
#   groupID = rse_filtered@rowRanges$gene_name
#   age_group_factor = relevel(as.factor(rse_filtered@colData$age_group), 
#                               ref = reference_condition)
#   design = model.matrix(~age_group_factor)
#   list(JunctExprfilt=JunctExprfilt, featureID=featureID, groupID=groupID, design=design)
# }
# 
# 
# runDJExpress = function(rse, tissue, reference_condition='adult',
#                         FDR_threshold=0, logFC_threashold=0){
#   filtered_rse = filterRse(rse, tissue)
#   prep_out = makePrepOutObj(filtered_rse, tissue, reference_condition)
#   reference_sample_ids = findConditionIds(filtered_rse, reference_condition)
#   
#   anlz_out <- DJEanalyze(prepare.out = prep_out,
#                          Group1 = reference_sample_ids,
#                          FDR = FDR_threshold,
#                          logFC = logFC_threashold)
#   anlz_out
# }


# ================================SAJR=========================================
makeSites = function(gene.jxns.info){ # junxtions of a gene
  gene.jxns.info = unique(as.data.frame(gene.jxns.info)) # choosing only unique rows nothing changes
  if(nrow(gene.jxns.info)==0)
    return(NULL) # if there are no junctions, return null

  # adding extra columns
  gene.jxns.info$rightmost = gene.jxns.info$leftmost =  NA  # adding 2 columns for rightmost and leftmost coordinates of a junction
  gene.jxns.info$id = rownames(gene.jxns.info) # adding a column with coordinates of each junction for gene i ("chrX:71910743-71913070:+" "chrX:71910743-71958862:+")

  # dublicating df
  jxns.same.start.info = jxns.same.end.info = gene.jxns.info  # duplicating jxns dataframe
  jxns.same.start.info$side = 'l' # adding a column 'side'
  jxns.same.end.info$side = 'r'

  # in duplicated dfs formatting rownames
  rownames(jxns.same.start.info) = paste0(rownames(jxns.same.start.info),':',jxns.same.start.info$gene_name,':l')
  rownames(jxns.same.end.info) = paste0(rownames(jxns.same.end.info),':',jxns.same.end.info$gene_name,':r')

  # forming a list, where key is a junction, and corresponding element is all junctions with same start/end coordinate
  jxns.same.coordinate.list = list()
  # filling in rigtmost and leftmost columns and the list
  for(junction in 1:nrow(gene.jxns.info)){ # number of rows in gene.jxns.info
    # same start (including the junction!)
    # finding junctions with the start coordinate same to the junction
    jxns.same.start.tf = gene.jxns.info$start == gene.jxns.info$start[junction]  # true/false vector for the junction
    # adding to the list [the junction coordinate] as a key, and all junctions that have the same start as the junction
    jxns.same.coordinate.list[[rownames(jxns.same.start.info)[junction]]] =
      rownames(gene.jxns.info)[jxns.same.start.tf] # assigning junction.info rownames! They fill further be used to select rows from counts df.

    # filling in leftmost column = start coordinate
    jxns.same.start.info$leftmost[junction]  = unique(gene.jxns.info$start[jxns.same.start.tf])
    # searching for the rightmost coordinate amongst junctions
    jxns.same.start.info$rightmost[junction] =  max(gene.jxns.info$end[jxns.same.start.tf])

    # same end (including the junction!)
    jxns.same.end.tf = gene.jxns.info$end == gene.jxns.info$end[junction]
    jxns.same.coordinate.list[[rownames(jxns.same.end.info)[junction]]] =
      rownames(gene.jxns.info)[jxns.same.end.tf]

    jxns.same.end.info$leftmost[junction]  = min(gene.jxns.info$start[jxns.same.end.tf])
    jxns.same.end.info$rightmost[junction] = unique(gene.jxns.info$end[jxns.same.end.tf])
  }

  # concatenating dataframes
  jxns.same.start.or.end.info = rbind(jxns.same.start.info, jxns.same.end.info)

  list(jxns.same.start.or.end.info = jxns.same.start.or.end.info,
       jxns.same.coordinate.list = jxns.same.coordinate.list[rownames(jxns.same.start.or.end.info)]) # elements to jxns.same.coordinate.list were asigned l after r. However, in df junctions.info.start.end.fixated first there is information for all l, than for all r. To make it corresponding, we will reorder the list
}


makeSAJRgene = function(rse.gene.filtered, min.cov=10){ # min.cov was set based on binomial dispersion
    jxns.sites.list = makeSites(rse.gene.filtered@rowRanges) # list, contaning jxns_s and js2j
    gene.jxns.counts = as.matrix(rse.gene.filtered@assays@data$counts)
    
    if(is.null(jxns.sites.list))
      return(NULL)
    jxns.same.start.or.end.info = jxns.sites.list$jxns.same.start.or.end.info
    jxns.same.coordinate.list = jxns.sites.list$jxns.same.coordinate.list
    
    # inclusion ratio calculation
    # -- creating count matrixes
  
    # inclusion segments - reads that map to the segment itself (cassette exon)
    # exclusion segments - reads that map to junction between upstream and downstream segments
    # all segments - reads that map to junction between upstream and downstream segments + reads that map to the junction
    jxns.incl.count.mtrx = jxns.incl.and.excl.count.mtrx = 
      matrix(0, nrow=nrow(jxns.same.start.or.end.info), ncol=ncol(gene.jxns.counts))
    colnames(jxns.incl.count.mtrx) = colnames(jxns.incl.and.excl.count.mtrx) = colnames(gene.jxns.counts) # samples
    rownames(jxns.incl.count.mtrx) = rownames(jxns.incl.and.excl.count.mtrx) = rownames(jxns.same.start.or.end.info)
  
    # -- filling in count matrixes
    for(junction in 1:nrow(jxns.same.start.or.end.info)){
      # counts for each sample are stored in gene_jxns_counts. New matrix's columns are asigned names same to colnames(gene_jxns_counts).
      # when filling in matrix, counts are taken from gene_jxns_counts df, so counts for samples in matrix and df gene_jxns_counts are filled in correctly.
  
      # selecting all raw of counts for junctions (number of inclusion)
      jxns.incl.count.mtrx[junction,] = gene.jxns.counts[jxns.same.start.or.end.info$id[junction], ] # id column contains same junction ids as in count matrix
      jxns.incl.and.excl.count.mtrx[junction,] = 
        apply(gene.jxns.counts[jxns.same.coordinate.list[[junction]], , drop=F], 2, sum)
      # drop=F prevents R from converting single column to vector. But here is not a single column?
      # rows from count.gene df are selected, and we sum counts for excluded (alternative) junctions for each column (sample) - 2 (apply to column), sum (function to apply)
    }
  
    # inclusion ratio - the proportion of transcripts that contains a segment
    jxn.count.incl.ratio.mtrx = jxns.incl.count.mtrx/jxns.incl.and.excl.count.mtrx # matrix
    jxn.count.incl.ratio.mtrx[jxns.incl.and.excl.count.mtrx < min.cov] = NA # frequency of inclusion is not defined (min.cov = 10, set based on binomial distribution dispersion)
  
    # exclusion counts. Substract included junctions from all
    jxn.count.excl.matr = jxns.incl.and.excl.count.mtrx - jxns.incl.count.mtrx
  
    sajr.gene = list(seg=jxns.incl.and.excl.count.mtrx, i=jxns.incl.count.mtrx, 
                     e=jxn.count.excl.matr, ir=jxn.count.incl.ratio.mtrx)
    class(sajr.gene)='sajr'
    sajr.gene
}


makeSAJR = function(rse.filtered){
  sajr = NULL
  gene.ids = unique(rse.filtered@rowRanges$gene_id)
  # for each cytoskeleton gene
  for(gene.id in gene.ids){
    # selecting information only for the gene
    rse.gene.filtered = rse.filtered[rse.filtered@rowRanges$gene_id==gene.id,]
    # making sajr object for the geme
    sajr.gene =  makeSAJRgene(rse.gene.filtered)
    
    if(is.null(sajr.gene))
      next
    
    # making combined sajr object for all genes
    if(is.null(sajr))
      sajr = sajr.gene
    else{
      for(element in names(sajr)){ # combining seg with seg, e with e, etc for each element, for every gene
        sajr[[element]] = rbind(sajr[[element]],sajr.gene[[element]])}
    }
  }
  sajr
}

calculateMetrics = function(sajr, fetus.samples.ids){
  fetus.indices = match(fetus.samples.ids, colnames(sajr$ir))
  dPSI = apply(sajr$ir, 1, function(x) mean(x[fetus.indices],na.rm=T)-mean(x[-fetus.indices],na.rm=T))
  #logFC = apply(sajr$ir, 1, function(x) log2(mean(x[adult.samples.ids],na.rm=T)/mean(x[fetus.samples.ids],na.rm=T)))
  dPSI #, logFC=logFC)
}

runSAJR = function(rse, tissue){
  rse.filtered = filterRse(rse, tissue)
  sajr.tissue = makeSAJR(rse.filtered)
  mod = rse.filtered@colData$age_group
  mod = list(age_group=factor(mod)) # ~ model data
  fetus.sample.ids = findConditionIds(rse.filtered, 'fetus')

  alt.glm = as.data.frame( fitSAGLM(sajr.tissue, terms(x ~ age_group, keep.order=T),mod,return.pv=T) )
  alt.glm$dPSI = calculateMetrics(sajr.tissue, fetus.sample.ids)
  alt.glm$FDR.sajr = p.adjust(alt.glm$age_group,m='BH')
  names(alt.glm)[names(alt.glm) == "age_group"] = "p.value.sajr"
  
  alt.glm
}

#======================

findSignJxns = function(rse, tissue){
  diego.output = makeDiegoInputFiles(rse, tissue)
  dje.output = runDJExpress(rse, tissue, reference_condition='adult', FDR_threshold=0, logFC_threashold=0)
  sajr.output = runSAJR(rse, tissue)
  list(diego.output=diego.output, dje.output=dje.output, sajr.output=sajr.output)
}

