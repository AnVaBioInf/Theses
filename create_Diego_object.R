#source("data.R")
#source("samples_ann_preprocessing.R")

# selecting only adult brain and testies samples

makeAFile = function(rse.filtered,tissue){
  # -- a file (table of splice junction supports per sample)
  counts = as.matrix(rse.filtered@assays@data$counts)
  counts = as.data.frame(counts)
  
  counts$chr_start_stop = 
    sapply(rownames(counts), function(x) sub('([0-9]+)(:)([-+*])', '\\1' , x))
  counts$type = 'N_w'
  counts$gene_identifier = rse.filtered@rowRanges$gene_id
  counts$gene_name = rse.filtered@rowRanges$gene_names

  a.input.file = counts[c('chr_start_stop', 'type', colnames(counts), 'gene_identifier', 'gene_name')]
  colnames(a.input.file) =  c('junction', 'type', colnames(counts), 'geneID', 'geneName')
  
  rownames(a.input.file) = NULL
  write.table(a.input.file, 
              paste0("/home/an/DIEGO_input_files/a.input.file",tissue), 
              sep = "\t", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}

makeBFile = function(rse.filtered, tissue){
  # --b file (condition to sample relation in the format: condition tab-delimiter sampleName)
  b_file = rse.filtered@colData[,'age_group',drop=F]
  b_file$count_column_name = rownames(b_file)

  # b_file = b_file[order(b_file$tissue), ]
  colnames(b_file) = c('group identifier', 'count_column_name')
  rownames(b_file) = NULL

  write.table(b_file,
              paste0("/home/an/DIEGO_input_files/b_file",tissue),
              sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

makeDiegoInputFiles = function(rse.filtered, tissue){
  makeAFile(rse.filtered, tissue)
  makeBFile(rse.filtered, tissue)
}

# -x file (specify base condition)
# Brain

# Есть ли какая-то фильтрация в DIEGO?
# to run DIEGO in terminal
# conda create -n DIEGO_1 numpy=1.9 scipy matplotlib
# Packages were reinstalled because -e (for drawing dendrograms) wasn't working (numpy 1.9 installed instead, and than other packages reinstalled, in the above code I specified numpy version, idk if it will help)
# conda activate DIEGO_1
# cp /home/an/DIEGO_input_files/a.input.file /home/an/DIEGO_input_files/b_file /home/an/anaconda3/envs/DIEGO_1
# wget http://legacy.bioinf.uni-leipzig.de/Software/DIEGO/DIEGO.tar.gz
# tar -xzf DIEGO.tar.gz
# python DIEGO/diego.py -a a.input.file -b b_file -x fetus --minsupp 1 -d 1 -q 1.0 -z 1.0  > DIEGO_output




# tissues=("Brain" "Heart" "Cerebellum" "Kidney" "Liver" "Testis")
# 
# for tissue in "${tissues[@]}"; do
# python DIEGO/diego.py -a "/home/an/DIEGO_input_files/a.input.file${tissue}" -b "/home/an/DIEGO_input_files/b_file${tissue}" -x fetus --minsupp 1 -d 1 -q 1.0 -z 1.0 > "/home/an/DIEGO_output_files/DIEGO_output${tissue}"
# done