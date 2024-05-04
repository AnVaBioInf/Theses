library(rtracklayer)

jxn = readRDS('rse.ERP109002.jxn.rds', refhook = NULL)
genes = readRDS('rse.ERP109002.rds', refhook = NULL)
genes.annotaion = readRDS('genes.annotaion.rds', refhook = NULL)

genes.annotaion


#--downloading bigwig files
# 
# for (sample in rownames(jxn@colData)[143:length(rownames(jxn@colData))]){
#   bigwig <- jxn[, sample]$BigWigURL
#   destfile <- paste0('/home/an/Manananggal/Input/bigWig/', sample, ".bw")
# 
#   max_retries <- 10  # Set maximum retry attempts
#   num_retries <- 0
# 
#   while (num_retries < max_retries) {
#     result <- tryCatch({
#       download.file(bigwig, destfile = destfile, timeout = 600, mode = "wb")
#       TRUE  # Download successful
#     }, error = function(e) {
#       FALSE  # Download failed
#     })
# 
#     if (result) {
#       message("Download successful for sample:", sample)
#       break  # Exit loop if successful
#     } else {
#       num_retries <- num_retries + 1
#       message("Download failed for sample:", sample, "- Retry attempt:", num_retries)
#     }
#   }
# 
#   if (num_retries == max_retries) {
#     message("Download failed for sample:", sample, "after maximum retries.")
#     # Implement additional error handling or notification
#   }
# }

# # chromosome size file
# #BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# library(BSgenome.Hsapiens.UCSC.hg38)
# 
# chrom_info <- seqlengths(Hsapiens)
# write.table(data.frame(chrom = names(chrom_info), size = chrom_info),
#             file = "chrom.sizes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#-- bedgraph files
gene.grange = genes@rowRanges[genes.annotaion$gene_id]
file_paths = list.files("/home/an/Manananggal/Input/bigWig", full.names = TRUE)
# i=0
# for (file in file_paths){
#   i=i+1
#   print(i)
#   print(file)
#   
#   sample.name = sub("\\.bw$", "", basename(file))
#   path = '/home/an/Manananggal/Input/begGraph/'
#   
#   print(sample.name)
#   
#   bw = rtracklayer::import.bw(file,which=gene.grange)
#   export(bw, paste0(path,sample.name,".bedGraph"), format = "bedGraph")
#   
#   print(bw)
# }

# checking
table(rownames(jxn@colData) %in% unname(sapply(file_paths, function(x) sub("\\.bw$", "", basename(x)))))
file_paths[!unname(sapply(file_paths, function(x) sub("\\.bw$", "", basename(x)))) %in%rownames(jxn@colData)]

# medgraph back to bigwig (in terminal)
# for file in *.bedGraph; do   '/home/an/bedGraphToBigWig'  "$file" '/home/an/chrom.sizes.txt'  "/home/an/Manananggal/Input/bigWigFiltered/${file%.bedGraph}.bw"; done


## -- count files
source('data.R')
counts = as.matrix(rse.ERP109002.jxn.cytosk.genes@assays@data$counts)
counts =  as.data.frame(counts)
counts$type = "Junction"

gene_ids = rse.ERP109002.jxn.cytosk.genes@rowRanges$gene_id

counts$gene_id = sub("\\..+$", '', gene_ids)
counts$strand = sub("(.)+(.$)",'\\2',  rownames(counts))
counts$forward_strand = ifelse(counts$strand == "+", "true",
                               ifelse(counts$strand == "-", "false", "unknown"))


chr_start_end = strsplit(rownames(counts), "[:-]")
counts$position = paste0(sapply(chr_start_end, "[[", 2),"-", sapply(chr_start_end, "[[", 3))
counts$reference <- sapply(chr_start_end, "[[", 1)
counts$start <- sapply(chr_start_end, "[[", 2)

counts$end <- sapply(chr_start_end, "[[", 3)
counts$place_holder = 0
print(counts)
print(counts[,c('strand', 'forward_strand')])
## -- project files

unique(counts$gene_id)

source('samples_ann_preprocessing.R')
# age+tissue 
ann = rse.ERP109002.jxn.cytosk.genes@colData
length(rownames(ann))



for (tissue in unique(ann$tissue)){
  sample.id = rownames(ann[ann$tissue==tissue,])
  print(sample.id)
  c = counts[rowSums(counts[,sample.id]) > 0,][sample.id]
  print(c)
  print(dim(c))
}
  
length('ENSG00000158195')

length(unique(counts$reference))

for (sample in sample.id){
  input.file = counts[, c('type', 'gene_id', 'forward_strand', 'position', 'reference', 'start', 'end', 'place_holder',sample, 'place_holder')]
  #colnames(input.file) =  c('reference', 'start', 'end', 'gene_id', 'count', 'strand')
  
  print(input.file)
  #  input.file <- input.file[input.file[,sample] != 0, ]
  print(input.file)
  
  
  write.table(input.file,
              paste0("/home/an/Manananggal/Input/counts/",sample, ".txt"),
              sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}



## -- project files
source('samples_ann_preprocessing.R')
# age+tissue 
ann = rse.ERP109002.jxn.cytosk.genes@colData

for (tissue in unique(ann$tissue)){
  project.file = ann[ann$tissue==tissue, c('age_group', 'tissue')]
  project.file$sample = rownames(project.file)
  project.file$size_factors = 0
  project.file$bigwig_files = paste0('/home/an/Manananggal/Input/bigWigFiltered/',
                                     rownames(project.file), '.bw')
  project.file$junction_count_files = paste0('/home/an/Manananggal/Input/counts/',
                                    rownames(project.file), '.txt')
  project.file = project.file[,c('sample',
                        'age_group',
                        'size_factors',
                        'bigwig_files',
                        'junction_count_files')]
  colnames(project.file) = c('sample',
                             'condition',
                             'size_factors',
                             'bigwig_files',
                             'junction_count_files')
  write.table(project.file,
              paste0("/home/an/Manananggal/Input/project_files/", tissue, ".project"),
              sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
} 



genes_to_keep = unique(rse.ERP109002.jxn.cytosk.genes@rowRanges$gene_id)






#-- filtering annotation
# Replace with the path to your GENCODE annotation file
file_path <- "/home/an/Manananggal/Input/ref_annotation/gencode.v26.annotation.gtf"

# Load required package
library(rtracklayer)


# Function to parse GENCODE annotation and select genes
parse_gencode <- function(file_path, genes_to_keep) {
  
  # Read the annotation file (assuming GTF format)
  annotation <- import(file_path)
  
  # Select genes based on gene_id
  selected_annotation <- annotation[annotation$gene_id %in% genes_to_keep,]
  
  # Return the selected annotation
  return(selected_annotation)
}


selected_annotation = parse_gencode(file_path, genes_to_keep)
export(selected_annotation, "/home/an/Manananggal/Input/ref_annotation/filtered_ann_v26.gtf")






#-- cross reference file
library(biomaRt)

ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")

# valid_attributes <- listAttributes(mart = ensembl)
# print(valid_attributes)

bm <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "refseq_mrna", "entrezgene_id", "hgnc_symbol", "description"), mart = ensembl)
bm=bm[bm$ensembl_gene_id %in% sub("\\..+$", '', genes_to_keep),]

names(genes_to_keep) = sub("\\..+$", '', genes_to_keep)
genes_to_keep

match(bm$ensembl_gene_id, names(genes_to_keep))

bm$ensembl_gene_id = unname(genes_to_keep[match(bm$ensembl_gene_id, names(genes_to_keep))])


# 
     
# write as tab-separated file

write.table(bm, "/home/an/Manananggal/Input/biomart/biomart_cross_ref_human.txt", row.names=FALSE, quote=FALSE, sep="\t")




java -cp '/home/an/Manananggal/WebContent/WEB-INF/lib/*':'/home/an/Manananggal/jar_files/Manananggal_v1.0.7.jar'  Manananggal.SplicingAnalyzer calculate_size_factors '/home/an/Manananggal/Input/ref_annotation/filtered_ann_v26.gtf'   '/home/an/Manananggal/Input/project_files/Brain.project' 

java -cp '/home/an/Manananggal/WebContent/WEB-INF/lib/*':'/home/an/Manananggal/jar_files/Manananggal_v1.0.7.jar' merge '/home/an/Manananggal/Input/project_files/Brain.project' 0 0 0 0




java -cp './WebContent/WEB-INF/lib/*':'./jar_files/Manananggal_v1.0.7.jar' Manananggal.SplicingAnalyzer './Input/ref_annotation/filtered_ann_v26.gtf' './Input/project_files/Brain.project' condition 3 5 0.7 0.05 './Input/biomart/biomart_cross_ref_human.txt'  -skipFirstAndLast -skipIntronDetection -threads=6
sudo su - tomcat -c "/opt/apache-tomcat-7.0.99/bin/startup.sh"
sdk use java 8.0.392.fx-zulu 

http://localhost:8080/Manananggal_v1.0.7/
  
sudo docker build -t manananggal .
sudo docker run -v ./Input/web:/data -p 8080:8080 manananggal

