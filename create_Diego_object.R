#install.packages('devtools')
#devtools::install_version("dbplyr", version = "2.3.4")

#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}
#BiocManager::install(c("recount3", "SRAdb"))
## Check that you have a valid Bioconductor installation
#BiocManager::valid()

#install.packages("tidyr")

library("recount3")
library("dbplyr")
library("tidyr")
library("readr")

#----------------downloading dataset-----------------------------
human_projects <- recount3::available_projects()

proj_info <- subset(
  human_projects,
  project == "ERP109002" & project_type == "data_sources"
)

# access the gene level expression data
# create RSE object 
rse_exon_ERP109002 <- create_rse(proj_info,
                                type="exon")
## !! Выбрала экзоны, потому что DIEGO с ними тоже работает и у них есть 
# генная аннотация, а у junction нет. Наверное, нужно подумать как переделать,
# чтобы входные данные между программами не различались?

#-------------------------gene metadata-------------------------------

# genes metadata matrix
gene_meta = rowRanges(rse_exon_ERP109002)@elementMetadata 
gene_meta@listData$type = as.character(gene_meta@listData$type)
gene_meta = gene_meta %>%
  as.matrix()

# leave only data needed in Diego input file
gene_ids = gene_meta[,c('recount_exon_id', 'type', 'gene_id', 'gene_name')]
head(gene_ids)


#-------------------------sample metadata-----------------------------

# https://bioconductor.org/packages/3.18/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# finding only heart and testicles samples
samples = colData(rse_exon_ERP109002)
sample_ids = samples[grep("Heart.10w.Male|Testis.10w.Male", samples$sra.library_name),] %>% 
  row.names()

#--------------------------counts-----------------------------

# exon expression data (heart, testicles) 
# converting sparse counts matrix to matrix
exon_counts = assays(rse_exon_ERP109002)[[1]][, sample_ids] %>%
  as.matrix()
exon_counts = cbind(recount_exon_id=rownames(exon_counts), exon_counts) # cbind converts all columns to character, since matrix can contain only one type of data
head(exon_counts)

#===================================================================================
#==============================DIEGO=================================================
#====================================================================================
#--------------------------diego imput-----------------------------
# merge counts and metadata matrixes
data = merge(exon_counts, gene_ids, by = "recount_exon_id", all.x = TRUE, all.y=FALSE, sort=FALSE)
head(data)

# -----------------------cytoskeleton genes
# filtering matrixes
# leaving onle cytosceleton genes
genes_citosceleton = c('WASF2', 
                       'KIAA1522', 
                       'ARPC5', 
                       'ABI1', 
                       'NCKAP1L', 
                       'ARPC3', 
                       'WASF3', 
                       'CYFIP1', 
                       'ABI3', 
                       'ACTR2', 
                       'ACTR3', 
                       'NCKAP1', 
                       'ABI2', 
                       'ARPC2',
                       'ARPC4', 
                       'BRK1', 
                       'CYFIP2', 
                       'WASF1', 
                       'NHSL1', 
                       'ARPC1A', 
                       'ARPC1B', 
                       'ACTR3B', 
                       'NHS', 
                       'NHSL2')

diego_count_input = data[data$gene_name %in% genes_citosceleton,]
diego_count_input

# check if all genes of interest are in a dataset
unique(diego_count_input$gene_name) %>%
  length() == length(genes_citosceleton)

# ----------------- filtering low covarage genes
# filtering out exonы with low expression
expressed_genes = diego_count_input[ , grepl( "ERR" , names( diego_count_input ) ) ] %>%
  sapply(as.numeric) %>%
  apply( 1, FUN = min) >2

diego_count_input = diego_count_input[expressed_genes,]


# ------------------formatting data for Diego
# editing matrix to fit Diego requirements for an input file
diego_count_input = diego_count_input %>%
  separate(recount_exon_id, c("chr", "start", "end", "strand"), "\\|")

diego_count_input[,"chr:start-stop"] = 
  paste0(diego_count_input[,'chr'], ":", diego_count_input[,'start'], "-", diego_count_input[,'end'])
head(diego_count_input)
c(c('1', 
  'type'), sample_ids)

colnames=c(
  c('chr:start-stop', 'type'), 
  sample_ids,
  c('gene_id', 'gene_name')
  )

diego_count_input = diego_count_input[, colnames]
colnames(diego_count_input)[ncol(diego_count_input)] = "gene name"
colnames(diego_count_input)[ncol(diego_count_input)-1] = "gene identifier"
head(diego_count_input)

# extracting to file with tab delim
write_delim(diego_count_input, "/home/an/diego_count_input.txt", delim = "\t")


# -------------- Diego input file #2 (group-sample)

condition = samples[grep("Heart.10w.Male|Testis.10w.Male", samples$sra.library_name),][c('sra.library_name')]

sub(".*?\\.", "", condition$sra.library_name) %>% length()
row.names(condition) %>% length()


diego_meta_input = data.frame(sub(".*?\\.", "", condition$sra.library_name),
           row.names(condition))
diego_meta_input

write_delim(diego_meta_input, "/home/an/diego_meta_input.txt", delim = "\t", col_names=FALSE)

