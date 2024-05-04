library("recount3")
library("dbplyr")
library("dplyr")
library("tidyr")
library("readr")
library('tibble')


#----------------downloading dataset-----------------------------
human_projects <- recount3::available_projects()

proj_info <- subset(
  human_projects,
  project == "ERP109002" & project_type == "data_sources"
)

# access the gene level expression data
# create RSE object 
rse_jxn_ERP109002 <- create_rse(proj_info,
                                type="jxn")
metadata(rse_jxn_ERP109002)
#  For human, we used the UCSC hg38 assembly, based on GRCh38.
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02533-6#Sec19

# --------------converting raw counts to CPM--------------------
assayNames(rse_jxn_ERP109002)
# ??? counts! not raw_counts. 

## Scale the counts using the AUC
#assays(rse_jxn_ERP109002)$counts <- transform_counts(rse_jxn_ERP109002)


#-------------------------sample metadata-----------------------------
# https://bioconductor.org/packages/3.18/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# finding only heart and testicles samples
samples_meta = colData(rse_jxn_ERP109002)['sra.library_name'] 
samples_meta

# data cleaning
# separating data to different collumns
samples_meta = as.data.frame(samples_meta) %>%
  separate(sra.library_name, c("sample_id", "organism", "tissue", "age", "gender"), "\\.")
head(samples_meta)

# unifying time units to years

# CS - Carnegie Stage
# data[grepl("^CS*", data$age),]['age'] %>% unique() 

# из вики взяты соответствия CS и дней
cs_to_days = c('CS13' = '32dpc',
               'CS14' = '33dpc',
               'CS15' = '36dpc',
               'CS16' = '39dpc',
               'CS17' = '41dpc',
               'CS18' = '44dpc',
               'CS19' = '46dpc',
               'CS20' = '49dpc',
               'CS21' = '51dpc',
               'CS22' = '53dpc',
               'CS23' = '56dpc')
names(cs_to_days)
# 4-20 wpc
# as.numeric(gsub("([0-9]+).*$", "\\1", cs_to_days)) / 4
#  8.00  8.25  9.00  9.75 10.25 11.00 11.50 12.25 12.75 13.25 14.00 wpc

# age = samples_meta$age %>% unique()
# age[grep("w$", age)] %>% sort()
# "10w" "11w" "12w" "13w" "16w" "18w" "19w" "20w" "4w"  "5w"  "6w"  "7w"  "8w"  "9w" 

# replacing CS with dpc
idx = which(samples_meta$age %in% names(cs_to_days))
samples_meta$age = replace(samples_meta$age, idx, cs_to_days[samples_meta$age[idx]])

# separating age and measure into two columns
samples_meta = samples_meta %>%
  separate(age, into=c("age", "measure"), sep="(?=[a-z +]+)(?<=[0-9])")

# setting column age as numeric
samples_meta = transform(samples_meta, age = as.numeric(age))
head(samples_meta)

# checking column data types
sapply(samples_meta, mode)


# time units used in the study
# data_copy$measure %>% unique()
# "w"   "dpc" "ypb" "dpb" "mpb"

# converting everything to days
# function that converts days/weeks/months into days from conception 
convert_to_ypc = function(x) {
  value = as.numeric(x[1])
  unit_of_time = x[2]
  if (unit_of_time == 'dpc') {
    value / 365.25
  } 
  else if (unit_of_time == 'dpb') {
    ( value + 280 ) / 365.25 # 280 days is the duration of pregnancy 
  }
  else if (unit_of_time == 'w') {
    ( value * 7 ) / 365.25 # weeks pc
  }
  else if (unit_of_time == 'mpb') {
    ( value * 30*4 + 280 ) / 365.25  # 30.4 - mean number of days in a month
  }
  else if (unit_of_time == 'ypb') {
    value  + 280 / 365.25
  }
}

samples_meta$age = apply(samples_meta[c('age', 'measure')], 1, convert_to_ypc)
samples_meta$measure = rep("ypc", length(samples_meta$measure))
samples_meta$sample_ids = rownames(samples_meta)
head(samples_meta)

# ---  sample-age-organ
#id_age = samples_meta$age
#names(id_age) = rownames(samples_meta)
#id_age

samples_meta = samples_meta[c('age', 'tissue')]
samples_meta

# --- organ-samples
# list if organs presented in the dataset
#organs = samples_meta$tissue %>% unique()
#organs

# list of ids for each organ
#organ_ids = list()
#for (i in organs){
#  organ_ids[i] = list(samples_meta[samples_meta$tissue==i,] %>% rownames() )
#}
#organ_ids

# ---------------------------- expression --------------------------------



#------------------------example Forebrain
organ = 'Forebrain'

# expressions
organ_samples_ids= samples_meta[samples_meta['tissue']==organ,] %>% rownames()

quant = assays(rse_jxn_ERP109002)[[1]][,organ_samples_ids] %>%
  as.matrix() %>%
  as.data.frame()
head(quant)


# ---------------------- gene annotation -----------------------------
jxt = rowRanges(rse_jxn_ERP109002)
head(jxt)

# creating granges object from genomic annotation Gencode v26
library(rtracklayer)
anno_gencode26 <- import('/home/an/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf')
head(anno_gencode26)

# subsetting gene (!) records (filtering out transcript etc records)
anno_gencode26_genes = subset(anno_gencode26,type %in% 'gene')[,'gene_name']
anno_gencode26_genes

# findOverlaps
# finding exons overlapping with genes from Uniprot v.26 annotation
overlaps = findOverlaps(jxt, 
             anno_gencode26_genes, 
             type='any')

# selecting only overlapped exons
quant_annotated = quant[queryHits(overlaps),]
quant_annotated

# assigning gene names to exons
quant_annotated$gene_name = anno_gencode26_genes[subjectHits(overlaps)]$gene_name

# leaving only cytosceleton genes
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

# length(genes_citosceleton)
# 24 citosceleton genes

quant_annotated_citosk = quant_annotated[which(
  quant_annotated$gene_name %in% genes_citosceleton),]
quant_annotated_citosk

# group by gene and sum expr
quant_annotated_citosk = group_by(quant_annotated_citosk,
         gene_name) %>%
  summarise_all(sum) %>%
  as.data.frame()
quant_annotated_citosk

# ----------------------------- plotting ----------------------------------
# tissue
library(ggplot2)
library(reshape2)

# replace ids with corresponding age

############################
samples_meta[samples_meta$tissue==organ,]


colnames(quant_annotated_citosk) = c("gene_name", samples_meta[colnames(quant_annotated_citosk)[-1],][['age']] %>%
                                       round(digits = 3))

quant_annotated_citosk
#melt data frame into long format
quant_annotated_citosk_melt = melt(
  quant_annotated_citosk,
  id = 'gene_name',     
  variable.name = 'age',
  value.name = 'counts',
  variable.factor=FALSE)
quant_annotated_citosk_melt$age = as.numeric(as.character(quant_annotated_citosk_melt$age))
head(quant_annotated_citosk_melt)
#class(quant_annotated_citosk_melt)
#sapply(quant_annotated_citosk_melt, mode)

#create line plot for each column in data frame
ggplot(quant_annotated_citosk_melt, aes(age, counts, colour = gene_name)) + 
  geom_line()+ 
  ggtitle(organ) +
  xlab("age, ypc") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
