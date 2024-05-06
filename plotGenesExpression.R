library("recount3")
library("dbplyr")
library("dplyr")
library("tidyr")
library("readr")
library('tibble')
library("gtools")
library('gridExtra')
library("tidyr")
library("ggplot2")


#----------------downloading dataset-----------------------------
human_projects <- recount3::available_projects()

proj_info <- subset(
  human_projects,
  project == "ERP109002" & project_type == "data_sources"
)

# access the gene level expression data
# create RSE object 
rse_genes_ERP109002 <- create_rse(proj_info)

#-------------------------sample metadata-----------------------------
# from sample metadata extracting sra.library_name, containing information on tissue, age etc
samples_meta_orig = colData(rse_genes_ERP109002)['sra.library_name'] 
head(samples_meta_orig)
# ERR2598051 5664sTS.Human.Forebr..

# --- data cleaning
# separating data to different columns
samples_meta_orig = as.data.frame(samples_meta_orig) %>%
  separate(sra.library_name, c("sample_id", "organism", "tissue", "age", "gender"), "\\.")   # \\. - delimiter dot
head(samples_meta_orig)
# ERR2598051   5664sTS    Human Forebrain 10w Female_s

# combining some tissues
samples_meta = as.data.frame(samples_meta_orig)
samples_meta$tissue = replace(samples_meta$tissue, samples_meta$tissue == 'Forebrain', "Brain")
samples_meta$tissue = replace(samples_meta$tissue, samples_meta$tissue == 'Hindbrain', "Cerebellum")
samples_meta$tissue = replace(samples_meta$tissue, samples_meta$tissue == 'KidneyTestis', 'Kidney')
samples_meta$tissue %>% unique()
# [1] "Brain"      "Heart"      "Cerebellum" "Kidney"     "Liver"      "Ovary"      "Testis"    

# CS - Carnegie Stage to weeks
cs_to_days = c('CS13' = '32',
               'CS14' = '33',
               'CS15' = '36',
               'CS16' = '39',
               'CS17' = '41',
               'CS18' = '44',
               'CS19' = '46',
               'CS20' = '49',
               'CS21' = '51',
               'CS22' = '53',
               'CS23' = '56') # dpc
# transferring it to weeks
cs_to_w = as.numeric(cs_to_days)/4
cs_to_w = paste(cs_to_w, "w", sep="") # adding w, that's how all values are stored in age column age+measure
names(cs_to_w) = names(cs_to_days) # asigning CS stages in correspondance to each week
cs_to_w
# CS13     CS14     CS15     CS16     CS17     CS18     CS19     CS20     CS21     CS22     CS23 
# "8w"  "8.25w"     "9w"  "9.75w" "10.25w"    "11w"  "11.5w" "12.25w" "12.75w" "13.25w"    "14w" 


# replacing CS with wpc
cs_idx = which(samples_meta$age %in% names(cs_to_w)) # index of cs values to replace
cs_values = samples_meta$age[cs_idx] # selecting all cs values
samples_meta$age = replace(samples_meta$age, cs_idx, cs_to_w[cs_values]) #  
samples_meta$age[cs_idx]
#  [1] "8 wpc"     "8.25 wpc"  "8 wpc"     "8 wpc"     "8 wpc"     "9.75 wpc"  "9.75 wpc"  "9.75 wpc"  "11 wpc"    "10.25 wpc"
# ...

# separating age and measure into two columns
# replace w (weeks) with wpc (weeks post conseption)
samples_meta$age = gsub('w$', ' wpc', samples_meta$age)
samples_meta$age %>% unique()
#  [1] "10 wpc"    "11 wpc"    "12 wpc"    "13 wpc"    "16 wpc"    "18 wpc"    "19 wpc"    "20 wpc"    "8 wpc"     "4 wpc"    
# ...

# -- splitting data into age groups
ages = samples_meta$age %>% unique()
ages
#  [1] "10 wpc"    "11 wpc"    "12 wpc"    "13 wpc"    "16 wpc"    "18 wpc"    "19 wpc"    "20 wpc"    "8 wpc"     "4 wpc"    
# ...

# order of columns (age)
ages_order = c()
for (i in c('dpc', 'wpc', 'dpb', 'mpb', 'ypb')){
  ages_order = c(ages_order, mixedsort(ages[endsWith(ages, i)])) # mixedsort function sorts values that are a mixture of numb and str
} # appending to vector
ages_order # in the next step age groups are formed manualy, using the output of this function
#  [1] "4 wpc"     "5 wpc"     "6 wpc"     "7 wpc"     "8 wpc"     "8.25 wpc"  "9 wpc"     "9.75 wpc"  "10 wpc"    "10.25 wpc"
# ...

# merging age groups
infant = c("0dpb","4dpb","6dpb","15dpb" ,"18dpb","19dpb" ,"34dpb","94dpb" ,"127dpb","221dpb","226dpb","270dpb","271dpb","6mpb","1ypb")
toddler = c('2ypb', '4ypb')
school = c("7ypb","8ypb"   )
teen = c("13ypb","14ypb","16ypb","17ypb" )
yo_25_35 = c("25ypb","28ypb","29ypb","32ypb" )
yo_36_45 = c("39ypb" )
yo_46_55 = c( "46ypb", "50ypb", "53ypb" , "54ypb", "55ypb")
yo_56_65 = c("58ypb")
all_age_gr = c(infant, toddler, school, teen, yo_25_35, yo_36_45, yo_46_55, yo_56_65)

samples_meta$age[samples_meta$age %in% infant] = "infant"
samples_meta$age[samples_meta$age %in% toddler] = "toddler"
samples_meta$age[samples_meta$age %in% school] = "school"
samples_meta$age[samples_meta$age %in% teen] = "teen"
samples_meta$age[samples_meta$age %in% yo_25_35] = "25-35 years"
samples_meta$age[samples_meta$age %in% yo_36_45] = "36-45 years"
samples_meta$age[samples_meta$age %in% yo_46_55] = "46-55 years"
samples_meta$age[samples_meta$age %in% yo_56_65] = "56-65 years"
head(samples_meta)
# ERR2598051   5664sTS    Human      Brain    10 wpc Female_s
# ...

order = c(ages_order[endsWith(ages_order, "wpc")], 
          c("infant", "toddler", "school", "teen", "25-35 years", "36-45 years", "46-55 years", "56-65 years")
          ) # setting new order
order

# checking the introduced changes
dim(samples_meta_orig)
dim(samples_meta)
head(samples_meta_orig)
head(samples_meta)

# -----------------------------------------------------------------------------------
# ------------------------- samples occurrence heatmap ------------------------------
# -----------------------------------------------------------------------------------
# occurrence of samples
sample_occurance = table(samples_meta$tissue, samples_meta$age) %>% as.data.frame.matrix() # as.data.frame
# making a column for tissues
sample_occurance$tissue = rownames(sample_occurance)
sample_occurance

# reshaping the df to plot heatmap
sample_occurance_long <- pivot_longer(
                            data = sample_occurance, 
                            cols = -length(colnames(sample_occurance)), # all, aside from the last one (tissue)
                            names_to = "age_group", 
                            values_to = "sampl_numb") # increases the number or rows and decreases the number of columns
sample_occurance_long

sample_occurance_heatmap <- ggplot(data = sample_occurance_long, mapping = aes(x = tissue,
                                                       y = age_group,
                                                       fill = sampl_numb)) + # fill - tile color depending on # of samples
  geom_tile() + # makes a 'tile' for each tissue+age
  scale_fill_gradient(low = "#FFFFFF",
                      high = "#FF2400") + # setting color palette
  geom_text(aes(label = sampl_numb)) + # labels tiles with numbers
  scale_y_discrete(limits = order) + # setting order of y values
  labs(title = 'age~tissue, number of samples',
       x = 'Tissue',
       y = 'Age',
       fill='# of samples')
sample_occurance_heatmap

# ----------------------------------------------------------------------------------
# ------------------------ gene expression vs time graphs -------------------------
# ----------------------------------------------------------------------------------
# expressions
quant = assays(rse_genes_ERP109002)[[1]] %>%
  as.matrix() %>%
  as.data.frame()
head(quant)

# normalisation (CPM)
# raw reads mapped to the transcript / sum counts per sample *10**6
quant_cpm = sweep(x = quant, MARGIN = 2, STATS = colSums(quant), FUN = `/`)*10**6
# MARGIN = 1 means row; MARGIN = 2 means column.
# STATS - the value(s) that should be added or subtracted
# FUN The operation that has to be done
head(quant_cpm)

# --- data preparation
# - gene annotation
# citoskeleton genes
genes = rbind(data.frame(names=c('CYFIP1','CYFIP2','NCKAP1','NCKAP1L','ABI1','ABI2','ABI3','WASF1','WASF2','WASF3','BRK1'),group='WAVE'),
              data.frame(names=c('NHS','NHSL1','NHSL2','KIAA1522'),group='NHS'),
              data.frame(names=c('ARPC1A','ARPC1B','ARPC2','ARPC3','ARPC4','ARPC5','ACTR2','ACTR3','ACTR3B'),group='Arp2/3'))
genes

# annotation for citoskeleton genes
rse_genes_ERP109002@rowRanges
#          ENSG00000278704.1 GL000009.2       56140-58376      - |  ENSEMBL     gene      2237      <NA>

annotation_citosk_genes = 
  rse_genes_ERP109002@rowRanges[rse_genes_ERP109002@rowRanges$gene_name %in% genes$names,] %>%
  as.data.frame() # from gene annotation leaving only citoskeleton genes
annotation_citosk_genes$gene_name %>% unique()
#  [1] "CYFIP1"   "WASF2"    "KIAA1522" "ARPC5"    "ABI1"     "NCKAP1L"  "ARPC3"    "WASF3"    "ABI3"    
# ...

# adding information about group = molecular complex
# matching genes and extracting their groups
annotation_citosk_genes$group = genes$group[match(annotation_citosk_genes$gene_name,genes$names)]
annotation_citosk_genes$group
# it matches first vector with second, order-1st, ids-2nd

annotation_citosk_genes$gene_name
genes$names

genes[match(annotation_citosk_genes$gene_name,genes$names),]

setdiff(genes$names, annotation_citosk_genes$gene_name) # checking if all genes were selected
# character(0)
dim(genes)
# [1] 24  2  
# 24 genes, 2 columns
dim(annotation_citosk_genes)
# [1] 25 15
# 1 gene more!
sort(table(annotation_citosk_genes$gene_name))
annotation_citosk_genes[annotation_citosk_genes$gene_name=='CYFIP1',]
# one of the genes belongs to scaffolds

# filtering out scaffolds from annotation
annotation_citosk_genes[!startsWith(as.character(annotation_citosk_genes$seqnames),'chr'),] #  scaffolds
# ENSG00000280618.2 protein_coding    CYFIP1     2 OTTHUMG00000189940.2 <NA>  WAVE
annotation_citosk_genes = annotation_citosk_genes[startsWith(as.character(annotation_citosk_genes$seqnames),'chr'),] # selecting only chromosomes, excluding scaffolds
annotation_citosk_genes # genes belonging to chr

# leaving only counts for citoskeleton genes on chr
quant_cpm = quant_cpm[rownames(annotation_citosk_genes),] # leaving only counts for chr genes (not scaffold genes)
# adding gene name
quant_cpm$gene_name = annotation_citosk_genes$gene_name 
quant_cpm$group = annotation_citosk_genes$group 
head(quant_cpm)

#melt data frame into long format
quant_cpm = reshape2::melt(
  quant_cpm,
  variable.name = 'sample_id',
  value.name = 'counts',
  variable.factor=FALSE)
# melt converts df from wide to long format 
# in long format values repeat in the first column
head(quant_cpm)


# adding info about samples from annotation
# sample tissues and age
quant_cpm$tissue = samples_meta[match(quant_cpm$sample_id, rownames(samples_meta)),]$tissue
quant_cpm$age = samples_meta[match(quant_cpm$sample_id, rownames(samples_meta)),]$age
quant_cpm[, 'age'] <- as.factor(quant_cpm[, 'age'])
head(quant_cpm)

# order - order of age, matching it with a number,
# this is needed to walk around the fact that loess needs numbers on
# both axis
names(order) = c(1:length(order))
quant_cpm$numb = match(quant_cpm$age, order) # adding a column that will indicate order of time
# output - idx of 1st vector, corresponding to second
head(quant_cpm)


# --- plotting
tissue.col=c(Brain="#3399CC",
             Cerebellum="#33CCFF",
             Heart="#CC0000",
             Kidney="#CC9900",
             Liver="#339900",
             Ovary="#CC3399",
             Testis="#FF6600")

tissues = unique(quant_cpm$tissue)

for (i in genes$group %>% unique()){
  gene_counts = quant_cpm[quant_cpm$group == i,] # selecting a genes belonging to the same group
  print(ggplot(data=gene_counts, 
       aes(x=as.numeric(numb), 
           y = counts, 
           colour = factor(tissue), 
           group = tissue)) +
  geom_point() +
  stat_smooth(method="loess", se=F) + # linear regression curve
  scale_color_manual(values = tissue.col) + # setting color pannel
  labs(title = i,
       x = "Age",
       y = "CPM") +
  scale_x_continuous(breaks = c(1:30),
                     labels = unname(order))+  # replacing numbers with age
  scale_y_continuous(limits = c(0,500)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)  )+ # rotating axes names
  facet_wrap(~gene_name)) # shared axes 
}
