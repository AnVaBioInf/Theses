library("recount3")
library("dbplyr")
library("dplyr")
library("tidyr")
library("readr")
library('tibble')
library("gtools")
library('gridExtra')
library("tidyr")
# рфзобраться, что из этого действительно нужно

source("data.R")

# import file data.R

# table with metadata
# from sample metadata extracting sra.library_name, containing information on tissue, age etc
samples_meta_orig = colData(rse.ERP109002.jxn.cytosk.genes)['sra.library_name'] 
head(samples_meta_orig)
# ERR2598051 5664sTS.Human.Forebr..

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
samples_meta

# -- adding column with age group
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
fetus = ages_order[grep("pc", ages_order)]
infant = c("0dpb","4dpb","6dpb","15dpb" ,"18dpb","19dpb" ,"34dpb","94dpb" ,"127dpb","221dpb","226dpb","270dpb","271dpb","6mpb","1ypb")
toddler = c('2ypb', '4ypb')
school = c("7ypb","8ypb"   )
teen = c("13ypb","14ypb","16ypb","17ypb" )
yo_25_35 = c("25ypb","28ypb","29ypb","32ypb" )
yo_36_45 = c("39ypb" )
yo_46_55 = c( "46ypb", "50ypb", "53ypb" , "54ypb", "55ypb")
yo_56_65 = c("58ypb")
all_age_gr = c(infant, toddler, school, teen, yo_25_35, yo_36_45, yo_46_55, yo_56_65)

samples_meta$age_group[samples_meta$age %in% fetus] = "fetus"
# samples_meta$age_group[samples_meta$age %in% infant] = "infant"
# samples_meta$age_group[samples_meta$age %in% toddler] = "toddler"
# samples_meta$age_group[samples_meta$age %in% school] = "school"
# samples_meta$age_group[samples_meta$age %in% teen] = "teen"
# samples_meta$age_group[samples_meta$age %in% c(yo_25_35, yo_36_45, yo_46_55, yo_56_65)] = "adult"
# 
# head(samples_meta)
# ERR2598051   5664sTS    Human      Brain 10 wpc Female_s before_birth

# проверить!

samples_meta$age_group[
  samples_meta$tissue != 'Testis' &
  samples_meta$age %in% c(toddler, school, teen, yo_25_35, yo_36_45, yo_46_55, yo_56_65)] = "adult"

samples_meta$age_group[
  samples_meta$tissue == 'Testis' &
    samples_meta$age %in% c(yo_25_35, yo_36_45, yo_46_55, yo_56_65)] = "adult"

samples_meta = samples_meta[samples_meta$age_group %in% c('fetus','adult'), ]

table(samples_meta$age_group)
table(samples_meta$tissue)

# ...
# 
# #  -----------------
# brain_adult_samples = samples_meta[samples_meta$tissue=="Brain" & samples_meta$age_group=="adult",]
# testies_adult_samples = samples_meta[samples_meta$tissue=="Testis" & samples_meta$age_group=="adult",]
# 
# 
# 
# #----------------??
# b = rownames(brain_adult_samples)
# t = rownames(testies_adult_samples)
# b
# t
