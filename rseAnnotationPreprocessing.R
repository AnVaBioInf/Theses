source("downloadRseData.R")
# metadata
# from sample metadata extracting sra.library_name, containing information on tissue, age etc
split.sra.library_name = strsplit(rse.jxn.cytosk@colData$sra.library_name, "\\.")
# create a matrix to store the extracted values
ncols = lengths(split.sra.library_name)[1]
matrix.sra.library_name = matrix(nrow = nrow(rse.jxn.cytosk@colData), ncol = ncols)
colnames(matrix.sra.library_name) = c("sample_id", "organism", "tissue", "age", "gender")

for (i in 1:nrow(matrix.sra.library_name)) {
  matrix.sra.library_name[i, ] = split.sra.library_name[[i]]
}
# add the extracted columns to the colData
rse.jxn.cytosk@colData <- cbind(rse.jxn.cytosk@colData, matrix.sra.library_name)
rm(matrix.sra.library_name, split.sra.library_name, ncols)

# ---
# combining some tissues
rse.jxn.cytosk@colData$tissue = sub("Forebrain", "Brain", rse.jxn.cytosk@colData$tissue)
rse.jxn.cytosk@colData$tissue = sub("Hindbrain", "Cerebellum", rse.jxn.cytosk@colData$tissue)
rse.jxn.cytosk@colData$tissue = sub("KidneyTestis", "Kidney", rse.jxn.cytosk@colData$tissue)
unique(rse.jxn.cytosk@colData$tissue)
# [1] "Brain"      "Heart"      "Cerebellum" "Kidney"     "Liver"      "Ovary"      "Testis"


# ---
# -- CS - Carnegie Stage to weeks
cs.to.days = c('CS13' = 32, 'CS14' = 33, 'CS15' = 36, 'CS16' = 39, 'CS17' = 41, 'CS18' = 44, 
               'CS19' = 46, 'CS20' = 49, 'CS21' = 51, 'CS22' = 53, 'CS23' = 56) # dpc
# CS to weeks
cs.to.w = round(cs.to.days/4)
cs.to.w = paste(cs.to.w, "w", sep="") # adding w, that's how all values are stored in age column age+measure
names(cs.to.w) = names(cs.to.days) # asigning CS stages in correspondance to each week

# replacing CS with w
cs.idx = which(rse.jxn.cytosk@colData$age %in% names(cs.to.w)) # index of cs values to replace
cs.values = rse.jxn.cytosk@colData$age[cs.idx] # selecting all cs values
rse.jxn.cytosk@colData$age = replace(rse.jxn.cytosk@colData$age, cs.idx, cs.to.w[cs.values]) 

# separating age and measure into two columns
# replace w (weeks) with wpc (weeks post conseption) (w is wpc according to figure 1 in doi: 10.1038/s41586-019-1338-5)
rse.jxn.cytosk@colData$age = gsub('w$', 'wpc', rse.jxn.cytosk@colData$age)
rm(cs.idx, cs.values, cs.to.w, cs.to.days)

# --
# getting ordered list of ages
measures_order <- c("wpc", "dpb", "mpb", "ypb")
unique_ages = unique(rse.jxn.cytosk@colData$age)
ages_order = c()
for (i in c("wpc", "dpb", "mpb", "ypb")){
  time_group_values = grep(paste0('^', i, '|', i, '$'), unique_ages, value=TRUE)
  period_numbers = as.numeric(gsub("\\D", "", time_group_values))
  ages_order = c(ages_order,
                 time_group_values[order(period_numbers)])
}
ages_order
rm(i, unique_ages)

# merging age groups
fetus = c("4wpc", "5wpc", "6wpc", "7wpc", "8wpc", "9wpc",  "10wpc", "11wpc", "12wpc", "13wpc",
          "14wpc", "16wpc", "18wpc", "19wpc", "20wpc")
newborn = c('0dpb', '4dpb', '6dpb', '15dpb', '18dpb', '19dpb', '34dpb') # typically refers to a baby from birth up to about 2 months of age
infant = c('94dpb', '127dpb', '221dpb', '226dpb', '270dpb', '271dpb', '6mpb', '1ypb') # infant - child before the age of 12 months
toddler = c('2ypb') # child between 2 and 3 years
school = c('4ypb', '7ypb','8ypb')
teen = c("13ypb","14ypb","16ypb","17ypb" ) # 13-19 years
legal.age = c("25ypb", "28ypb", "29ypb", "32ypb", "39ypb",  "46ypb",  "50ypb",  "53ypb", 
          "54ypb", "55ypb",  "58ypb")
adult = c(toddler, school, teen, legal.age)

# -- Adding column age group
# removing ovary samples, since there are not enough samples for statistical analysis
rse.jxn.cytosk = rse.jxn.cytosk[,rse.jxn.cytosk@colData$tissue != 'Ovary']

# fetus
# heart develops earlier than other organs, so period before birth was limited to first 10 weeks
rse.jxn.cytosk@colData$age_group[rse.jxn.cytosk@colData$tissue != 'Heart' &
                                 rse.jxn.cytosk@colData$age %in% fetus] = "fetus"

rse.jxn.cytosk@colData$age_group[rse.jxn.cytosk@colData$tissue == 'Heart' &
                                rse.jxn.cytosk@colData$age %in% 
      c("4wpc", "5wpc",  "6wpc" , "7wpc",  "8wpc",  "9wpc",  "10wpc")] = "fetus"

# adult
rse.jxn.cytosk@colData$age_group[!(rse.jxn.cytosk@colData$tissue %in% c('Testies', 'Kidney')) &
                                 rse.jxn.cytosk@colData$age %in% adult] = "adult"

# including infant in adult group for kidney, since otherwise there are not enough of samples
rse.jxn.cytosk@colData$age_group[rse.jxn.cytosk@colData$tissue == 'Kidney' &
    rse.jxn.cytosk@colData$age %in% c(infant, adult)] = "adult"

# excluding teen from Testis, considering developmental peculiarities of testis development
rse.jxn.cytosk@colData$age_group[rse.jxn.cytosk@colData$tissue == 'Testis' &
                                 rse.jxn.cytosk@colData$age %in% legal.age] = "adult"

rse.jxn.cytosk@colData[100,c('tissue', 'age', "age_group")]
# leaving only fetus and adult age groups
rse.jxn.cytosk = rse.jxn.cytosk[,rse.jxn.cytosk@colData$age_group %in% c('fetus','adult')]

# save rse file with formatted annotation
saveRDS(rse.jxn.cytosk,'rse.jxn.cytosk.rds')

# read rse file
rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)
