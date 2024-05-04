source("data.R")
# rse.ERP109002.jxn.cytosk.genes = readRDS('rse.ERP109002.jxn.cytosk.genes.rds', refhook = NULL)
# 
# # metadata
# # renaming to shorten the name
# rse = rse.ERP109002.jxn.cytosk.genes
# 
# # from sample metadata extracting sra.library_name, containing information on tissue, age etc
# split_cols = strsplit(rse@colData$sra.library_name, "\\.")
# rse@colData = within(rse@colData, {
#   sample_id <- sapply(split_cols, "[[", 1)
#   organism <- sapply(split_cols, "[[", 2)
#   tissue <- sapply(split_cols, "[[", 3)
#   age <- sapply(split_cols, "[[", 4)
#   gender <- sapply(split_cols, "[[", 5)
# })
# 
# # combining some tissues
# rse@colData$tissue = sub("Forebrain", "Brain", rse@colData$tissue)
# rse@colData$tissue = sub("Hindbrain", "Cerebellum", rse@colData$tissue)
# rse@colData$tissue = sub("KidneyTestis", "Kidney", rse@colData$tissue)
# unique(rse@colData$tissue)
# # [1] "Brain"      "Heart"      "Cerebellum" "Kidney"     "Liver"      "Ovary"      "Testis"    
# 
# # age measures
# #unique(gsub("[^[:alpha:]]", "", rse@colData$age))
# 
# # order of columns (age)
# measures_order <- c("CS", "w", "dpb", "mpb", "ypb")
# unique_ages = unique(rse@colData$age)
# ages_order = c()
# for (i in c("CS", "w", "dpb", "mpb", "ypb")){
#   time_group_values = grep(paste0('^', i, '|', i, '$'), unique_ages, value=TRUE)
#   period_numbers = as.numeric(gsub("\\D", "", time_group_values))
#   ages_order = c(ages_order, 
#                  time_group_values[order(period_numbers)])
# }
# #print(paste0("'",ages_order, "'", collapse = ", "))
# #length(ages_order) == length(unique_ages)
# 
# # merging age groups
# # CS - Carnegie Stage 
# fetus = c('CS13', 'CS14', 'CS16', 'CS17', 'CS18', 'CS19', 'CS20', 'CS21', 'CS22', 'CS23', '4w', '5w', '6w', '7w', '8w', '9w', '10w', '11w', '12w', '13w', '16w', '18w', '19w', '20w')
# infant = c('0dpb', '4dpb', '6dpb', '15dpb', '18dpb', '19dpb', '34dpb', '94dpb', '127dpb', '221dpb', '226dpb', '270dpb', '271dpb', '6mpb', '1ypb') # infant - child before the age of 12 months
# toddler = c('2ypb') # child between 2 and 3 years
# school = c('4ypb', '7ypb','8ypb')
# teen = c("13ypb","14ypb","16ypb","17ypb" ) # 13-19 years
# yo_25_35 = c("25ypb","28ypb","29ypb","32ypb" )
# yo_36_45 = c("39ypb" )
# yo_46_55 = c("46ypb", "50ypb", "53ypb" , "54ypb", "55ypb")
# yo_56_65 = c("58ypb")
# all_age_gr = c(infant, toddler, school, teen, yo_25_35, yo_36_45, yo_46_55, yo_56_65)
# 
# # -- Adding column age group
# rse@colData$age_group[
#   rse@colData$tissue != 'Heart' &
#   rse@colData$age %in% fetus] = "fetus"
# 
# rse@colData$age_group[
#     rse@colData$tissue == 'Heart' &
#     rse@colData$age %in% c('CS13', 'CS14', 'CS16','4w', '5w', '6w', '7w', '8w', '9w', '10w')] = "fetus"
# 
# rse@colData$age_group[
#   !(rse@colData$tissue %in% c('Testies', 'Kidney')) &
#   rse@colData$age %in% c(toddler, school, teen, yo_25_35, yo_36_45, yo_46_55, yo_56_65)] = "adult"
# 
# # including infant into adult group for kidney, since otherwise there is not enoug of samples 
# rse@colData$age_group[
#   rse@colData$tissue == 'Kidney' &
#   rse@colData$age %in% c(infant, toddler, school, teen, yo_25_35, yo_36_45, yo_46_55, yo_56_65)] = "adult"
# 
# # excluding teen from Testies
# rse@colData$age_group[
#   rse@colData$tissue == 'Testis' &
#   rse@colData$age %in% c(yo_25_35, yo_36_45, yo_46_55, yo_56_65)] = "adult"
# 
# # removing ovary samples, since there are not enough samples for statistical analysis
# rse = rse[,rse@colData$tissue != 'Ovary']
# rse = rse[,rse@colData$age_group %in% c('fetus','adult')]
# 
# # renaming object back
# rse.ERP109002.jxn.cytosk.genes = rse
# 
# rm(rse, split_cols)
# 
# saveRDS(rse.ERP109002.jxn.cytosk.genes,'rse.ERP109002.jxn.cytosk.genes.rds')
rse.ERP109002.jxn.cytosk.genes = readRDS('rse.ERP109002.jxn.cytosk.genes.rds', refhook = NULL)
#rse@assays@data$counts) 