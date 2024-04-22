#source('ages_comparison.R')


cancer.sign.jxns = c('chr13:26671990-26681053',
                     'chr12:110436205-110437083',
                     'chr7:152801732-152814549',
                     'chr7:152759927-152800530',
                     'chr7:152820348-152823341',
                     'chr10:26751784-26765217',
                     'chr10:26770346-26777064',
                     'chr15:22914883-22917787')

'chr13:26671990-26681053' %in% cancer.sign.jxns

#output
for (tissue in unique.tissues){
  print(tissue)
  sign.all = na.omit(output[[tissue]]$sign.all.tools)
  #print(sign.all)
  all = sign.all[sign.all$unifiedJxnID %in% cancer.sign.jxns,]
  if (nrow(all)>0){
    print('all')
    print(all)
  }
  
  print("")
  sign.sajr.dje = na.omit(output[[tissue]]$sajr.dje.significant.events)
  sajr.dje = sign.sajr.dje[sign.sajr.dje$unifiedJxnID %in% cancer.sign.jxns,]
  if (nrow(sajr.dje)>0){
    print('sajr.dje')
    print(sajr.dje)
  }
  
  print("")
  sign.dje.diego = na.omit(output[[tissue]]$dje.diego.significant.events)
  dje.diego = sign.dje.diego[sign.dje.diego$unifiedJxnID %in% cancer.sign.jxns,]
  if (nrow(dje.diego)>0){
    print('dje.diego')
    print(dje.diego)
  }
  
  print("")
  sign.sajr = na.omit(output[[tissue]]$sajr.only.significant.events)
  sajr = sign.sajr[sign.sajr$unifiedJxnID %in% cancer.sign.jxns,]
  if (nrow(sajr)>0){
    print('sajr.only')
    print(sajr)
  }
  
  print("")
  sign.dje = na.omit(output[[tissue]]$dje.only.significant.events)
  dje = sign.dje[sign.dje$unifiedJxnID %in% cancer.sign.jxns,]
  if (nrow(dje)>0){
    print('dje.only')
    print(dje)
  }
  
  print("")
  sign.diego = na.omit(output[[tissue]]$diego.only)
  diego = sign.diego[sign.diego$unifiedJxnID %in% cancer.sign.jxns,]
  if (nrow(diego)>0){
    print('diego.only')
    print(diego)
  }
}


sign.tissue = lapply(output, function(tissue) print(na.omit(tissue$sign.all.tools)))
sign.tissue

sign.tissue.filt <- Filter(function(df) nrow(df) > 0, sign.tissue)

# 1. Add key column to each dataframe
s <- lapply(names(sign.tissue.filt), function(tissue) {
  df <- sign.tissue.filt[[tissue]]
  df$tissue <- tissue
  return(df)
})

s

# 2. Merge dataframes using full join
merged_df <- Reduce(function(x, y) merge(x, y, all = TRUE), s)


# 3. (Optional) Move the key column to the first position
merged_df <- merged_df[, c("GeneID", setdiff(names(merged_df), "GeneID"))]

df = merged_df[,c('GeneID', 'unifiedJxnID', 'dPSI', 'FDR.sajr', 'logFC', 'FDR', 'abundance_change', 'p_val', 'tissue')]

df[order(-df$dPSI, -df$logFC), ] # Sort descending

write.csv(df, file = "my_data.csv", row.names = FALSE)
