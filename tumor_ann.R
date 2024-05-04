source('ages_comparison.R')

# cancer significant events found by SAJR
cancer.sign.jxns = c('chr13:26671990-26681053',
                     'chr12:110436205-110437083',
                     'chr7:152801732-152814549',
                     'chr7:152759927-152800530',
                     'chr7:152820348-152823341',
                     'chr10:26751784-26765217',
                     'chr10:26770346-26777064',
                     'chr15:22914883-22917787')

# signigicant events in tissue and development

sign.jxns.info.dev.and.cancer = list()
for (tissue in names(output)){
  output.tissue = output[[tissue]]
  print(tissue)
  for (tools.found in names(output.tissue)){
    print(tools.found)
    output.tissues.sign = output.tissue[[tools.found]]
    sign = output.tissues.sign[output.tissues.sign$unifiedJxnID %in% cancer.sign.jxns,]
    if (nrow(sign)==0){
      next
    }
    if (tools.found == 'all.common.jxns'){
      print(table(cancer.sign.jxns %in% unique(output.tissues.sign$unifiedJxnID) ))
      print(cancer.sign.jxns[!(cancer.sign.jxns %in% unique(output.tissues.sign$unifiedJxnID))])
      next
    }
    sign$tissue = tissue
    sign$tools_found = tools.found

    sign <- sign[order(abs(sign$dPSI), decreasing = TRUE), ]
    sign <- sign[!duplicated(sign$unifiedJxnID), ]

    
    sign.jxns.info.dev.and.cancer[[paste0(tissue, ".", tools.found)]] = sign
  }
}
sign.jxns.info.dev.and.cancer = do.call(rbind, sign.jxns.info.dev.and.cancer)
sign.jxns.info.dev.and.cancer = 
  sign.jxns.info.dev.and.cancer[order(sign.jxns.info.dev.and.cancer$GeneID, 
                                      sign.jxns.info.dev.and.cancer$unifiedJxnID,
                                      sign.jxns.info.dev.and.cancer$tissue), ]

#write.csv(sign.dev.and.cancer, file = "my_data.csv", row.names = FALSE)
