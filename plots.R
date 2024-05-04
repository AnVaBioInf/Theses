source('ages_comp_trying2.0.R')

prepareDataBarplot = function(all.common.jxns, 
                              sign.both.tools, 
                              dje.sign.common, sajr.sign.common,
                              dje.not.common.jxns, sajr.not.common.jxns,
                              dje.sign.not.common, sajr.sign.not.common){
  data.frame(
    software = c("DJexpress", "SAJR"),
    # -- common junctions
    # junctions, significant in both tools (green)
    sign.in.both.tools = c(nrow(sign.both.tools), 
                           nrow(sign.both.tools)),
    
    # not significant from common junctions (dark grey)
    dje.sajr.not.sign.common = c(nrow(all.common.jxns)-nrow(dje.sign.common),
                                 nrow(all.common.jxns)-nrow(sajr.sign.common)),
    
    # significant junctions only in one method from common junctions (red)
    dje.sajr.sign.common = c(nrow(dje.sign.common) - nrow(sign.both.tools), # all signifficant in dje in common, minus sign in both tools
                             nrow(sajr.sign.common) - nrow(sign.both.tools)), # all signifficant in dje in common, minus sign in both tools
    
    # -- not common junctions
    # not significant, not from common junctions
    dje.sajr.not.sign.not.common = c(nrow(dje.not.common.jxns) - nrow(dje.sign.not.common),
                                     nrow(sajr.not.common.jxns) - nrow(sajr.sign.not.common)),
    
    # sign only in one method, not from common junctions (red)
    dje.sajr.sign.not.common = c(nrow(dje.sign.not.common), nrow(sajr.sign.not.common))
  )
}


makeBarplot = function(tissue,df){
  # Create the vector of colors for each quarter
  colors <- c("green", "grey", "red","grey","red")
  barplot(
    t(as.matrix(df[, -1])), 
    # main = pair.name,
    ylab = "# jxns",
    col = colors,
    beside = F,
    #log='y',
 #   ylim = c(0, max(nrow(df$all.sajr),nrow(df$all.dje))),
    #names.arg = sales$software,  # Specify the labels for each subplot
    legend = FALSE  # Disable the legend creation within each subplot
  )
  if (par("mfg")[1]==6){
    axis(1, at = c(1,2), labels = df$software)
  }
  # assigns tissue name for each row of graphs
  mtext(tissue, side = 2, line = 3, las = 0)  
}

setPlotParameters = function(){
  # stacked barplot
  par(mar = c(1, 3, 1, 3) + 0.1)
  # par(mar = c(bottom_margin, left_margin, top_margin, right_margin))
  par(oma = c(3, 3, 1, 1) + 0.1)  
  # par(mar = c(bottom_margin, left_margin, top_margin, right_margin))
  par(mgp = c(2, 1, 0)) 
  #par(mgp = c(title_distance, label_distance, line_distance))
  #par(xpd = TRUE)
  par(mfrow = c(length(unique.tissues)+1, 4))
}

plot_graphs = function(all.common.jxns, dje.only.sign, sajr.only.sign, sign.both.tools, par.1, par.2){
  # dpsi vs logFC
  dict = list('dPSI' = c(-1,1), 'logFC' = c(-4,4), 'logFC.sajr'=c(-4,4))
  dict.ticks = list('dPSI' = seq(-1,1,by=0.5), 'logFC.sajr'=seq(-4,4,by=2))
  # ??? как то переписать красиво
  if (par.1 %in% names(dict.ticks)) log='' else log='xy'
  
  # как от этого избавиться?
  plot(all.common.jxns[,par.1], all.common.jxns[,par.2],
       xaxt = "n",        # Suppress x-axis ticks and labels
       xlab = '', ylab = par.2,
       xlim=dict[[par.1]], ylim=dict[[par.2]],
       log=log
       ) 
  
  # axis. Only for last row of graphs x axis is assigned with values and axis name
  # axcept for p.values - x axis values are assigned for each row, but not the name
  if(par.1 %in% names(dict.ticks) & par("mfg")[1]!=6){
    axis(1, at=dict.ticks[[par.1]], labels = FALSE)
    if (par("mfg")[1]==6){
      axis(1, at=dict.ticks[[par.1]], labels = TRUE)
      title(xlab = par.1, xpd = NA)}
  } else {
    axis(1, labels = TRUE)
    if (par("mfg")[1]!=6) {
      title(xlab = par.1, xpd = NA)}
    
    # plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    # plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    # legend("center", legend = c("Common", "SAJR Only", "DJexpress Only"),
    #        col = c('green','pink','violet'), pch = 22,
    #        bty = "n",
    #        xpd = TRUE, horiz = TRUE, #inset = c(0.1, 0.1),
    #        cex=0.8)  # Adjust positioning as needed
  }
  
  # Highlight points belonging to "common"
  points(sign.both.tools[,par.1], 
         sign.both.tools[,par.2],
         col = "green",
         pch = 16)
  
  # only sajr significant events
  points(
    all.common.jxns[all.common.jxns$junctionID %in% sajr.only.sign$junctionID,][,par.1],
    all.common.jxns[all.common.jxns$junctionID %in% sajr.only.sign$junctionID,][,par.2],
    col = "pink",
    pch = 16)
  
  # only dje significant events
  points(
    all.common.jxns[all.common.jxns$junctionID %in% dje.only.sign$junctionID,][,par.1],
    all.common.jxns[all.common.jxns$junctionID %in% dje.only.sign$junctionID,][,par.2],
    col = "violet",
    pch = 16)
}

#makePlots = function(){
setPlotParameters()

Map(function(dje.output.tissue, sajr.output.tissue, significant.events.tissue, tissue) {
    print(tissue)
  
    all.sign.dje.tissue = significant.events.tissue$dje.significant.events
    all.sign.sajr.tissue = significant.events.tissue$sajr.significant.events

    print(c("all.sign.sajr.tissue",all.sign.sajr.tissue$junctionID))
    print(c("all.sign.dje.tissue",all.sign.dje.tissue$junctionID))
    
    all.common.jxns = merge(sajr.output.tissue, dje.output.tissue, by = "junctionID", all = FALSE)
    dje.not.common.jxns = dje.output.tissue[!dje.output.tissue$junctionID %in% all.common.jxns$junctionID,]
    sajr.not.common.jxns = sajr.output.tissue[!sajr.output.tissue$junctionID %in% all.common.jxns$junctionID,]
    print(c("all.common.jxns",all.common.jxns$junctionID))
    # print(c("sajr.not.common.jxns",head(sajr.not.common.jxns)))
    # print(c("dje.not.common.jxns",head(dje.not.common.jxns)))
    
    sign.both.tools = merge(all.sign.sajr.tissue, all.sign.dje.tissue, by = "junctionID", all = FALSE)
    dje.only.sign = all.sign.dje.tissue[!all.sign.dje.tissue$junctionID %in% all.sign.sajr.tissue$junctionID, ]
    sajr.only.sign = all.sign.sajr.tissue[!all.sign.sajr.tissue$junctionID %in% all.sign.dje.tissue$junctionID, ]
    print(c("sign.both.tools", sign.both.tools$junctionID))
    # print(c("sajr.only.sign",head(sajr.only.sign)))
    # print(c("dje.only.sign",head(dje.only.sign)))

    dje.sign.common = all.sign.dje.tissue[all.sign.dje.tissue$junctionID %in% all.common.jxns$junctionID,]
    sajr.sign.common = all.sign.sajr.tissue[all.sign.sajr.tissue$junctionID %in% all.common.jxns$junctionID,]
    dje.sign.not.common = all.sign.dje.tissue[all.sign.dje.tissue$junctionID %in% dje.not.common.jxns$junctionID,]
    sajr.sign.not.common = all.sign.sajr.tissue[all.sign.sajr.tissue$junctionID %in% sajr.not.common.jxns$junctionID,]
    print(c("sajr.sign.common",sajr.sign.common$junctionID))
    print(c("dje.sign.common",dje.sign.common$junctionID))
    # print(c("sajr.sign.not.common",head(sajr.sign.not.common)))
    # print(c("dje.sign.not.common",head(dje.sign.not.common)))
    # print(nrow(dje.sign.not.common))
    
    
    df = prepareDataBarplot(all.common.jxns, 
                                  sign.both.tools, 
                                  dje.sign.common, sajr.sign.common,
                                  dje.not.common.jxns, sajr.not.common.jxns,
                                  dje.sign.not.common, sajr.sign.not.common)
    #print(prepareDataBarplot)
  
    makeBarplot(tissue,df)
    plot_graphs(all.common.jxns, dje.only.sign, sajr.only.sign, sign.both.tools, 'dPSI', 'logFC')   #, xlim==c(-1,1), ylim=c(-4,4))
    plot_graphs(all.common.jxns, dje.only.sign, sajr.only.sign, sign.both.tools, 'logFC.sajr', 'logFC')   #, xlim==c(-1,1), ylim=c(-4,4))
    plot_graphs(all.common.jxns, dje.only.sign, sajr.only.sign, sign.both.tools, 'age', 'P.Value')   #, xlim==c(-1,1), ylim=c(-4,4))
  },
  dje.output,
  sajr.output,
  significant.events,
  names(significant.events)
  )
#}

setPlotParameters()
makePlots()

# добавить пересекающиеся
# исправить оси у барплота
# дать название всему графику
# проверить отметки (розовые итд)
# добавить легенды