#source('ages_comp_trying2.0.R')
#install.packages("RColorBrewer")
#dev.off()
library(RColorBrewer)
# display.brewer.pal(n = 8, name = "Spectral")
# set2_colors <- brewer.pal(n = 8, name = "Spectral")
# set2_colors

# легенда
# корреляционная кривая ?
# группы в djexpress?
# дубликаты значимые только в djexpress? раздваиваются?
# corr coef для p value
# дать название всему графику

# сделать таблицу значимых во всех методах
# сделать реверсед таблицу

prepareDataBarplot = function(all.common.jxns, sign.all.tools, 
                              sajr.dje.sign, sajr.diego.sign, dje.diego.sign,
                              sajr.only.sign, dje.only.sign, diego.only.sign){
  data.frame(
    software = c("DJexpress", "SAJR", "DIEGO"),
    
    # -- common junctions
    # junctions, significant in both tools (green)
    sign.in.all.tools = c(nrow(sign.all.tools), 
                           nrow(sign.all.tools),
                           nrow(sign.all.tools)),
    
    # significant junctions only in one method from common junctions (red)
    sign.in.sajr.dje = c(nrow(sajr.dje.sign), 
                         nrow(sajr.dje.sign),
                           0),    
    
    sign.in.sajr.diego = c(0, 
                           nrow(sajr.diego.sign),
                           nrow(sajr.diego.sign)),  
    
    sign.in.dje.diego = c(nrow(dje.diego.sign), 
                           0,
                           nrow(dje.diego.sign)),    
    
    
    dje.sajr.only.sign=  c(nrow(dje.only.sign),
                           nrow(sajr.only.sign),
                           nrow(diego.only.sign)) # all signifficant in dje in common, minus sign in both tools
  )
}

setPlotParameters <- function() {
  #Units: Margins are specified in lines of text, which can vary depending on the font size and graphical device.
  par(tcl=-0.2)
  
  # Purpose: mar controls the margins around the plot region itself. These margins define the empty space between the plot area (where data is plotted) and the edges of the plotting device or figure region.
  par(mar = c(2, 0, 0, 0) + 0.1)  # mar = c(bottom, left, top, right)
  
  #Purpose: par(oma) controls the outer margins around the entire figure region, including all plots and axes when you have multiple plots arranged in a grid.
  par(oma = c(0, 2, 1, 0) + 0.1)  # Bottom, left, top, right
  
  #Purpose: mgp controls the distance between the axis title, axis labels, and the axis line. This indirectly affects the space around the plot area, especially when you have long axis labels or titles.
  # Values: mgp takes a numeric vector of length 3:
  #  mgp = c(axis title distance, axis label distance, axis line distance)
  par(mgp = c(2, 1, 0))  # Title, label, and line distances 
  
  # Purpose: plt defines the plot region within the figure region as fractions of the device's width and height. By adjusting these fractions, you can indirectly control the spacing around the plot area.
  # plt = c(x1, x2, y1, y2)
  # where x1 and x2 are the horizontal limits (left and right) and y1 and y2 are the vertical limits (bottom and top) of the plot region.
  
  # Set up multi-plot layout
  #par(mfrow = c(n_rows, n_cols))
  par(xpd=TRUE)
  
  par(pty = "s")
  
  par(cex.axis = 0.7)  # Increase font size by 50%
  
  # Create the initial matrix
  layout_matrix = matrix(1:42, nrow = 6, ncol = 7, byrow = TRUE)
  # Combine the matrix and the new row
  layout_matrix = rbind(layout_matrix, c(43, 43, 43))
  layout(mat = layout_matrix)
  
}

makeBarplot = function(tissue,df){
  # Create the vector of colors for each quarter
  data = t(as.matrix(df[, -1]))
  
  print(cbind(data[,1], c(0,0), c(0,0)))
  
  bp = barplot(
    cbind(data[,1], c(0,0), c(0,0)), 
    # main = pair.name,
    ylab = "# jxns",
    col = c("#4DAF4A", "#377EB8", '#F781BF', "#A65628","#984EA3"),
    beside = F,
    ylim = c(0, 250),
             # max(sapply(output, 
             #            function(tissue) {
             #              sapply(
             #                tissue[c('sajr.only.significant.events', 
             #                              'dje.only.significant.events')], 
             #                     nrow
             #                )}
             #              ))),
             # 
    legend = FALSE  # Disable the legend creation within each subplot
  )
  
  barplot(
    cbind(c(0,0), data[,2], c(0,0)), 
    # main = pair.name,
    col = c("#4DAF4A", "#377EB8", '#F781BF', "#A65628", "#FF7F00"),
    beside = F,
    add=TRUE,
  )
  
  barplot(
    cbind(c(0,0), c(0,0), data[,3]), 
    # main = pair.name,
    col = c("#4DAF4A", "#377EB8", '#F781BF', "#A65628",'#FFFF33' ),
    beside = F,
    add=TRUE,
  )
  
  if (par("mfg")[1]==6){
    axis(1, at = bp, labels = df$software, cex=0.8,las=3)
  }
  # assigns tissue name for each row of graphs
  mtext(tissue, side = 2, line = 3, las = 0)  
}


plot_graphs = function(all.common.jxns, sign.all.tools, 
                       sajr.dje.sign, sajr.diego.sign, dje.diego.sign,
                       sajr.only.sign, dje.only.sign, diego.only.sign,
                       par.1, par.2){
  # dpsi vs logFC
  dict = list('dPSI' = c(-1,1), 
              'logFC' = c(-4,4), 
              'abundance_change' = c(-3,3))
  
  dict.ticks = list('dPSI' = seq(-1,1,by=0.5), 
                    'logFC'=seq(-4,4,by=2),
                    'abundance_change' = seq(-3,3,by=1.5))
  
  # ??? как то переписать красиво
  if (par.1 %in% names(dict.ticks)) log='' else log='xy'
  
  # как от этого избавиться?
  plot(all.common.jxns[,par.1], all.common.jxns[,par.2],
       xaxt = "n", yaxt='n',    # Suppress x-axis ticks and labels
       xlab = '', ylab = par.2,
       xlim=dict[[par.1]], ylim=dict[[par.2]],
       
       log=log,
       type = 'p', col = 'lightgrey',
       pch = 16 , cex = 0.7, lwd = 1,
       bty = "L" # Set bty = "L" for left and bottom lines only
       ) 
  
  # Calculate Spearman correlation
  # result = cor.test(all.common.jxns[,par.1], all.common.jxns[,par.2], method = "spearman")
  # corr.coef = round(result$estimate, digits=2)
  # text(1, -3, corr.coef, col = "black", cex=0.8)
  
  
  # Order data by x ranks
  #ordered_data <- data.frame(x = all.common.jxns[,par.1], 
                            # y = all.common.jxns[,par.2])[order(all.common.jxns[,par.1]), ]
  #print(ordered_data)
  
  #lines(lowess(ordered_data$x, ordered_data$y), col = "red")
  

  # axis. Only for last row of graphs x axis is assigned with values and axis name
  # axcept for p.values - x axis values are assigned for each row, but not the name
  if(par.1 %in% names(dict.ticks)){
    axis(1, at=dict.ticks[[par.1]], labels = FALSE)
    axis(2, at=dict.ticks[[par.2]], labels = TRUE)
    if (par("mfg")[1]==6){
      axis(1, at=dict.ticks[[par.1]], labels = TRUE)
      axis(2, at=dict.ticks[[par.2]], labels = TRUE)
      title(xlab = par.1, xpd = NA)
      }
  } else {
    axis(1, labels = TRUE)
    axis(2, labels = TRUE)
    if (par("mfg")[1]==6) {
      title(xlab = par.1, xpd = NA)
      }
  }
  
  # Highlight points belonging to "common"
  
  # only sajr significant events
  points(
    sajr.only.sign[,par.1],
    sajr.only.sign[,par.2],
    col = "#FF7F00",
    pch = 16)

  # only dje significant events
  points(
    dje.only.sign[,par.1],
    dje.only.sign[,par.2],
    col = "#984EA3",
    pch = 16)
  
  points(
    diego.only.sign[,par.1],
    diego.only.sign[,par.2],
    col = '#FFFF33',
    pch = 16)
  
  points(
    sajr.dje.sign[,par.1],
    sajr.dje.sign[,par.2],
    col = "#377EB8",
    pch = 16)
  
  # only dje significant events
  points(
    sajr.diego.sign[,par.1],
    sajr.diego.sign[,par.2],
    col = "#F781BF",
    pch = 16)
  
  points(
    dje.diego.sign[,par.1],
    dje.diego.sign[,par.2],
    col = '#A65628',
    pch = 16)
  
  points(sign.all.tools[,par.1], 
         sign.all.tools[,par.2],
         col = "#4DAF4A",
         pch = 16)
  
  
  labels <- sign.all.tools$GeneID
  # Add text labels to the points
  n_labels <- min(5, nrow(sign.all.tools))
  
  # if (n_labels > 0) {
  #   # Determine label positions based on row number (odd/even)
  #   label_positions <- ifelse(seq_len(n_labels) %% 2 == 1, 4, 3)  # 4 for right, 3 for above
  #   text(sign.both.tools[1:n_labels, par.1], sign.both.tools[1:n_labels, par.2], 
  #        labels = labels[1:n_labels], pos = label_positions, cex = 0.8)
  # }
  
  if (n_labels > 0) {
    # Calculate label positions based on modulo 3
    label_positions <- (seq_len(n_labels) - 1) %% 3 + 1
    label_positions <- c(4, 3, 1)[label_positions]  # Map 1 to right, 2 to below, 3 to above

    text(sign.all.tools[1:n_labels, par.1], sign.all.tools[1:n_labels, par.2],
         labels = labels[1:n_labels], pos = label_positions, cex = 0.8)
  }

}

# SUBSETTING?


# addGeneNames = function(all.common.jxns){
#   print(all.common.jxns)
#   
# }


#makePlots = function(){
setPlotParameters()

Map(function(output.tissue, tissue) {
    all.common.jxns = output.tissue$all.common.jxns
    sajr.dje.sign = output.tissue$sajr.dje.significant.events
    sajr.diego.sign = output.tissue$sajr.diego.significant.events
    dje.diego.sign = output.tissue$dje.diego.significant.events
    sajr.only.sign = output.tissue$sajr.only.significant.events
    dje.only.sign = output.tissue$dje.only.significant.events
    diego.only.sign = output.tissue$diego.only.significant.events
    sign.all.tools = output.tissue$sign.all.tools
    
    print('here')
    print(max(all.common.jxns$abundance_change, na.rm = TRUE))
    print('end')
  
    df = prepareDataBarplot(all.common.jxns, sign.all.tools, 
                            sajr.dje.sign, sajr.diego.sign, dje.diego.sign,
                            sajr.only.sign, dje.only.sign, diego.only.sign)
    #print(all.common.jxns)
    
    print(df)
    makeBarplot(tissue,df)
    plot_graphs(all.common.jxns, sign.all.tools, 
                sajr.dje.sign, sajr.diego.sign, dje.diego.sign,
                sajr.only.sign, dje.only.sign, diego.only.sign, 'dPSI', 'logFC')   #, xlim==c(-1,1), ylim=c(-4,4))
    plot_graphs(all.common.jxns, sign.all.tools, 
                sajr.dje.sign, sajr.diego.sign, dje.diego.sign,
                sajr.only.sign, dje.only.sign, diego.only.sign, 'dPSI', 'abundance_change')   #, xlim==c(-1,1), ylim=c(-4,4))
    plot_graphs(all.common.jxns, sign.all.tools, 
                sajr.dje.sign, sajr.diego.sign, dje.diego.sign,
                sajr.only.sign, dje.only.sign, diego.only.sign, 'logFC', 'abundance_change')   #, xlim==c(-1,1), ylim=c(-4,4))
    plot_graphs(all.common.jxns, sign.all.tools, 
                sajr.dje.sign, sajr.diego.sign, dje.diego.sign,
                sajr.only.sign, dje.only.sign, diego.only.sign, 'p.value.sajr', 'P.Value')   #, xlim==c(-1,1), ylim=c(-4,4))
    plot_graphs(all.common.jxns, sign.all.tools, 
                sajr.dje.sign, sajr.diego.sign, dje.diego.sign,
                sajr.only.sign, dje.only.sign, diego.only.sign, 'p.value.sajr', 'p_val')   #, xlim==c(-1,1), ylim=c(-4,4))
    plot_graphs(all.common.jxns, sign.all.tools, 
                sajr.dje.sign, sajr.diego.sign, dje.diego.sign,
                sajr.only.sign, dje.only.sign, diego.only.sign, 'P.Value', 'p_val')   #, xlim==c(-1,1), ylim=c(-4,4))


    
  },
  output,
  names(output)
  )

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend('center', inset = c(0, 0),  # Adjust inset as needed
         legend = c(
                    "signigicant in DJE and SAJR", "signigicant SAJR and DIEGO","signigicant DJE and DIEGO", 
                    'significant only in DJE', 'significant only in SAJR', 'significant only in DIEGO',
                    "significant in all tools", 
                    "number - spearman correlation"),
         xpd=TRUE,
         col = c("#377EB8", '#F781BF', "#A65628","#984EA3", "#FF7F00",'#FFFF33',"#4DAF4A",  NA),
         pch=20,
         bty='n',
       horiz = FALSE,
       ncol=3,
       pt.cex=3)

#}

# подписать гены
# подписать корреляции
# подписать число джанкшенов в барплотах
# добавить легенды
# можно ли в барплотах сделать разные раскраски для разных столбцов?
# проверить что с дубликатами в значимых


# setPlotParameters()
# makePlots()

