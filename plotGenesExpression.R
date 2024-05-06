source("rseAnnotationPreprocessing.R")

# -----------------------------------------------------------------------------------
# ------------------------- samples occurrence heatmap ------------------------------
# -----------------------------------------------------------------------------------
plotHeatmapSamples = function(){
  par(oma = c(2, 2, 0, 0))  # bottom, left, top, right.
  # occurrence of samples
  sample_occurance = as.data.frame.matrix(
    table(rse.gene.cytosk@colData$tissue, rse.gene.cytosk@colData$age_group_specific) )
  # The table() function takes these two vectors as input and creates a contingency table. This table shows the frequency distribution of cells across different combinations of tissue types and age groups.
  sample_occurance = sample_occurance[,names(order)]
  sample_occurance = t(sample_occurance)
  # # making a column for tissues
  # sample_occurance$tissue = rownames(sample_occurance)
  # sample_occurance
  
  # Set color palette (e.g., blue to red)
  colors <- colorRampPalette(c("white", "red"))(256)
  
  # Create the heatmap
  image(1:ncol(sample_occurance), 1:nrow(sample_occurance), 
        t(sample_occurance), col = colors, axes = FALSE, xlab = "", ylab = "")
  
  # Add text labels with data values
  text(x = col(sample_occurance), 
       y = row(sample_occurance), 
       labels = sample_occurance, 
       col = "black")
  
  # Add axes and labels
  axis(1, at = 1:ncol(sample_occurance), labels = colnames(sample_occurance), las = 2)
  axis(2, at = 1:nrow(sample_occurance), labels = rownames(sample_occurance), las = 2)
}

# ----------------------------------------------------------------------------------
# ------------------------ gene expression vs time graphs -------------------------
# ----------------------------------------------------------------------------------

findYlim = function(rse.gene.cytosk){
  cpm = as.data.frame(t(rse.gene.cytosk@assays@data$cpm))
  cpm$tissue = rse.gene.cytosk@colData$tissue
  
  a = by(cpm[, sapply(cpm, is.numeric)], cpm$tissue, 
        function(x) as.data.frame(
          apply(x, 2, quantile, probs = c(0.25, 0.75)), simplify = FALSE))
  a = do.call(cbind, a)
  max.up.quantile = max.col(a)[2]
  cpm.max.quantiles = a[,max.up.quantile]
  ylim.max = cpm.max.quantiles[2]+1.5*(cpm.max.quantiles[2]-cpm.max.quantiles[1])
  c(0,ylim.max)
}

fitCurve = function(gene.id, age.group.specific, gene.tissue.counts, col){
  # Fit a curve (example using loess)
  model = loess(gene.tissue.counts[[gene.id]] ~ as.numeric(age.group.specific))
  # Generate x-values for prediction
  x_pred = seq(min(age.group.specific), max(age.group.specific), length.out = 100)
  # Predict y-values using the model
  y_pred = predict(model, newdata = data.frame(age.group.specific = x_pred))
  # Add the curve to the plot
  lines(x_pred, y_pred, col = col, lwd = 2)  # Adjust lwd for line thickness
}

createGrid = function(gene.rse){
  numb.graphs = round(sqrt(length(gene.rse@rowRanges)))
  numb.empty = numb.graphs**2 - length(gene.rse@rowRanges)
  list(numb.graphs=numb.graphs, numb.empty=numb.empty)
}

setParams = function(gene.rse){
  numb = createGrid(gene.rse)
  par(mfrow = c(numb$numb.graphs, numb$numb.graphs),
      oma = c(4, 3, 1, 1),  # bottom, left, top, right.
      mar = c(1, 2, 1, 1),  # bottom, left, top, right.
      cex.axis = 0.7,  # Adjust the value as needed
      bty="l")
}

setAxis = function(x.value, gene.name, gene.rse){
  numb = createGrid(gene.rse)
  
  title(main = gene.name, line = 0.2, cex=0.8)  # Place title 3 lines above the plot
  axis(1, at = 1:length(x.value), labels = FALSE)  # x-axis ticks
  axis(2, labels = FALSE)  # y-axis ticks
  
  boxplot.coord = par("mfg")
  # Add y-axis only for rightmost plots
  if (boxplot.coord[2] == 1) {
    axis(2, las = 1)  
    mtext("CPM", side = 2, line = 2.2, las = 0, cex = par("cex.axis"))
  }
  # Add x-axis only for bottom plots
  if ((boxplot.coord[1] == 5) |
      (boxplot.coord[2] == numb$numb.graphs & boxplot.coord[1] == (numb$numb.graphs-numb$numb.empty)) ){
    axis(1, at = 1:length(x.value), labels = x.value, las = 2) 
  }
  grid(nx = NULL, ny = NULL)
}


plotScatterplotExpression = function(gene.rse, tissue.col){
  setParams(gene.rse)
  
  gene.ids.ordered = gene.rse@rowRanges[order(gene.rse@rowRanges$gene_name),
                                        c('gene_id', 'gene_name')]
  gene.ids.ordered = setNames(gene.ids.ordered$gene_id, gene.ids.ordered$gene_name)
  tissues = unique(gene.rse@colData$tissue)
  
  for (gene.name in names(gene.ids.ordered)){
    # Create a new plot for each gene
    gene.id = gene.ids.ordered[gene.name]
    plot(
      0, 0,  # Placeholder values, will be replaced by actual data
      xlim = c(1, length(order)),  # Set x-axis limits based on number of tissues
      ylim = findYlim(rse.gene.cytosk),  # Set y-axis limits based on gene expression
      type = "n",  # Start with an empty plot
      xaxt = "n",  # Suppress default x-axis
      yaxt = "n",  # Suppress default x-axis
      # main = gene.rse@rowRanges[gene.rse@rowRanges$gene_id==gene.id,]$gene_name  # Set title to gene name
    )
    setAxis(names(order), gene.name, gene.rse)
    
    for (tissue in tissues){
      samples = rownames(gene.rse@colData[gene.rse@colData$tissue==tissue,])
      col = tissue.col[tissue]
      age.group.specific = 
        gene.rse@colData[gene.rse@colData$tissue==tissue,]$age_group_specific
      age.group.specific = order[age.group.specific]
      
      gene.tissue.counts = as.data.frame(t(gene.rse@assays@data$cpm[gene.id,samples,drop=F]))
      gene.tissue.counts$age.group.specific = age.group.specific
      
      points(age.group.specific, gene.tissue.counts[[gene.id]], col=col)
      fitCurve(gene.id, age.group.specific, gene.tissue.counts, col=col)
    }
  }
  # Plot the legend in the last cell
  plot(x=0, y=0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend('bottom', legend = tissues, col = tissue.col, pch = 16,
         bty = "n", y.intersp = 0.6, xpd = TRUE,
         inset = c(0, -0.2))
}


# boxplots
plotBarplotExpression = function(gene.rse, tissue.col...){
  # Box: The box represents the interquartile range (IQR), which contains the middle 50% of the data. The bottom and top edges of the box correspond to the first quartile (Q1) and third quartile (Q3), respectively.
  # Median Line: A horizontal line inside the box that marks the median (Q2) of the data.
  # Whiskers: Lines extending from the box that represent the range of the data, excluding outliers.
  
  cpm = as.data.frame(t(gene.rse@assays@data$cpm))
  cpm$tissue = gene.rse@colData$tissue
  cpm$tissue = factor(cpm$tissue)
  
  gene.ids.ordered = gene.rse@rowRanges[order(gene.rse@rowRanges$gene_name),
                                        c('gene_id', 'gene_name')]
  gene.ids.ordered = setNames(gene.ids.ordered$gene_id, gene.ids.ordered$gene_name)
  tissues = unique(gene.rse@colData$tissue)
  
  # Create boxplot
  for (gene.name in names(gene.ids.ordered)){
    gene.id = gene.ids.ordered[gene.name]
    boxplot(cpm[[gene.id]] ~ tissue, data = cpm,
            xlab = "Tissue", ylab = "CPM",
            ylim = findYlim(rse.gene.cytosk),
            xaxt = "n", yaxt = "n", 
            col=tissue.col)
    setAxis(tissues, gene.name, gene.rse)
  }
}


order = c("4wpc", "5wpc", "6wpc", "7wpc", "8wpc", "9wpc",  "10wpc", "11wpc", "12wpc", "13wpc",
          "14wpc", "16wpc", "18wpc", "19wpc", "20wpc", "newborn", "infant", "toddler", "school", 
          "teen", "25-35 y.o.", "36-45 y.o.", "46-45 y.o.", "56-55 y.o.") # setting new order
order = setNames(1:length(order), order)

# --- plotting
tissue.col=c(Brain="#3399CC",
             Cerebellum="#33CCFF",
             Heart="#CC0000",
             Kidney="#CC9900",
             Liver="#339900",
             Ovary="#CC3399",
             Testis="#FF6600")

plotHeatmapSamples()
plotScatterplotExpression(rse.gene.cytosk, tissue.col)
plotBarplotExpression(rse.gene.cytosk, tissue.col)

# gene.ids.ordered = gene.info[order(gene.info$gene_name),]$gene_id
# tissues = unique(gene.rse@colData$tissue)
# ylim=findYlim(rse.gene.cytosk)
# 
# cpm = as.data.frame(t(gene.rse@assays@data$cpm))
# cpm$tissue = gene.rse@colData$tissue
# cpm$tissue = factor(cpm$tissue)

