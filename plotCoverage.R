library(recount3)
library(rtracklayer)

bigWig2Cov = function(bw){
  bw = as.data.frame(bw) # columns: seqnames, start, end, width, strand, score
  bw = bw[order(bw$start),] # ordering bw df my start coordinate column
  start = bw$start[1] # starting coordinate of the gene
  stop = bw$end[nrow(bw)] # end coordinate of the gene
  cov = rep(0,stop-start+1) # vector with 0s of length the gene

  for(i in 1:nrow(bw)){ # for every record in bw
    cov[(bw$start[i]:bw$end[i])-start+1] = bw$score[i] # assigning score to every position of gene
    # every position of a range is asigned with the same score
  }
  list(cov=cov,x=start:stop) # coverage, start and end gene coordinates on a chromosome
}

# sid - sample id
# jxn - rse object
# gene.grange - granges of genes of interest
import = function(sample.id, gene.grange){
  folder.path = '/home/an/Manananggal/Input/bigWig/'
  sample.path = paste0(folder.path, sample.id, '.bw')
  # download and subset bw file
  bw = rtracklayer::import.bw(sample.path, which=gene.grange) 
  bw
}

get.counts = function(sample.id, rse.gene.tissue){
  # rse for one sample
  rse.tissue.samples = rse.gene.tissue[, sample.id]
  all.genejxn.info = cbind(as.data.frame(rse.gene.tissue@rowRanges)[,c('start','end','strand')],
                      counts = rse.gene.tissue@assays@data$counts[,sample.id]) # ?? зачем еще раз выбирать sample.id?
  all.genejxn.info
}

getRecountCov = function(sample.id, rse.gene.tissue, gene.grange){
  bw = import(sample.id, gene.grange)
  # coverage on a gene, start:stop 
  cov.list.sample = bigWig2Cov(bw) 
  cov.list.sample$juncs =  get.counts(sample.id, rse.gene.tissue)
  # coverages and counts for the entire gene (all jxns)
  cov.list.sample
}

sumCovs = function(gene.cov.samples.list){
  # launched for each of the conditions (2ce)
  # launched for every tissue where significant junction was found
  # creating a template of merged object
  cov.merged = gene.cov.samples.list[[1]] # read coverages for the first sample
  
  # Extract "val" elements from each sublist
  cov.list <- lapply(gene.cov.samples.list, `[[`, "cov")
  # Sum corresponding elements using Reduce
  cov.merged$cov <- Reduce("+", cov.list)

  juncs.list <- lapply(gene.cov.samples.list, `[[`, "juncs")
  counts.list <- lapply(juncs.list, `[[`, "counts")
  cov.merged$juncs$counts <- Reduce("+", counts.list)
  cov.merged
}


#' Plots read coverage
#'
#' @param r read coverage; output of \code{\link{getReadCoverage}}
#' @param min.junc.cov numeric, plots only junctions (introns) with coverage not less than \code{min.junc.cov}
#' @param min.junc.cov.f numeric, plots only junctions (introns) with coverage not less than \code{min.junc.cov.f} of maximal coverage in the region
#' @param plot.junc.only.within logical, plot only juction with both ends within the region, FALSE plots all junctions with at least one end within region. NA plot all junctions overlapping the region.
#' @param reverse reverse x coordinates
#' @param junc.col colour for junction line. Individual color could be specified for each junction
#' @param junc.lwd line width for jucntion line
#' @param ... other parameters for plot function
#'
#' @export

filter.data = function(condition.cov.list, xlim, min.junc.cov,plot.junc.only.within, min.junc.cov.f){
    print(min.junc.cov.f)
    #print(condition.cov.list)
    # choosing only x inside the range of interest
    print(xlim)
    x.in.range.tf = condition.cov.list$x >= xlim[1]-500 &
                    condition.cov.list$x <= xlim[2]+500
    #??? x is from 1 to what?
    
    condition.cov.list$x = condition.cov.list$x[x.in.range.tf]
    condition.cov.list$cov = condition.cov.list$cov[x.in.range.tf]
    
    condition.cov.list$juncs =
      condition.cov.list$juncs[(condition.cov.list$juncs$start == xlim[1] |
                               condition.cov.list$juncs$end == xlim[2]) &
                               condition.cov.list$juncs$counts > min.junc.cov.f,]
    
    condition.cov.list$cov[c(1,length(condition.cov.list$cov))] = 0 # assigning cov 0 to first and last elements of cov
    condition.cov.list
  }


plotReadCov = function(condition.cov.list,
                       min.junc.cov=0,
                       min.junc.cov.f=0,
                       plot.junc.only.within=FALSE,
                       xlim,
                       reverse=FALSE,
                   #    junc.col='blue',
                       junc.lwd=3,
                       bottom.mar=0,...){
  
  condition.cov.list = filter.data(condition.cov.list, xlim, min.junc.cov,plot.junc.only.within, min.junc.cov.f)
  print('filtered')
  print(condition.cov.list$juncs)
  
  # create a graph
  plot(condition.cov.list$x,
       condition.cov.list$cov,
       t='n',
       xlim=c(xlim[1]-500, xlim[2]+500),
       ...)

  # plot verticale lines
  polygon(condition.cov.list$x,
          condition.cov.list$cov,
          col = 'gray',
          border=NA)
  
  if(nrow(condition.cov.list$juncs)>0) {
    for(i in 1:nrow(condition.cov.list$juncs)){
        print(c('i', i,
        condition.cov.list$juncs$start[i],
        condition.cov.list$juncs$end[i],
        condition.cov.list$juncs$cols[i]
       ))
      
      plotArc(condition.cov.list$juncs$start[i],
              condition.cov.list$juncs$end[i],
              condition.cov.list$juncs$counts[i],
              col=condition.cov.list$juncs$cols[i],
              lwd=junc.lwd)
    }
  }
}

#' Plots parabolic arc
#'
#' @param from,to x coordinates of arc
#' @param top highest point of arc
#' @param n number of points
#' @param y.base bottom coordinate of arc
#' @param ... other parameters of lines functoin
plotArc = function(from,to,top,n=100,y.base=0,...){
  len = to - from
  x = seq(from=0,to=len,length.out = n)
  y = x*4*top/len - x^2*(4*top/len^2) 
  # This equation represents a downward-facing parabola that starts at 0, reaches its peak at top, and ends at 0 again.
  lines(x+from,y+y.base,...)
}