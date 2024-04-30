####################################################################
# All helper functions in the utils/plotseg-utils folder are from the copynumber package. In particular, the plotGenome function has been slightly modified to allow CNA state specific coloring of segments.
# Reference: 
## Publication: Nilsen and Liestøl et al. (2012), BMC Genomics
## Package: 10.18129/B9.bioc.copynumber  
## License: Artistic 2.0
####################################################################


plotseglabel <- function(filepath=NULL, filename=NULL, data,assembly='hg38',no_label=F,ylim=NULL,...){
  colnames(data)[1] <- 'sampleID'
  data <- cbind(data[1:2], arm=rep('.',nrow(data)), data[3:ncol(data)])
  if (no_label){
    col = '#d69f7e'
  } else{
    col <- rep(NA, length(data$label))
    col[data$label == '0'] <- '#90be6d'
    col[data$label == '+1'] <- '#f8961e'
    col[data$label == '+2'] <- 'red'
    col[data$label == '-1'] <- '#8ecae6'
    col[data$label == '-2'] <- '#014f86'
  }

  
  if (!is.null(filepath) & !is.null(filename)){
    pdf(paste0(filepath,'/',filename), width = 10, height = 5)
  } else if (!is.null(filepath) & is.null(filename)){
    pdf(paste0(filepath,data[,1][1],".pdf"), width = 10, height = 5)
  }
  plotGenome(segments = data,connect=FALSE,seg.col=col,assembly = assembly,ylim=ylim,...)
  if (!is.null(filepath)){
    dev.off()
  }
}




