quiet <- function(...){
  nullfile <- file("/dev/null", "w")
  sink(nullfile,type = "output")
  out <- eval(...)
  # return output to the console and free up the sink stack
  sink(type = "output")
  # remove the connection object creacted by eval()
  rval <- NULL
  close(nullfile)
  out
}

checkfile <- function(filepath,filetype,pattern){
  filepath <- list.files(filepath,full.names=T, pattern=pattern)
  if(length(filepath)==0) stop(filetype, " is not found")
  if (length(filepath) > 1) stop("more than one ",filetype)
  return(filepath)
}

get.labelseg <- function(output_dir,filename,sample=NULL){
  files <- list.files(output_dir,recursive = T)
  files <- files[grep(filename,files)]
  if (!is.null(sample)){
    files <- as.character(unlist(sapply(sample,FUN = function(x){files[grep(x,files)]})))
  }
  data <- list()
  for ( i in seq_len(length(files)) ){
    result <- read.table(file.path(output_dir , files[i]) , sep = '\t',header=T,colClasses = c("label"="character"))
    data[[i]] <- result
  }
  data <- do.call(rbind,data)
  return(data)
}

seg.cov <- function(seg){
  if (is.null(seg)){
    lowCov <- 0
    highCov <- 0
    normalCov <- 0
    lowCov_ratio <- 0
  } else{
    segNum <- dim(seg)[1]
    segLen <- seg[,4]-seg[,3]
    total.call.length <- sum(segLen)
    lowCov <- round(sum(segLen[seg$label %in% c('+1','-1')])/total.call.length,2)
    highCov <- round(sum(segLen[seg$label %in% c('+2','-2')])/total.call.length,2)
    normalCov <- round(sum(segLen[seg$label == '0'])/total.call.length,2)
    lowdupCov <- round(sum(segLen[seg$label %in% c('+1')])/total.call.length,2)
    lowdelCov <- round(sum(segLen[seg$label %in% c('-1')])/total.call.length,2)
    lowCov_ratio <- lowdupCov/lowdelCov 
    lowCov_ratio[is.na(lowCov_ratio)] <- 1
  }
  score <- data.frame('normal_cov'=normalCov,'lowCNA_cov'=lowCov,'lowCNA_ratio'=lowCov_ratio,'highCNA_cov'=highCov)
  return(score)
}

update.seg <- function(output_dir,sample,oriseg, relabelSeg,report,score,problem,genome='hg38',mergeseg=F, outputfilename, outputplotname){
  mergseg_path <- file.path(output_dir,sample,'merged_segments.seg.txt')
  if (is.null(relabelSeg)){
    if (file.exists(mergseg_path)){
      relabelSeg <- read.table(mergseg_path,sep = '\t',header = T)
    } else{
      relabelSeg <- oriseg
    }
    plotseglabel(filepath = file.path(output_dir,sample),filename=outputplotname,data = relabelSeg,assembly = genome,no_label = T)
    problem <- paste0(problem,';failed-to-label')
  } else{
    experiment = unique(relabelSeg[,1])
    write.table(relabelSeg, file=file.path(output_dir,sample, outputfilename),sep = '\t',quote=F,row.names = F)
    plotseglabel(filepath = file.path(output_dir,sample),filename=outputplotname,data = relabelSeg,assembly = genome)
    if (file.exists(mergseg_path)) system(sprintf("rm  %s", mergseg_path))
  }
  #update report
  idx <- which(report$sample_name == sample)
  report[idx,c('normal_cov','lowCNA_cov','lowCNA_ratio','highCNA_cov')] <- c(score$normal_cov,score$lowCNA_cov,score$lowCNA_ratio,score$highCNA_cov)
  report$note[idx] <- paste0(report$note[idx],";", problem)
  
  if (mergeseg){
    lrrsegsd <- round(sd(rep(relabelSeg[,6],abs(relabelSeg[,5]))),3)
    report[idx,c('segment_num',	'LLR_segsd')] <- c(dim(relabelSeg)[1],lrrsegsd)
  }
  
  return(report)
}

